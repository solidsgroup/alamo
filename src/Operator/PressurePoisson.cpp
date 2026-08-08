#include "PressurePoisson.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include "BC/Constant.H"

namespace Operator
{
void
PressurePoisson::Parse(PressurePoisson& value, IO::ParmParse& pp)
{
    pp.query_default("tol_rel", value.tolerance_relative, 1.0e-11);
    pp.query_default("tol_abs", value.tolerance_absolute, 1.0e-12);
    pp.query_default("verbose", value.verbose, 0);
    pp.query_default("max_order", value.max_order, 2);
}

void
PressurePoisson::SetLayout(
    const amrex::Vector<amrex::Geometry>& a_geometry,
    const amrex::Vector<amrex::IntVect>& a_refinement_ratio,
    const Set::Field<Set::Scalar>& layout,
    int number_of_levels)
{
    Util::Assert(INFO, TEST(number_of_levels > 0));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(a_geometry.size())));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(layout.size())));
    Util::Assert(INFO, TEST(number_of_levels == 1 ||
                            number_of_levels - 1 <= static_cast<int>(a_refinement_ratio.size())));

    bool layout_changed = nlevels != number_of_levels;
    if (!layout_changed)
    {
        for (int lev = 0; lev < number_of_levels; ++lev)
        {
            if (solution[lev]->boxArray() != layout[lev]->boxArray() ||
                solution[lev]->DistributionMap() != layout[lev]->DistributionMap())
            {
                layout_changed = true;
                break;
            }
        }
    }

    nlevels = number_of_levels;
    geometry.assign(a_geometry.begin(), a_geometry.begin() + nlevels);
    refinement_ratio.assign(a_refinement_ratio.begin(), a_refinement_ratio.begin() + nlevels - 1);
    if (!layout_changed) return;

    grids.resize(nlevels);
    distribution_mapping.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        grids[lev] = layout[lev]->boxArray();
        distribution_mapping[lev] = layout[lev]->DistributionMap();
    }

    solution.Define(nlevels, grids, distribution_mapping, 1, 1);
    rhs.Define(nlevels, grids, distribution_mapping, 1, 0);
    coefficient.Define(nlevels, grids, distribution_mapping, 1, 1);
    divergence.Define(nlevels, grids, distribution_mapping, 1, 0);
    cell_velocity_predictor.Define(
        nlevels, grids, distribution_mapping, AMREX_SPACEDIM, 1);
    face_coefficient.resize(nlevels);
    face_velocity.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        coefficient[lev]->setVal(0.0);
        divergence[lev]->setVal(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_grids = grids[lev];
            face_grids.surroundingNodes(d);
            face_coefficient[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 0);
            face_velocity[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 1);
        }
    }
}

void
PressurePoisson::PrepareCoefficients(Set::Scalar time)
{
    BC::Constant::ZeroNeumann coefficient_bc(1);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        coefficient_bc.define(geometry[lev]);
        if (lev == 0)
        {
            coefficient_bc.FillBoundary(*coefficient[lev], 0, 1, time, 0);
            coefficient[lev]->FillBoundary(geometry[lev].periodicity());
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> coarse{coefficient[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fine{coefficient[lev].get()};
            amrex::Vector<amrex::Real> coarse_time{time};
            amrex::Vector<amrex::Real> fine_time{time};
            amrex::Vector<amrex::BCRec> bcs(1, coefficient_bc.GetBCRec());
            amrex::FillPatchTwoLevels(
                *coefficient[lev], time, coarse, coarse_time, fine, fine_time,
                0, 0, 1, geometry[lev - 1], geometry[lev],
                coefficient_bc, 0, coefficient_bc, 0,
                refinement_ratio[lev - 1], &amrex::cell_cons_interp, bcs, 0);
        }
    }

    for (int lev = 0; lev < nlevels; ++lev)
    {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> face_ptr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            face_ptr[d] = &face_coefficient[lev][d];
        amrex::average_cellcenter_to_face(
            face_ptr, *coefficient[lev], geometry[lev], 1, true, 0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            face_coefficient[lev][d].OverrideSync(
                geometry[lev].periodicity());
    }
}

void
PressurePoisson::PrepareRHS(
    int lev, const amrex::MultiFab& velocity, Set::Scalar dt,
    const FaceField* face_capillary_acceleration)
{
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> face_velocity_const_ptr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        face_velocity_const_ptr[d] = &face_velocity[lev][d];
        const int di = d == 0;
        const int dj = d == 1;
        const int dk = d == 2;
        for (amrex::MFIter mfi(face_velocity[lev][d], amrex::TilingIfNotGPU());
            mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto u = velocity.const_array(mfi);
            const auto face = face_velocity[lev][d].array(mfi);
            amrex::Array4<const Set::Scalar> capillary_acceleration;
            const bool has_capillary_acceleration =
                face_capillary_acceleration != nullptr;
            if (has_capillary_acceleration)
                capillary_acceleration =
                    (*face_capillary_acceleration)[d]->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                face(i,j,k) = 0.5 *
                    (u(i-di,j-dj,k-dk,d) + u(i,j,k,d));
                if (has_capillary_acceleration)
                    face(i,j,k) += dt * capillary_acceleration(i,j,k);
            });
        }
        // Face-centered BoxArrays overlap at patch boundaries and contain
        // two valid representations of a periodic seam.  FillBoundary alone
        // updates ghost cells but does not reconcile those valid nodal
        // values.  Keep one flux on every geometric face before taking its
        // divergence so the periodic cells see equal and opposite fluxes.
        face_velocity[lev][d].FillBoundaryAndSync(
            geometry[lev].periodicity());
    }
    amrex::computeDivergence(
        *divergence[lev], face_velocity_const_ptr, geometry[lev]);

    const Set::Scalar inv_dt = 1.0 / dt;
    for (amrex::MFIter mfi(*rhs[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        const auto source = rhs[lev]->const_array(mfi);
        const auto div = divergence[lev]->const_array(mfi);
        const auto projection_rhs = rhs[lev]->array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            projection_rhs(i,j,k) = (div(i,j,k) - source(i,j,k)) * inv_dt;
        });
    }
}

void
PressurePoisson::Solve(Set::Scalar /*time*/, const amrex::BCRec& pressure_bc)
{
    BL_PROFILE("Operator::PressurePoisson::Solve");

    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>> face_coefficient_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            face_coefficient_ptr[lev][d] = &face_coefficient[lev][d];
    }

    amrex::LPInfo info;
    amrex::MLABecLaplacian poisson(geometry, grids, distribution_mapping, info);
    poisson.setMaxOrder(max_order);
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> boundary_lo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> boundary_hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        if (geometry[0].isPeriodic(d))
        {
            boundary_lo[d] = amrex::LinOpBCType::Periodic;
            boundary_hi[d] = amrex::LinOpBCType::Periodic;
        }
        else
        {
            boundary_lo[d] = BC::BCUtil::IsDirichlet(pressure_bc.lo(d)) ?
                amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
            boundary_hi[d] = BC::BCUtil::IsDirichlet(pressure_bc.hi(d)) ?
                amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
        }
    }
    poisson.setDomainBC(boundary_lo, boundary_hi);
    poisson.setScalars(0.0, -1.0);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        poisson.setLevelBC(lev, nullptr);
        poisson.setACoeffs(lev, 0.0);
        poisson.setBCoeffs(lev, face_coefficient_ptr[lev]);
    }

    amrex::Vector<amrex::MultiFab*> solution_ptr(nlevels);
    amrex::Vector<amrex::MultiFab const*> rhs_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        solution_ptr[lev] = solution[lev].get();
        rhs_ptr[lev] = rhs[lev].get();
    }

    amrex::MLMG solver(poisson);
    solver.setVerbose(verbose);
    solver.setFinalFillBC(true);
    solver.solve(solution_ptr, rhs_ptr,
                tolerance_relative, tolerance_absolute);
}

void
PressurePoisson::ApplyCorrection(
    int lev, amrex::MultiFab& velocity, Set::Scalar dt)
{
    const auto dx = geometry[lev].CellSizeArray();
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        const int di = d == 0;
        const int dj = d == 1;
        const int dk = d == 2;
        for (amrex::MFIter mfi(face_velocity[lev][d], amrex::TilingIfNotGPU());
            mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto phi = solution[lev]->const_array(mfi);
            const auto beta = face_coefficient[lev][d].const_array(mfi);
            const auto face = face_velocity[lev][d].array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                face(i,j,k) -= dt * beta(i,j,k) *
                    (phi(i,j,k) - phi(i-di,j-dj,k-dk)) / dx[d];
            });
        }
        // The pressure correction is also evaluated independently on each
        // overlapping face.  Synchronize the valid copies as well as the
        // ghosts before reconstructing the cell-centered increment.
        face_velocity[lev][d].FillBoundaryAndSync(
            geometry[lev].periodicity());
    }

    // Retain the pre-projection cell field so the pressure-induced face
    // increment can be applied without modifying the predictor itself.
    amrex::MultiFab& predictor = *cell_velocity_predictor[lev];
    amrex::MultiFab::Copy(
        predictor, velocity, 0, 0, AMREX_SPACEDIM, 1);
    for (amrex::MFIter mfi(velocity, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        const auto u = velocity.array(mfi);
        const auto u_predictor = predictor.const_array(mfi);
        amrex::GpuArray<amrex::Array4<const Set::Scalar>, AMREX_SPACEDIM>
            projected_face;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            projected_face[d] = face_velocity[lev][d].const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                const int di = d == 0;
                const int dj = d == 1;
                const int dk = d == 2;
                const Set::Scalar predictor_face_lo = 0.5 *
                    (u_predictor(i-di,j-dj,k-dk,d) +
                     u_predictor(i,j,k,d));
                const Set::Scalar predictor_face_hi = 0.5 *
                    (u_predictor(i,j,k,d) +
                     u_predictor(i+di,j+dj,k+dk,d));
                const Set::Scalar face_increment = 0.5 *
                    ((projected_face[d](i,j,k) - predictor_face_lo) +
                     (projected_face[d](i+di,j+dj,k+dk) -
                      predictor_face_hi));

                u(i,j,k,d) = u_predictor(i,j,k,d) + face_increment;
            }
        });
    }
}
}
