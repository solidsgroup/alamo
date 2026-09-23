#include "PressurePoisson.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include "BC/Constant.H"
#include "Unit/Unit.H"

namespace Operator
{
void
PressurePoisson::Parse(PressurePoisson& value, IO::ParmParse& pp)
{
    pp.query_default("tol_rel", value.tolerance_relative, 1.0e-11);
    pp.query_default("tol_abs", value.tolerance_absolute,
                     "1.0e-12_1/s", 1.0 / Unit::Time());
    pp.query_default("verbose", value.verbose, 0);
    pp.query_default("max_order", value.max_order, 2);
    if (!(value.tolerance_relative > 0.0) ||
        !(value.tolerance_absolute >= 0.0))
        Util::Exception(INFO,
            "projection.tol_rel must be positive and projection.tol_abs "
            "must be nonnegative");
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

    const bool had_mixed_velocity_state = mixed_velocity_state_initialized;
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
    occupancy.Define(nlevels, grids, distribution_mapping, 1, 1);
    divergence.Define(nlevels, grids, distribution_mapping, 1, 0);
    cell_velocity_predictor.Define(
        nlevels, grids, distribution_mapping, AMREX_SPACEDIM, 1);
    cell_velocity_reference.Define(
        nlevels, grids, distribution_mapping, AMREX_SPACEDIM, 1);
    face_coefficient.resize(nlevels);
    face_velocity.resize(nlevels);
    face_velocity_base.resize(nlevels);
    face_momentum_exchange.resize(nlevels);
    mixed_velocity_state_initialized = had_mixed_velocity_state;
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        coefficient[lev]->setVal(0.0);
        occupancy[lev]->setVal(1.0);
        divergence[lev]->setVal(0.0);
        cell_velocity_predictor[lev]->setVal(0.0);
        cell_velocity_reference[lev]->setVal(0.0);
        if (had_mixed_velocity_state)
            amrex::MultiFab::Copy(*cell_velocity_reference[lev], *layout[lev],
                0, 0, AMREX_SPACEDIM, 1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_grids = grids[lev];
            face_grids.surroundingNodes(d);
            face_coefficient[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 0);
            face_velocity[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 1);
            face_velocity_base[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 1);
            face_momentum_exchange[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 0);
            face_velocity[lev][d].setVal(0.0);
            face_velocity_base[lev][d].setVal(0.0);
            face_momentum_exchange[lev][d].setVal(0.0);
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
    int lev, const amrex::MultiFab& cell_velocity, Set::Scalar dt,
    const FaceField* face_acceleration)
{
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> face_ptr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        face_ptr[d] = &face_velocity[lev][d];
        const int di = d == 0;
        const int dj = d == 1;
        const int dk = d == 2;
        for (amrex::MFIter mfi(face_velocity[lev][d],
                amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto u = cell_velocity.const_array(mfi);
            const auto face = face_velocity[lev][d].array(mfi);
            const auto base = face_velocity_base[lev][d].array(mfi);
            const bool volume = volume_velocity;
            const auto weight = occupancy[lev]->const_array(mfi);
            amrex::Array4<const Set::Scalar> acceleration;
            const bool has_face_acceleration = face_acceleration != nullptr;
            if (has_face_acceleration)
                acceleration = (*face_acceleration)[d]->const_array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    base(i,j,k) = 0.5 *
                        (u(i-di,j-dj,k-dk,d) + u(i,j,k,d));
                    if (volume) base(i,j,k) = 0.5 *
                        (weight(i-di,j-dj,k-dk) * u(i-di,j-dj,k-dk,d) +
                         weight(i,j,k) * u(i,j,k,d));
                    face(i,j,k) = base(i,j,k);
                    if (has_face_acceleration)
                        face(i,j,k) += dt * acceleration(i,j,k);
                });
        }
        // A periodic seam is one geometric face.  Reconcile overlapping valid
        // copies before taking the divergence; this is synchronization, not a
        // mean-flow or pressure-gradient correction.
        face_velocity[lev][d].FillBoundaryAndSync(
            geometry[lev].periodicity());
        face_velocity_base[lev][d].FillBoundaryAndSync(
            geometry[lev].periodicity());
    }
    amrex::computeDivergence(*divergence[lev], face_ptr, geometry[lev]);

    const Set::Scalar inverse_dt = 1.0 / dt;
    for (amrex::MFIter mfi(*rhs[lev], amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        const auto source = rhs[lev]->const_array(mfi);
        const auto div = divergence[lev]->const_array(mfi);
        const auto projection_rhs = rhs[lev]->array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                projection_rhs(i,j,k) =
                    (div(i,j,k) - source(i,j,k)) * inverse_dt;
            });
    }
}

void
PressurePoisson::Solve(Set::Scalar /*time*/, const amrex::BCRec& pressure_bc)
{
    BL_PROFILE("Operator::PressurePoisson::Solve");

    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>>
        face_coefficient_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            face_coefficient_ptr[lev][d] = &face_coefficient[lev][d];

    amrex::LPInfo info;
    amrex::MLABecLaplacian poisson(
        geometry, grids, distribution_mapping, info);
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

    // A single pressure outlet carries the net volume production. Seed the
    // solve with a quadratic pressure lifting that supplies that
    // boundary flux. Its value vanishes at the outlet and its normal gradient
    // vanishes at the opposite wall. MLMG still solves the original RHS; no
    // source mean is discarded and no extra momentum flux is introduced.
    int outlets = 0, outlet_direction = 0;
    bool outlet_high = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        if (!geometry[0].isPeriodic(d))
            for (int side = 0; side < 2; ++side)
                if (BC::BCUtil::IsDirichlet(side ? pressure_bc.hi(d) : pressure_bc.lo(d)))
                {
                    ++outlets;
                    outlet_direction = d;
                    outlet_high = side;
                }
    Set::Scalar curvature = 0.0;
    if (outlets == 1)
    {
        for (int lev = nlevels - 1; lev > 0; --lev)
            amrex::average_down(*rhs[lev], *rhs[lev-1],
                geometry[lev], geometry[lev-1], 0, 1, refinement_ratio[lev-1]);
        const int d = outlet_direction;
        const auto domain = geometry[0].Domain();
        const int boundary = outlet_high ? domain.bigEnd(d) + 1 : domain.smallEnd(d);
        Set::Scalar boundary_beta = 0.0;
        // Iterate cells, so overlapping face boxes cannot count an outlet
        // face twice. The face coefficient uses the same layout and index.
        for (amrex::MFIter mfi(*rhs[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto beta = face_coefficient[0][d].const_array(mfi);
            amrex::ReduceOps<amrex::ReduceOpSum> op;
            amrex::ReduceData<Set::Scalar> data(op);
            using Tuple = typename decltype(data)::Type;
            const bool high = outlet_high;
            op.eval(mfi.tilebox(), data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> Tuple
                {
                    const int index = (d == 0 ? i : (d == 1 ? j : k)) + high;
                    return {index == boundary ?
                        beta(i + (high && d == 0), j + (high && d == 1),
                             k + (high && d == 2)) : 0.0};
                });
            boundary_beta += amrex::get<0>(data.value());
        }
        amrex::ParallelDescriptor::ReduceRealSum(boundary_beta);
        if (boundary_beta > 0.0)
            curvature = rhs[0]->sum(0) / (domain.length(d) * boundary_beta);
    }

    amrex::Vector<amrex::MultiFab*> solution_ptr(nlevels);
    amrex::Vector<amrex::MultiFab const*> rhs_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        if (outlets == 1)
        {
            const int d = outlet_direction;
            const auto domain = geometry[lev].Domain();
            const Set::Scalar spacing = geometry[lev].CellSize(d);
            const Set::Scalar length = geometry[lev].ProbLength(d);
            const bool high = outlet_high;
            for (amrex::MFIter mfi(*solution[lev], amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const auto phi = solution[lev]->array(mfi);
                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        const int index = d == 0 ? i : (d == 1 ? j : k);
                        const Set::Scalar x = spacing * (high ?
                            index - domain.smallEnd(d) + 0.5 :
                            domain.bigEnd(d) - index + 0.5);
                        // The cell-centered Dirichlet face gradient uses an
                        // odd ghost value; this offset gives the exact net
                        // flux with the second-order boundary stencil.
                        phi(i,j,k) = 0.5 * curvature *
                            (x*x - length*length - 0.25*spacing*spacing);
                    });
            }
        }
        solution_ptr[lev] = solution[lev].get();
        rhs_ptr[lev] = rhs[lev].get();
    }

    amrex::MLMG solver(poisson);
    // Thin periodic strips can leave an elongated bottom grid. The pressure
    // operator is symmetric; CG avoids BiCGStab breakdown on nearly planar
    // residuals, while allowing the bottom solve to resolve the long axis.
    solver.setBottomSolver(amrex::MLMG::BottomSolver::cg);
    solver.setBottomMaxIter(2000);
    solver.setVerbose(verbose);
    solver.setFinalFillBC(true);
    solver.solve(solution_ptr, rhs_ptr,
                tolerance_relative, tolerance_absolute);
}

void
PressurePoisson::ApplyCorrection(
    int lev, amrex::MultiFab& cell_velocity, Set::Scalar dt)
{
    const auto dx = geometry[lev].CellSizeArray();
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        const int di = d == 0;
        const int dj = d == 1;
        const int dk = d == 2;
        for (amrex::MFIter mfi(face_velocity[lev][d],
                amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto phi = solution[lev]->const_array(mfi);
            const auto beta = face_coefficient[lev][d].const_array(mfi);
            const auto face = face_velocity[lev][d].array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    face(i,j,k) -= dt * beta(i,j,k) *
                        (phi(i,j,k) - phi(i-di,j-dj,k-dk)) / dx[d];
                });
        }
        face_velocity[lev][d].FillBoundaryAndSync(
            geometry[lev].periodicity());
    }

    amrex::MultiFab& predictor = *cell_velocity_predictor[lev];
    amrex::MultiFab::Copy(predictor, cell_velocity, 0, 0,
        AMREX_SPACEDIM, 1);
    for (amrex::MFIter mfi(cell_velocity, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        const auto u = cell_velocity.array(mfi);
        const auto u_predictor = predictor.const_array(mfi);
        amrex::GpuArray<amrex::Array4<const Set::Scalar>, AMREX_SPACEDIM>
            projected_face;
        amrex::GpuArray<amrex::Array4<const Set::Scalar>, AMREX_SPACEDIM>
            base_face;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            projected_face[d] = face_velocity[lev][d].const_array(mfi);
            base_face[d] = face_velocity_base[lev][d].const_array(mfi);
        }
        const auto weight = occupancy[lev]->const_array(mfi);
        const bool weighted = volume_velocity;
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    const int di = d == 0;
                    const int dj = d == 1;
                    const int dk = d == 2;
                    u(i,j,k,d) = u_predictor(i,j,k,d) +
                        0.5 *
                            ((projected_face[d](i,j,k) -
                              base_face[d](i,j,k)) +
                             (projected_face[d](i+di,j+dj,k+dk) -
                              base_face[d](i+di,j+dj,k+dk))) /
                            (weighted ? weight(i,j,k) : 1.0);
                }
            });
    }
}

void
PressurePoisson::ReconcileCellVelocity(
    Set::Field<Set::Scalar>& velocity,
    const Set::Field<Set::Scalar>& density,
    bool projection_increment,
    const CompositeFaceField* face_weight)
{
    // Arithmetic cell-to-face interpolation has an exact alternating-cell
    // nullspace.  Exchange only the newly accumulated normal momentum between
    // adjacent cells.  On a uniform field this is the compact [1/4,1/2,1/4]
    // compatible reconstruction; writing it as a shared face flux preserves
    // linear momentum, angular momentum, periodic seams, and coarse/fine
    // conservation without damping the previously projected velocity state.
    // An optional composite face weight localizes that same conservative
    // exchange without changing its equal-and-opposite momentum transfer.
    for (int lev = 0; lev < nlevels; ++lev)
    {
        const auto dx = geometry[lev].CellSizeArray();
        const amrex::Box domain = geometry[lev].Domain();
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            const int di = d == 0;
            const int dj = d == 1;
            const int dk = d == 2;
            const int domain_face_lo = domain.smallEnd(d);
            const int domain_face_hi = domain.bigEnd(d) + 1;
            const bool periodic = geometry[lev].isPeriodic(d);
            const bool weighted = face_weight != nullptr;
            amrex::MultiFab& exchange = face_momentum_exchange[lev][d];
            for (amrex::MFIter mfi(exchange, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                const auto cell_velocity = velocity[lev]->const_array(mfi);
                const auto reference = projection_increment ?
                    cell_velocity_predictor[lev]->const_array(mfi) :
                    cell_velocity_reference[lev]->const_array(mfi);
                const auto rho = density[lev]->const_array(mfi);
                const auto flux = exchange.array(mfi);
                amrex::Array4<const Set::Scalar> weight;
                if (weighted)
                    weight = (*face_weight)[lev][d]->const_array(mfi);
                const bool has_mixed_velocity_state =
                    projection_increment || mixed_velocity_state_initialized;
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        const int face_index = d == 0 ? i :
                            (d == 1 ? j : k);
                        if (!periodic &&
                            (face_index == domain_face_lo ||
                             face_index == domain_face_hi))
                        {
                            flux(i,j,k) = 0.0;
                            return;
                        }
                        const Set::Scalar increment_lo =
                            cell_velocity(i-di,j-dj,k-dk,d) -
                            (has_mixed_velocity_state ?
                                reference(i-di,j-dj,k-dk,d) : 0.0);
                        const Set::Scalar increment_hi =
                            cell_velocity(i,j,k,d) -
                            (has_mixed_velocity_state ?
                                reference(i,j,k,d) : 0.0);
                        const Set::Scalar rho_lo =
                            rho(i-di,j-dj,k-dk);
                        const Set::Scalar rho_hi = rho(i,j,k);
                        const Set::Scalar face_density =
                            2.0 * rho_lo * rho_hi / (rho_lo + rho_hi);
                        const Set::Scalar localization = weighted ?
                            weight(i,j,k) : 1.0;
                        flux(i,j,k) = 0.25 * localization * dx[d] *
                            face_density *
                            (increment_hi - increment_lo);
                    });
            }
            exchange.FillBoundaryAndSync(geometry[lev].periodicity());
        }
    }

    for (int lev = nlevels - 1; lev > 0; --lev)
    {
        amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> fine;
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> coarse;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            fine[d] = &face_momentum_exchange[lev][d];
            coarse[d] = &face_momentum_exchange[lev-1][d];
        }
        amrex::average_down_faces(
            fine, coarse, refinement_ratio[lev-1], geometry[lev-1]);
    }

    for (int lev = 0; lev < nlevels; ++lev)
    {
        const auto dx = geometry[lev].CellSizeArray();
        for (amrex::MFIter mfi(*velocity[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto u = velocity[lev]->array(mfi);
            const auto rho = density[lev]->const_array(mfi);
            amrex::GpuArray<amrex::Array4<const Set::Scalar>, AMREX_SPACEDIM>
                exchange;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                exchange[d] =
                    face_momentum_exchange[lev][d].const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        const int di = d == 0;
                        const int dj = d == 1;
                        const int dk = d == 2;
                        u(i,j,k,d) +=
                            (exchange[d](i+di,j+dj,k+dk) -
                             exchange[d](i,j,k)) /
                            (dx[d] * rho(i,j,k));
                    }
                });
        }
    }
}

void
PressurePoisson::CommitVelocityState(
    const Set::Field<Set::Scalar>& velocity)
{
    for (int lev = 0; lev < nlevels; ++lev)
        amrex::MultiFab::Copy(*cell_velocity_reference[lev], *velocity[lev],
            0, 0, AMREX_SPACEDIM, 1);
    mixed_velocity_state_initialized = true;
}
}
