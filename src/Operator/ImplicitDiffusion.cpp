#include "ImplicitDiffusion.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include "BC/Constant.H"

namespace Operator
{
void
ImplicitDiffusion::Parse(ImplicitDiffusion& value, IO::ParmParse& pp)
{
    pp.query_default("tol_rel", value.tolerance_relative, 1.0e-10);
    pp.query_default("tol_abs", value.tolerance_absolute, 0.0);
    pp.query_default("verbose", value.verbose, 0);
    pp.query_default("max_order", value.max_order, 2);
}

void
ImplicitDiffusion::SetLayout(
    const amrex::Vector<amrex::Geometry>& a_geometry,
    const amrex::Vector<amrex::IntVect>& a_refinement_ratio,
    const Set::Field<Set::Scalar>& layout,
    int number_of_levels, int number_of_components)
{
    Util::Assert(INFO, TEST(number_of_levels > 0));
    Util::Assert(INFO, TEST(number_of_components > 0));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(a_geometry.size())));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(layout.size())));
    Util::Assert(INFO, TEST(number_of_levels == 1 ||
                            number_of_levels - 1 <=
                                static_cast<int>(a_refinement_ratio.size())));

    bool layout_changed = nlevels != number_of_levels ||
                            ncomponents != number_of_components;
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
    ncomponents = number_of_components;
    geometry.assign(a_geometry.begin(), a_geometry.begin() + nlevels);
    refinement_ratio.assign(
        a_refinement_ratio.begin(), a_refinement_ratio.begin() + nlevels - 1);
    if (!layout_changed) return;

    grids.resize(nlevels);
    distribution_mapping.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        grids[lev] = layout[lev]->boxArray();
        distribution_mapping[lev] = layout[lev]->DistributionMap();
    }

    solution.Define(nlevels, grids, distribution_mapping, ncomponents, 1);
    rhs.Define(nlevels, grids, distribution_mapping, ncomponents, 0);
    a_coefficient.Define(nlevels, grids, distribution_mapping, 1, 0);
    b_coefficient.Define(nlevels, grids, distribution_mapping, ncomponents, 1);
    face_b_coefficient.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        a_coefficient[lev]->setVal(0.0);
        b_coefficient[lev]->setVal(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_grids = grids[lev];
            face_grids.surroundingNodes(d);
            face_b_coefficient[lev][d].define(
                face_grids, distribution_mapping[lev], ncomponents, 0);
        }
    }
}

void
ImplicitDiffusion::Solve(Set::Scalar time, Set::Scalar dt,
                        const amrex::BCRec& bc, bool harmonic_averaging)
{
    Solve(time, dt, amrex::Vector<amrex::BCRec>(ncomponents, bc),
            harmonic_averaging);
}

void
ImplicitDiffusion::Solve(Set::Scalar time, Set::Scalar dt,
                        const amrex::Vector<amrex::BCRec>& bc,
                        bool harmonic_averaging)
{
    BL_PROFILE("Operator::ImplicitDiffusion::Solve");
    Util::Assert(INFO, TEST(nlevels > 0));
    Util::Assert(INFO, TEST(dt > 0.0));
    Util::Assert(INFO, TEST(static_cast<int>(bc.size()) == ncomponents));

    BC::Constant::ZeroNeumann coefficient_bc(ncomponents);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        coefficient_bc.define(geometry[lev]);
        if (lev == 0)
        {
            coefficient_bc.FillBoundary(
                *b_coefficient[lev], 0, ncomponents, time, 0);
            b_coefficient[lev]->FillBoundary(geometry[lev].periodicity());
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> coarse{b_coefficient[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fine{b_coefficient[lev].get()};
            amrex::Vector<amrex::Real> coarse_time{time};
            amrex::Vector<amrex::Real> fine_time{time};
            amrex::Vector<amrex::BCRec> bcs(
                ncomponents, coefficient_bc.GetBCRec());
            amrex::FillPatchTwoLevels(
                *b_coefficient[lev], time, coarse, coarse_time, fine, fine_time,
                0, 0, ncomponents, geometry[lev - 1], geometry[lev],
                coefficient_bc, 0, coefficient_bc, 0,
                refinement_ratio[lev - 1], &amrex::cell_cons_interp, bcs, 0);
        }
    }

    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>>
        face_b_coefficient_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> face_ptr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            face_ptr[d] = &face_b_coefficient[lev][d];
            face_b_coefficient_ptr[lev][d] = &face_b_coefficient[lev][d];
        }
        amrex::average_cellcenter_to_face(
            face_ptr, *b_coefficient[lev], geometry[lev], ncomponents,
            harmonic_averaging, 0);
    }

    amrex::LPInfo info;
    amrex::MLABecLaplacian diffusion(
        geometry, grids, distribution_mapping, info, {}, ncomponents);
    diffusion.setMaxOrder(max_order);
    amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
        boundary_lo(ncomponents);
    amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
        boundary_hi(ncomponents);
    for (int n = 0; n < ncomponents; ++n)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            if (geometry[0].isPeriodic(d))
            {
                boundary_lo[n][d] = amrex::LinOpBCType::Periodic;
                boundary_hi[n][d] = amrex::LinOpBCType::Periodic;
            }
            else
            {
                boundary_lo[n][d] = BC::BCUtil::IsDirichlet(bc[n].lo(d)) ?
                    amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
                boundary_hi[n][d] = BC::BCUtil::IsDirichlet(bc[n].hi(d)) ?
                    amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
            }
        }
    }
    diffusion.setDomainBC(boundary_lo, boundary_hi);
    diffusion.setScalars(1.0, dt);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        diffusion.setLevelBC(lev, solution[lev].get());
        diffusion.setACoeffs(lev, *a_coefficient[lev]);
        diffusion.setBCoeffs(lev, face_b_coefficient_ptr[lev]);
    }

    amrex::Vector<amrex::MultiFab*> solution_ptr(nlevels);
    amrex::Vector<amrex::MultiFab const*> rhs_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution_ptr[lev] = solution[lev].get();
        rhs_ptr[lev] = rhs[lev].get();
    }

    amrex::MLMG solver(diffusion);
    solver.setVerbose(verbose);
    solver.setFinalFillBC(true);
    solver.solve(solution_ptr, rhs_ptr,
                tolerance_relative, tolerance_absolute);
}
}
