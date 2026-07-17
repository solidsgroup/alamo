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
    for (int lev = 0; lev < nlevels; ++lev)
    {
        solution[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        coefficient[lev]->setVal(0.0);
    }
}

void
PressurePoisson::Solve(Set::Scalar time)
{
    BL_PROFILE("Operator::PressurePoisson::Solve");

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
                coefficient_bc, 0, coefficient_bc, 0, refinement_ratio[lev - 1],
                &amrex::cell_cons_interp, bcs, 0);
        }
    }

    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> face_coefficient(nlevels);
    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>> face_coefficient_ptr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> face_ptr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_grids = grids[lev];
            face_grids.surroundingNodes(d);
            face_coefficient[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 0);
            face_ptr[d] = &face_coefficient[lev][d];
            face_coefficient_ptr[lev][d] = &face_coefficient[lev][d];
        }
        amrex::average_cellcenter_to_face(
            face_ptr, *coefficient[lev], geometry[lev], 1, true, 0);
    }

    amrex::LPInfo info;
    amrex::MLABecLaplacian poisson(geometry, grids, distribution_mapping, info);
    poisson.setMaxOrder(max_order);
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> boundary_lo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> boundary_hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        boundary_lo[d] = geometry[0].isPeriodic(d) ?
            amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
        boundary_hi[d] = geometry[0].isPeriodic(d) ?
            amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
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
}
