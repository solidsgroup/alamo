#ifndef OPERATOR_PRESSUREPOISSON_H
#define OPERATOR_PRESSUREPOISSON_H

#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

#include "IO/ParmParse.H"
#include "Set/Set.H"

namespace Operator
{
/// Composite AMR variable-coefficient Poisson solve used by pressure projections.
///
/// The caller supplies the cell-centered right-hand side and coefficient. This
/// class owns the solve scratch data and handles coarse/fine coefficient fills,
/// cell-to-face conversion, and the AMReX MLMG setup.
class PressurePoisson
{
public:
    PressurePoisson() = default;
    static void Parse(PressurePoisson& value, IO::ParmParse& pp);

    void SetLayout(const amrex::Vector<amrex::Geometry>& geometry,
                    const amrex::Vector<amrex::IntVect>& refinement_ratio,
                    const Set::Field<Set::Scalar>& layout,
                    int number_of_levels);

    amrex::MultiFab& RHS(int lev) { return *rhs[lev]; }
    amrex::MultiFab& Coefficient(int lev) { return *coefficient[lev]; }
    const amrex::MultiFab& Solution(int lev) const { return *solution[lev]; }

    void PrepareRHS(int lev, const amrex::MultiFab& velocity, Set::Scalar dt);
    void Solve(Set::Scalar time, const amrex::BCRec& pressure_bc);
    void ApplyCorrection(int lev, amrex::MultiFab& velocity, Set::Scalar dt);

private:
    int nlevels = 0;
    Set::Scalar tolerance_relative = 1.0e-11;
    Set::Scalar tolerance_absolute = 1.0e-12;
    int verbose = 0;
    int max_order = 2;
    int max_coarsening_level = -1;
    bool abort_on_fail = true;
    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratio;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> distribution_mapping;
    Set::Field<Set::Scalar> solution;
    Set::Field<Set::Scalar> rhs;
    Set::Field<Set::Scalar> coefficient;
    Set::Field<Set::Scalar> divergence;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> face_coefficient;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> face_velocity;
};
}

#endif
