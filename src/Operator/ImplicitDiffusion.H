#ifndef OPERATOR_IMPLICITDIFFUSION_H
#define OPERATOR_IMPLICITDIFFUSION_H

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "IO/ParmParse.H"
#include "Set/Set.H"

namespace Operator
{
/// Composite AMR backward-Euler solve for
/// (a - dt div(b grad)) phi = rhs.
class ImplicitDiffusion
{
public:
    ImplicitDiffusion() = default;
    static void Parse(ImplicitDiffusion& value, IO::ParmParse& pp);

    void SetLayout(const amrex::Vector<amrex::Geometry>& geometry,
                    const amrex::Vector<amrex::IntVect>& refinement_ratio,
                    const Set::Field<Set::Scalar>& layout,
                    int number_of_levels, int number_of_components);

    amrex::MultiFab& Solution(int lev) { return *solution[lev]; }
    amrex::MultiFab& RHS(int lev) { return *rhs[lev]; }
    amrex::MultiFab& ACoefficient(int lev) { return *a_coefficient[lev]; }
    amrex::MultiFab& BCoefficient(int lev) { return *b_coefficient[lev]; }

    void Solve(Set::Scalar time, Set::Scalar dt, const amrex::BCRec& bc,
                bool harmonic_averaging);
    void Solve(Set::Scalar time, Set::Scalar dt,
                const amrex::Vector<amrex::BCRec>& bc,
                bool harmonic_averaging);

private:
    int nlevels = 0;
    int ncomponents = 0;
    Set::Scalar tolerance_relative = 1.0e-10;
    Set::Scalar tolerance_absolute = 0.0;
    int verbose = 0;
    int max_order = 2;
    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratio;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> distribution_mapping;
    Set::Field<Set::Scalar> solution;
    Set::Field<Set::Scalar> rhs;
    Set::Field<Set::Scalar> a_coefficient;
    Set::Field<Set::Scalar> b_coefficient;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
        face_b_coefficient;
};
}

#endif
