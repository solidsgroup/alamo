#ifndef OPERATOR_DIFFUSION_H
#define OPERATOR_DIFFUSION_H

#include <map>
#include <memory>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "BC/BC.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"

namespace Operator
{
/// Composite AMR solve for variable-coefficient backward-Euler diffusion.
///
/// This class solves
///
///   a (q - q_old) - dt div(B grad(q)) = dt source
///
/// for one or more cell-centered components sharing the coefficients. Scalar
/// mobility uses b I; tensor mobility adds a full symmetric tensor B. This
/// class owns the AMReX solve data and all coarse/fine bookkeeping.
class Diffusion
{
public:
    Diffusion() = default;
    static void Parse(Diffusion& value, IO::ParmParse& pp);

    void SetLayout(const amrex::Vector<amrex::Geometry>& geometry,
                    const amrex::Vector<amrex::IntVect>& refinement_ratio,
                    const Set::Field<Set::Scalar>& layout,
                    int number_of_levels, int number_of_components);

    amrex::MultiFab& State(int lev, int ncomp);
    amrex::MultiFab& Source(int lev, int ncomp);
    amrex::MultiFab& Mass(int lev, int ncomp);
    amrex::MultiFab& Mobility(int lev, int ncomp);
    amrex::MultiFab& TensorMobility(int lev, int ncomp);

    void Solve(Set::Scalar time, Set::Scalar dt,
                const amrex::BCRec& boundary_condition, int ncomp,
                bool tensor_mobility = false,
                bool include_source = false);
    void Solve(Set::Scalar time, Set::Scalar dt,
                const amrex::Vector<amrex::BCRec>& boundary_conditions,
                int ncomp, bool tensor_mobility = false,
                bool include_source = false);
    void FillBoundary(Set::Field<Set::Scalar>& field,
                        BC::BC<Set::Scalar>& boundary_condition,
                        Set::Scalar time, int ncomp) const;
    void Synchronize(Set::Field<Set::Scalar>& field, int ncomp) const;

private:
    struct System
    {
        Set::Field<Set::Scalar> state;
        Set::Field<Set::Scalar> source;
        Set::Field<Set::Scalar> rhs;
        Set::Field<Set::Scalar> mass;
        Set::Field<Set::Scalar> mobility;
        Set::Field<Set::Scalar> tensor_mobility;
        amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
            face_mobility;
        amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
            face_tensor_mobility;
    };

    System& GetSystem(int ncomp);

    int nlevels = 0;
    Set::Scalar tolerance_relative = 1.0e-11;
    Set::Scalar tolerance_absolute = 1.0e-12;
    int verbose = 0;
    int max_order = 2;
    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratio;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> distribution_mapping;
    std::map<int, std::unique_ptr<System>> systems;
};
}

#endif
