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
/// The caller supplies the cell-centered volume source and inverse-density
/// coefficient.  Each newly accumulated cell-centered momentum increment is
/// reconciled by conservative face exchange before the normal face velocity is
/// formed and projected.  This removes the odd-even nullspace of a purely
/// collocated projection while preserving composite linear and angular
/// momentum.  It does not constrain or modify the domain-mean flow.
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
    const amrex::MultiFab& FaceCoefficient(int lev, int direction) const
    {
        return face_coefficient[lev][direction];
    }
    const amrex::MultiFab& FaceVelocity(int lev, int direction) const
    {
        return face_velocity[lev][direction];
    }
    const amrex::MultiFab& Solution(int lev) const { return *solution[lev]; }

    using FaceField =
        amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>;
    using CompositeFaceField = amrex::Vector<FaceField>;

    void PrepareCoefficients(Set::Scalar time);
    void PrepareRHS(int lev, const amrex::MultiFab& velocity, Set::Scalar dt,
                    const FaceField* face_acceleration = nullptr);
    void Solve(Set::Scalar time, const amrex::BCRec& pressure_bc);
    void ApplyCorrection(int lev, amrex::MultiFab& velocity, Set::Scalar dt);
    void ReconcileCellVelocity(Set::Field<Set::Scalar>& velocity,
                               const Set::Field<Set::Scalar>& density,
                               bool projection_increment = false,
                               const CompositeFaceField* face_weight = nullptr);
    void CommitVelocityState(const Set::Field<Set::Scalar>& velocity);

private:
    int nlevels = 0;
    Set::Scalar tolerance_relative = 1.0e-11;
    Set::Scalar tolerance_absolute = 1.0e-12;
    int verbose = 0;
    int max_order = 2;
    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratio;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> distribution_mapping;
    Set::Field<Set::Scalar> solution;
    Set::Field<Set::Scalar> rhs;
    Set::Field<Set::Scalar> coefficient;
    Set::Field<Set::Scalar> divergence;
    Set::Field<Set::Scalar> cell_velocity_predictor;
    Set::Field<Set::Scalar> cell_velocity_reference;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
        face_coefficient;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
        face_velocity;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
        face_velocity_base;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
        face_momentum_exchange;
    bool mixed_velocity_state_initialized = false;
};
}

#endif
