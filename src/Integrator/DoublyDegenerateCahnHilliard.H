#ifndef INTEGRATOR_CAHNHILLIARD_H
#define INTEGRATOR_CAHNHILLIARD_H

#include <AMReX.H>
#include <AMReX_MLMG.H>

#include "BC/BC.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

namespace Integrator
{

class DoublyDegenerateCahnHilliard : public Integrator
{
public:
    static constexpr const char *name = "cahnhilliard";

    /// Basic constructor (don't use)
    DoublyDegenerateCahnHilliard();

    /// Destroy pointers defined in Parse
    ~DoublyDegenerateCahnHilliard();

    /// Use this constructor
    DoublyDegenerateCahnHilliard(IO::ParmParse &pp) : DoublyDegenerateCahnHilliard()
    {
        Parse(*this, pp);
    }

    /// Scan input values and initialize fields
    static void Parse(DoublyDegenerateCahnHilliard &value, IO::ParmParse &pp);

protected:
    /// Set values in fields
    void Initialize(int lev) override;
    /// Integrate eta over one timestep on lev
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    /// Mark any cells that need to be refined
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow) override;

private:
    Set::Field<Set::Scalar> etanew_mf;    /// The new value for eta this timestep
    Set::Field<Set::Scalar> etaold_mf;    /// Last timestep's value for eta
    Set::Field<Set::Scalar> grad_etaold_mf;
    Set::Field<Set::Scalar> lap_etaold_mf;

    BC::BC<Set::Scalar> *bc; /// eta's bc object
    IC::IC<Set::Scalar> *ic; /// eta's ic object

    Set::Scalar epsilon = NAN;
    Set::Scalar omega = NAN;
    Set::Scalar gamma = NAN;
    // Set::Scalar p = NAN;
    Set::Scalar alpha = NAN;
    Set::Scalar mu = NAN;
    Set::Scalar refinement_threshold = NAN;
};
}
#endif
