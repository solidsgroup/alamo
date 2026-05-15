#ifndef INTEGRATOR_DOUBLYDEGENERATECAHNHILLIARD_H
#define INTEGRATOR_DOUBLYDEGENERATECAHNHILLIARD_H

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
    static constexpr const char *name = "doublydegeneratecahnhilliard";

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
    /// Integrate initial thermo output
    void TimeStepBegin(Set::Scalar time, int step) override;
    /// Integrate thermo output after each completed timestep
    void TimeStepComplete(Set::Scalar time, int step) override;
    /// Integrate eta over one timestep on lev
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    /// Integrate eta volume for thermo output
    void Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &box) override;
    /// Mark any cells that need to be refined
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow) override;

private:
    Set::Field<Set::Scalar> etanew_mf;                  /// The new value for eta this timestep
    Set::Field<Set::Scalar> etaold_mf;                  /// Last timestep's value for eta
    Set::Field<Set::Scalar> deta_mf;                    /// Change in eta from the previous timestep
    Set::Field<Set::Scalar> driving_force_mf;           /// The variational derivative of the free energy functional
    Set::Field<Set::Vector> evolution_driving_force_mf; /// The product of the degenerate mobility with the gradient of the driving force

    BC::BC<Set::Scalar> *bc; /// eta's bc object
    IC::IC<Set::Scalar> *ic; /// eta's ic object

    Set::Scalar epsilon = NAN;
    Set::Scalar omega = NAN;
    Set::Scalar gamma = NAN;
    // Set::Scalar p = NAN;
    Set::Scalar alpha = NAN;
    Set::Scalar mu = NAN;
    Set::Scalar volume = 0.0;
};
}
#endif
