#ifndef INTEGRATOR_NARROWBANDLEVELSET_H_
#define INTEGRATOR_NARROWBANDLEVELSET_H_

// Standard includes
#include <memory>
#include <vector>

// Alamo includes
#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"
#include "IC/IC.H"

// Fowarard declarations
struct SimulationDomain;
struct LevelSetField;


namespace Integrator
{

class NarrowBandLevelset : virtual public Integrator
{
public:
    // --------------------------
    // Integrator identification
    // --------------------------
    static constexpr const char* name = "narrowbandlevelset";

    // Constructors and destructor
    NarrowBandLevelset();
    NarrowBandLevelset(IO::ParmParse& pp);

    static void Parse(NarrowBandLevelset& value, IO::ParmParse& pp);

    virtual ~NarrowBandLevelset();

protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;

    // --------------------------
    // Multiphysics coupling API
    // --------------------------

    // // Level set value φ
    // const Set::Field<Set::Scalar>& phi(int field = 0) const;

    // // Extension / object velocity component
    // const Set::Field<Set::Scalar>& velocity(int dir, int field = 0) const;

    // // Object ID per cell (-1 = no object)
    // const amrex::iMultiFab& object_id(int field = 0) const;

private:
    // ============================================================
    // Internal data structures (implementation details)
    // ============================================================

    // Helper functions
    void RegisterOneField(
        LevelSetField& ls,
        const std::string& prefix
    );

    // ============================================================
    // Ownership
    // ============================================================

    // Grid, geometry, ghosting, field ownership
    std::unique_ptr<SimulationDomain> domain_;

    // CFL for time computation
    Set::Scalar cfl = NAN;
};

} // namespace Integrator

#endif // INTEGRATOR_NARROWBANDLEVELSET