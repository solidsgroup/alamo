#ifndef INTEGRATOR_NARROWBANDLEVELSET_H // Include guards
#define INTEGRATOR_NARROWBANDLEVELSET_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <memory>

// Alamo includes
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"

// Numeric includes
//#include "Numeric/FluxHandler.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
class NarrowBandLevelset : virtual public Integrator
{
public:

    // Constructors
    NarrowBandLevelset(int a_nghost = 2); // Default constructor
    NarrowBandLevelset(IO::ParmParse& pp); // Constructor that triggers parse
    static void Parse(NarrowBandLevelset& value, IO::ParmParse& pp);

    // Destructor
    virtual ~NarrowBandLevelset()
    {
        delete ls_ic;
        delete velocity_ic;
        delete ls_bc;
    }
    
    // Public getters to access proteted/private members in Numeric templates
    int GetNumGhostCells() const { return number_of_ghost_cells; }
    const Set::Field<Set::Scalar>& GetLevelSetField() const { return ls_mf; }
    const Set::Field<Set::Scalar>& GetLevelSetFieldOld() const { return ls_old_mf; }
    const Set::Field<Set::Scalar>& GetVelocityField() const { return velocity_mf; }
    
protected:
    // Variables (multi-fabs)
    Set::Field<Set::Scalar>    ls_mf;                 // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    ls_old_mf;             // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar>    velocity_mf;           // Velocity field variable (Assumed constant for now)

    // Solver setup and control functions (override)
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;
    
    // Additional functions
    void Advect(int lev, Set::Scalar dt);
    void UpdateInterfaceVelocity(int lev);
    void Reinitialize(int lev);

// Declare private functions such as Parse 
private:
    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    int current_timestep            = 0;
    const int number_of_components  = 1;     // Number of components (1-levelset value)
    const int number_of_ghost_cells = 2;     // Number of ghost cells

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC* ls_ic = nullptr;                     // Object used to initialize levelset field
    IC::IC* velocity_ic = nullptr;               // Object used to initialize velocity field
    BC::BC<Set::Scalar>* ls_bc = nullptr;        // Object used to update temp levelset field boundary ghost cells
    BC::BC<Set::Scalar>* velocity_bc = nullptr;  // Object used to update temp velocity field boundary ghost cells

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif
