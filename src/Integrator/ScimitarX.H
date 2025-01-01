#ifndef INTEGRATOR_SCIMITARX_H // Include guards
#define INTEGRATOR_SCIMITARX_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/SodShock.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"
#include "Util/Util.H"
#include "Util/ScimitarX_SetIndex.H"


namespace Integrator
{
class ScimitarX : virtual public Integrator
{
public:


    ScimitarX() {} ; 
    ScimitarX(IO::ParmParse& pp);   
    static void Parse(ScimitarX& value, IO::ParmParse& pp);
    virtual ~ScimitarX()
    {
        delete bc_PVec;
        delete bc_Pressure;

        delete ic_PVec;
        delete ic_Pressure;

    }

    std::string constructKey(const std::string& prefix, const std::string& component_name, const std::string& property, const std::string& side) {
    return prefix + component_name + "." + property + "." + side;
    }

    
   static IO::ParmParse setupPVecBoundaryConditions(IO::ParmParse& pp, Util::SolverType solverType, Util::ScimitarX_SetIndex::VariableIndicesResult variableIndicesResult);
protected:

    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;

protected:
    Set::Field<Set::Scalar> QVec_mf;          // Q Vector holding all the field variables in Conservative Form (current timestep)
    Set::Field<Set::Scalar> QVec_old_mf;      // Q Vector holding all the  field variable in Conservative Form (previous timestep)
    Set::Field<Set::Scalar> XFlux_mf;         // X-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> YFlux_mf;         // Y-dir Flux Vector holding all the  field variable (current timestep)
    Set::Field<Set::Scalar> ZFlux_mf;         // Z-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> Source_mf;         // Source term Vector for all the field variables (current timestep)
    Set::Field<Set::Scalar> PVec_mf;      // Primitive Vector holding all the  field variables  in primitive form (current timestep)
    Set::Field<Set::Scalar> Pressure_mf;  // Multifab for pressure variable   

private:
    int number_of_components;            // Number of components
    int number_of_ghost_cells = 4;           // Number of ghost cells

    Set::Scalar refinement_threshold;

    IC::IC* ic_PVec = nullptr;            // Object used to initialize Q Vector
    IC::IC* ic_Pressure = nullptr;       //  Object used to initialize Pressure
                                                      //  
    BC::BC<Set::Scalar>* bc_PVec = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::BC<Set::Scalar>* bc_Pressure = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::Nothing bc_nothing;                                                       //
                                                            

  }; //ScimitarX Class
} // namespace Integrator
#endif

