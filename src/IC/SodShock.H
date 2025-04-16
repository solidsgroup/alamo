#ifndef IC_SODSHOCK_H_
#define IC_SODSHOCK_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/ScimitarX_Util.H"
#include <vector> // Include the standard vector header

namespace IC {

class SodShock : public IC<Set::Scalar> {
public:
    // Constructor requiring variable indices
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
             const Util::ScimitarX_Util::getVariableIndex& precomputed_indices);

    // Constructor for pressure (no variable indices required)
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name);

    // Add field variables or pressure
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);

protected:
    // Initialize the SodShock properties
    void initialize(IO::ParmParse& pp, const std::string& name);

    // MultiFab type name
    std::string mf_name;

    // Variable indices (only required for qvec)
    const Util::ScimitarX_Util::getVariableIndex* variable_indices;
    bool requires_variable_indices;

    // Number of shocks (zones)
    int number_of_shocks;   

    // Shock positions in x-direction
    std::vector<Set::Scalar> shock_positions;
    Set::Scalar shock_xpos = 0.5; // Shock position in x-direction
    Set::Scalar shock_ypos = 0.0; // Shock position in y-direction
    Set::Scalar shock_zpos = 0.0; // Shock position in z-direction

private:
    // QVec initialization variables
    std::vector<Set::Scalar> density_zone;
    std::vector<Set::Scalar> uvel_zone;
    std::vector<Set::Scalar> vvel_zone;
    std::vector<Set::Scalar> wvel_zone;

    // Pressure initialization variables
    std::vector<Set::Scalar> pressure_zone;
};

} // namespace IC

#endif // IC_SODSHOCK_H_
