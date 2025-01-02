#ifndef IC_SODSHOCK_H_
#define IC_SODSHOCK_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/ScimitarX_Util.H"


namespace IC {
class SodShock : public IC {
public:
    // Constructor requiring variable indices
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
             const Util::ScimitarX_Util::getVariableIndex& precomputed_indices);

    // Constructor for pressure (no variable indices required)
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name);

    // Add field variables or pressure
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);

private:
    // Initialize the SodShock properties
    void initialize(IO::ParmParse& pp, const std::string& name);

    // MultiFab type name
    std::string mf_name;

    // Variable indices (only required for qvec)
    const Util::ScimitarX_Util::getVariableIndex* variable_indices;
    bool requires_variable_indices;

    // QVec initialization variables
    Set::Scalar rho_left = 1.0;  // Left density
    Set::Scalar rho_right = 0.125; // Right density
    Set::Scalar u_left = 0.0;     // Left velocity (x-direction)
    Set::Scalar u_right = 0.0;    // Right velocity (x-direction)
    Set::Scalar v_left = 0.0;     // Left velocity (y-direction)
    Set::Scalar v_right = 0.0;    // Right velocity (y-direction)
    Set::Scalar w_left = 0.0;     // Left velocity (z-direction)
    Set::Scalar w_right = 0.0;    // Right velocity (z-direction)

    // Pressure initialization variables
    Set::Scalar p_left = 1.0;     // Left pressure
    Set::Scalar p_right = 0.1;    // Right pressure

    // Shock position
    Set::Scalar shock_xpos = 0.5; // Shock position in x-direction
    Set::Scalar shock_ypos = 0.0; // Shock position in y-direction
    Set::Scalar shock_zpos = 0.0; // Shock position in z-direction
};

} // namespace IC

#endif // IC_SODSHOCK_H_

