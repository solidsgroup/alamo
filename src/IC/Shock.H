#ifndef IC_SHOCK_H_
#define IC_SHOCK_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/ScimitarX_Util.H"
#include <vector> 
#include "Numeric/IntegratorVariableAccessLayer.H"


namespace IC {
class Shock : public IC<Set::Scalar> {
public:
    // Constructor requiring variable indices
    Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
            const Numeric::GenericVariableAccessor::VariableIndices& precomputed_indices);

    // Constructor for pressure (no variable indices required)
    Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name);

    // Add field variables or pressure
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);

protected:

    // Initialize the Shock properties
    void initialize(IO::ParmParse& pp, const std::string& name);

    // MultiFab type name
    std::string mf_name;

    // Variable indices (only required for qvec)
    const Numeric::GenericVariableAccessor::VariableIndices* variable_indices;
    bool requires_variable_indices;

    // Number of shocks (zones)
    int number_of_shocks;   

    // Shock positions in x-direction
    std::vector<Set::Scalar> shock_positions;

    // Shock directions (xdir, ydir and zdir)
    std::string shock_direction;

    int direction_index;  // Direction indices (0=x, 1=y, 2=z)
    std::vector<Set::Scalar> zone_lower_bounds;  // Pre-computed lower bounds for each zone
    std::vector<Set::Scalar> zone_upper_bounds;  // Pre-computed upper bounds for each zone


private:
    // pvec initialization variables
    std::vector<Set::Scalar> density_zone;
    std::vector<Set::Scalar> uvel_zone;
    std::vector<Set::Scalar> vvel_zone;
    std::vector<Set::Scalar> wvel_zone;
    // pressure initialization variables
    std::vector<Set::Scalar> pressure_zone;
};

} // namespace IC

#endif // IC_SHOCK_H_
