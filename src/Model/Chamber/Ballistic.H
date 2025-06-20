#ifndef MODEL_CHAMBER_BALLISTIC_H_
#define MODEL_CHAMBER_BALLISTIC_H_

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"
#include <vector>
#include <cmath>

namespace Model
{
namespace Chamber
{
class Ballistic
{
public:
    static constexpr const char* name = "ballistic";
    
    Ballistic() = default;
    Ballistic(IO::ParmParse &pp, std::string name) { pp_queryclass(name,*this); }
    virtual ~Ballistic() = default;


public:
    // Function that computes dpdt using the mass flow rate (mdot)
    [[nodiscard]] Set::Scalar 
    Advance(Set::Scalar dt, Set::Scalar mdot, Set::Scalar volume, Set::Scalar pressure)
    {
        // Calculate mass flow rates
        Set::Scalar mdot_in = mdot;
        
        Util::Message(INFO,"mdot_in",mdot_in);
        Util::Message(INFO,"pressure",pressure);
        Util::Message(INFO,"At",At);
        Util::Message(INFO,"volume",volume);
        Util::Message(INFO,"gamma",gamma);
        Util::Message(INFO,"R",R);
        
        // Set::Scalar c_star = sqrt(gamma * R * T0);
        Set::Scalar c_star = sqrt((gamma / (R * T0))*(pow(2 / (gamma - 1),(gamma + 1.0) / (gamma - 1.0))));
        Set::Scalar mdot_out = (pressure * At) / c_star;
        
        Util::Message(INFO,"c*", c_star);
        Util::Message(INFO,"mdot_out", mdot_out);
        Util::Message(INFO,"mass-flow",mdot_in-mdot_out);
        
        // Compute dp/dt
        Set::Scalar current_dpdt = (mdot_in - mdot_out) / (volume / (R * T0));
        Util::Message(INFO,"dp/dt", current_dpdt);
    
        // Update pressure using Euler integration: p(t+dt) = p(t) + dp/dt * dt
        Set::Scalar new_p = pressure + current_dpdt * dt;
    
        // Store the updated dp/dt and update chamber pressure
        dpdt.push_back(current_dpdt);

        return new_p;
    }

    static void Parse (Ballistic &value, IO::ParmParse &pp)
    { 
        // Chamber pressure [MPa]
        pp.query_default("pressure0",value.pressure0,4.0);
        
        // Throat area [m^2]
        pp.query_default("At",value.At,0.1);
        
        // pp.query_default("vol",value.vol,1.0);
        
        // Specific gas constant for air [J/kg/K]
        pp.query_default("R",value.R, 0.000287);
        
        // Initial temperature [K]
        pp.query_default("T0",value.T0,300.0);
        
        // Specific heat ratio
        pp.query_default("gamma",value.gamma,1.4);
    }


    // Class members for chamber properties
    Set::Scalar pressure0=NAN;  
    Set::Scalar At=NAN;         // Throat area
    Set::Scalar gamma=NAN;      // Specific heat ratio
    Set::Scalar R=NAN;          // Gas constant
    Set::Scalar T0=NAN;         // Initial temperature
    Set::Scalar c_star=NAN;     // Characteristic Velocity
    // Set::Scalar vol=NAN;     // Chamber volume

    std::vector<Set::Scalar> dpdt;  // dp/dt vector to store rate of change of pressure over time
};

}
}

#endif
