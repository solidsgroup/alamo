#ifndef MODEL_BALLISTICS_BALLISTIC_NEW_H_
#define MODEL_BALLISTICS_BALLISTIC_NEW_H_

#include "IO/ParmParse.H"
//#include "Integrator/Flame.cpp"
#include <vector>      // Include vector for dynamic array
#include <cmath>       // Include cmath for sqrt

namespace Model
{
namespace Ballistics
{
class Ballistic_new 
{
public:
    Ballistic_new()
    {
        // Initialize constants and initial values
        DefineConstants();
        pressure0=4.0;
        dpdt.push_back(0.0);  // Initialize dp/dt as zero
    }

    virtual ~Ballistic_new() {}

    // Main function to execute the calculation at each time step
    void Execute(Set::Scalar mdot, Set::Scalar dt)
    {
        Compute(mdot, dt);  // Perform calculations using mdot and time step dt
    }

public:
    // Function to define constants like gamma, At, R, T0, etc.
    void DefineConstants()
    {
        // Define constants here
        gamma = 1.4;    // Example value for the specific heat ratio
        At = 1.0;       // Example value for throat area (m^2)
        R = 287.0;        // Specific gas constant for air in J/(kg·K)
        T0 = 300.0;       // Initial temperature in K
        volume = 1.0;   // Volume in cubic meters
    }

    // Function that computes dpdt using the mass flow rate (mdot)
    void Compute(Set::Scalar mdot, Set::Scalar dt)
{
    // Get the current pressure                
    Set::Scalar current_p = pressure0;

    // Calculate mass flow rates
    Set::Scalar mdot_in = mdot;
    
    Set::Scalar mdot_out = (current_p * At) / sqrt(gamma * R * T0);

    // Compute dp/dt
    Set::Scalar current_dpdt = (mdot_in - mdot_out) / (volume / (R * T0));

    // Update pressure using Euler integration: p(t+dt) = p(t) + dp/dt * dt
    Set::Scalar new_p = current_p + current_dpdt * dt;

    // Store the updated dp/dt and update chamber pressure
    dpdt.push_back(current_dpdt);
    pressure0 = new_p;
    
}


    // Class members for chamber properties
    Set::Scalar pressure0=NAN;  // Chamber pressure
    Set::Scalar At=NAN;         // Throat area
    Set::Scalar gamma=NAN;      // Specific heat ratio
    Set::Scalar R=NAN;          // Gas constant
    Set::Scalar T0=NAN;         // Initial temperature
    Set::Scalar volume=NAN;     // Chamber volume

    std::vector<Set::Scalar> dpdt;  // dp/dt vector to store rate of change of pressure over time
};

} // namespace Ballistics
} // namespace Model

#endif // MODEL_BALLISTICS_BALLISTIC_NEW_H_
