#ifndef MODEL_CHAMBER_BALLISTIC_H_
#define MODEL_CHAMBER_BALLISTIC_H_


// TODO: Imports
#include "Util/Util.H"
#include "IO/ParmParse.H"
#include <vector>
#include <cmath>
#include <utility>


namespace Model
{
namespace Chamber
{
class Ballistic
{
public:
    static constexpr const char* name = "ballistic";

    Ballistic() = default;
    Ballistic(IO::ParmParse &pp, std::string name) { pp_queryclass(name, *this); }
    virtual ~Ballistic() = default;

public:
    // Advances chamber pressure by one timestep.
    // Implements Eq. 13 from Meier et al. (2024):
    //
    //   dp0/dt = (R*T0/V) * [mdot_in - p0*At*sqrt(gamma/(R*T0) * (2/(gamma+1))^((gamma+1)/(gamma-1)))]
    //
    // Caller must compute: mdot_in = r * Ab * (rho_s - rho_0)
    // where rho_0 = p0/(R*T0)
    [[nodiscard]] std::pair<Set::Scalar, Set::Scalar>
    Advance(Set::Scalar dt, Set::Scalar mdot_in, Set::Scalar volume, Set::Scalar pressure)
    {
        // inverse C star: sqrt(gamma/(R*T0) * (2/(gamma+1))^((gamma+1)/(gamma-1)))
        Set::Scalar exponent = (gamma + 1.0) / (gamma - 1.0);
        Set::Scalar inv_c_star = std::sqrt((gamma / (R * T0)) * std::pow(2.0 / (gamma + 1.0), exponent));

        // Mass flow out through choked nozzle
        Set::Scalar mdot_out = pressure * At * inv_c_star;

        // Pressure rate of change
        Set::Scalar current_dpdt = (R * T0 / volume) * (mdot_in - mdot_out);

        // Forward Euler integration
        Set::Scalar new_pressure = pressure + current_dpdt * dt;

        dpdt.push_back(current_dpdt);
        return {new_pressure, current_dpdt};
    }



    static void Parse(Ballistic &value, IO::ParmParse &pp)
    {
        // Initial chamber pressure [Pa]
        pp.query_default("pressure", value.pressure, "4.0e6_Pa", Unit::Pressure());

        // Nozzle throat area [m^2]
        pp.query_default("At", value.At, "0.01_m^2", Unit::Area());

        // Specific gas constant [J/(kg*K)]
        pp.query_default("R", value.R, "287.0_J/kg/K",
                         Unit::Energy() / Unit::Mass() / Unit::Temperature());

        // Chamber stagnation temperature [K]
        pp.query_default("T0", value.T0, "300.0_K", Unit::Temperature());

        // Ratio of specific heats [-]
        pp.query_default("gamma", value.gamma, "1.4", Unit::Less());
    }

    Set::Scalar pressure = NAN;
    Set::Scalar At = NAN;
    Set::Scalar gamma = NAN;
    Set::Scalar R = NAN;
    Set::Scalar T0 = NAN;

    std::vector<Set::Scalar> dpdt;
};

} // namespace Chamber
} // namespace Model

#endif