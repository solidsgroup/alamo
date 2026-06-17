#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#include "Util/Util.H"
#include "Unit.H"
#include <stdexcept>

namespace UnitTest
{
int Equivalence(int verbose = 0)
{

    std::vector<std::pair<std::string,std::string>> unit_tests = {

        // Length
        {"1_m", "1000_mm"},
        {"1_m", "1e9_nm"},
        {"1_cm", "10_mm"},
        {"1_km", "1000_m"},
        {"1_in", "0.0254_m"},
        {"1_ft", "0.3048_m"},

        // Area
        {"1_m^2", "1e6_mm^2"},
        {"1_cm^2", "100_mm^2"},

        // Volume
        {"1_m^3", "1e9_mm^3"},
        {"1_cm^3", "1000_mm^3"},
        {"1_L", "1e-3_m^3"},
        {"1_mL", "1e-6_m^3"},
        {"1_cc", "1e-6_m^3"},

        // Time
        {"1_s", "1000_ms"},
        {"1_s", "1e9_ns"},
        {"1_min", "60_s"},
        {"1_hr", "3600_s"},

        // Mass
        {"1_kg", "1000_g"},
        {"1_g", "1000_mg"},
        {"1_tonne", "1000_kg"},
        {"1_lb", "0.45359237_kg"},
        {"1_oz", "0.0283495_kg"},

        // Force
        {"1_N", "1_kg*m/s^2"},
        {"1_kN", "1000_N"},
        {"1_lbf", "4.44822162_N"},

        // Pressure
        {"1_Pa", "1_N/m^2"},
        {"1_kPa", "1000_Pa"},
        {"1_MPa", "1e6_Pa"},
        {"1_bar", "1e5_Pa"},
        {"1_psi", "6894.76_Pa"},

        // Energy
        {"1_J", "1_N*m"},
        {"1_mJ", "1e-3_J"},
        {"1_kJ", "1e3_J"},
        {"1_eV", "1.602176634e-19_J"},

        // Power
        {"1_W", "1_J/s"},
        {"1_kW", "1e3_W"},
        {"1_MW", "1e6_W"},

        // Mobility / Phase-field (example)
        {"1_m/Pa/s", "1e+3_mm/kPa/ms"},   
        {"1_m^3/J/s", "1e+9_mm^3/kJ/ms"}, 
        {"1_m^3/Pa/s", "1e+9_mm^3/kPa/ms"}, 

        // Thermal
        {"1_W/m/K", "1e3_mW/m/K"}, // thermal conductivity
        {"1_m^2/s", "1e6_mm^2/s"}, // diffusivity

        // Velocity / Kinematics
        {"1_m/s", "1e3_mm/s"},
        {"1_m/s^2", "1e3_mm/s^2"},
        {"1_km/hr", "0.27777778_m/s"},

        // Frequency
        {"1_Hz", "1/s"},
        {"1_kHz", "1e3_Hz"},
        {"1_MHz", "1e6_Hz"},

        // Charge / Current
        {"1_C", "1_A*s"},
        {"1_mAh", "3.6_C"},
        {"1_A", "1_C/s"},

        // Pressure/Stress alternative forms
        {"1_Pa", "1_N/m^2"},
        {"1_N/m^2", "1_kg/m/s^2"},

        // Force / Energy conversions
        {"1_J", "1_kg*m^2/s^2"},
        {"1_N", "1_kg*m/s^2"},

        // Phase-field driving forces
        {"1_J/m^3", "1_Pa"}, // energy density = stress
        {"1_N/m^2", "1_Pa"},

        // Example: mobility, dimensionally correct
        {"1E-12_m/Pa/s", "1E-9_mm/kPa/ms"},
        {"1E-15_m^3/J/s", "1E-15_mm^3/J/ms"},

    };


    int fails = 0;
    for (auto unit_test : unit_tests)
    {
        Unit a=Unit::Less(), b = Unit::Less();
        try
        {
            a = Unit::Parse(unit_test.first);
        }
        catch(std::runtime_error &e)
        {
            if (verbose)
            {
                Util::Message(INFO,"Error parsing ", unit_test.first);
                Util::Message(INFO,e.what());
            }
            fails++;
            continue;
        }
        try
        {
            b = Unit::Parse(unit_test.second);
        }
        catch(std::runtime_error &e)
        {
            if (verbose)
            {
                Util::Message(INFO,"Error parsing ", unit_test.second);
                Util::Message(INFO,e.what());
            }
            fails++;
            continue;
        }

        if (! a.isType(b))
        {
            if (verbose)
            {
                Unit::Parse(unit_test.first,true);
                Unit::Parse(unit_test.second,true);
                Util::Message(INFO,unit_test.first," different type than ",unit_test.second);
                Util::Message(INFO,"   ",unit_test.first," ==> ",a);
                Util::Message(INFO,"   ",unit_test.second," ==> ",b);
            }
            fails++;
            continue;
        }
        if (fabs(a.normalized_value()-b.normalized_value()) > 1E-8)
        {
            if (verbose)
            {
                Unit::Parse(unit_test.first,true);
                Unit::Parse(unit_test.second,true);
                Util::Message(INFO,unit_test.first," resolved value does not match ",unit_test.second);
                Util::Message(INFO,"   ",unit_test.first," ==> ",a);
                Util::Message(INFO,"   ",unit_test.second," ==> ",b);
            }
            fails++;
            continue;
        }
    }

    return fails;
}
}

#endif
