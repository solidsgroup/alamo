#ifndef UNIT_UNIT_H
#define UNIT_UNIT_H


#include <array>
#include <charconv>
#include <cstdlib>
#include <iostream>
#include <system_error>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <cmath>
#include <regex>

//#include "Set/Base.H"
//#include "String.H"

struct Unit : std::pair<double,std::array<int,7>>
{

    //
    // STATIC CONST DEFINITIONS, TRUE EVERYWHERE
    // 
    enum class Type {
        Length=0, Time=1, Mass=2, Temperature=3, Current=4, Amount=5, LuminousIntensity=6
    };    
    // Fundamental units
    inline static const std::map<std::string, std::pair<double,std::array<int,7>>> base_units =
    {
        {"m",  {1.0, {1,0,0,0,0,0,0}}},
        {"s",  {1.0, {0,1,0,0,0,0,0}}},
        {"kg", {1.0, {0,0,1,0,0,0,0}}},
        {"K",  {1.0, {0,0,0,1,0,0,0}}},
        {"A",  {1.0, {0,0,0,0,1,0,0}}},
        {"mol",{1.0, {0,0,0,0,0,1,0}}},
        {"cd", {1.0, {0,0,0,0,0,0,1}}}
    };
    //
    // Compound unit relationships and nomenclature
    // Any desired unit relationship can be added to this list.
    // The unit conversion can be specified in any desired format, as long as
    // all of the units on the right hand side are resolvable.
    //
    inline static const std::map<std::string, std::pair<double,std::string>> compound = {
        // Lengths
        {"km",        {1e+3,    "m"        }},
        {"hm",        {1e+2,    "m"        }},
        {"dam",       {1e+1,    "m"        }},
        {"dm",        {1e-1,    "m"        }},
        {"cm",        {1e-2,    "m"        }},
        {"mm",        {1e-3,    "m"        }},
        {"um",        {1e-6,    "m"        }},
        {"nm",        {1e-9,    "m"        }},
        {"pm",        {1e-12,   "m"        }},
        {"in",        {0.0254,  "m"        }},
        {"ft",        {0.3048,  "m"        }},
        {"yd",        {0.9144,  "m"        }},
        {"mi",        {5280.0,  "ft"       }},
        {"furlong",   {660.0,   "ft"       }},
        {"au",        {1.495978707e11, "m" }}, // astronomical unit
        {"ly",        {9.4607304725808e15, "m"}}, // lightyear
        {"pc",        {3.08567758149137e16, "m"}}, // parsec

        // Times
        {"ms",        {1e-3,    "s"        }},
        {"us",        {1e-6,    "s"        }},
        {"ns",        {1e-9,    "s"        }},
        {"ps",        {1e-12,   "s"        }},
        {"min",       {60.0,    "s"        }},
        {"hr",        {60.0,    "min"      }},
        {"day",       {24.0,    "hr"       }},
        {"week",      {7.0,     "day"      }},
        {"fortnight", {14.0,    "day"      }},
        {"year",      {365.25,  "day"      }},
        {"decade",    {10.0,    "year"     }},

        // Mass
        {"g",         {1e-3,    "kg"       }},
        {"mg",        {1e-6,    "kg"       }},
        {"ug",        {1e-9,    "kg"       }},
        {"tonne",     {1000.0,  "kg"       }},
        {"lb",        {0.45359237, "kg"    }},
        {"oz",        {0.0283495, "kg"     }},
        {"slug",      {14.5939, "kg"       }},

        // Speed
        {"mph",       {1.0,     "mi/hr"    }},
        {"kph",       {1.0,     "km/hr"    }},
        {"mps",       {1.0,     "m/s"      }},
        {"fps",       {1.0,     "ft/s"     }},
        {"knot",      {1852.0 / 3600.0, "m/s"}},

        // Area
        {"ha",        {1e4,     "m^2"      }},
        {"acre",      {4046.85642, "m^2"   }},
        {"sqft",      {1.0,     "ft^2"     }},
        {"sqin",      {1.0,     "in^2"     }},
        {"sqkm",      {1.0,     "km^2"     }},
        {"sqmi",      {1.0,     "mi^2"     }},

        // Volume
        {"L",         {1e-3,    "m^3"      }},
        {"mL",        {1e-6,    "m^3"      }},
        {"cc",        {1e-6,    "m^3"      }},
        {"gal",       {3.78541e-3, "m^3"   }},
        {"qt",        {0.946353e-3, "m^3"  }},
        {"pt",        {0.473176e-3, "m^3"  }},
        {"cup",       {0.24e-3, "m^3"      }},
        {"floz",      {29.5735e-6, "m^3"   }},
        {"tbsp",      {14.7868e-6, "m^3"   }},
        {"tsp",       {4.92892e-6, "m^3"   }},

        // Force
        {"N",         {1.0,     "kg*m/s^2" }},
        {"kN",        {1e3,     "N"        }},
        {"lbf",       {4.44822162, "N"     }},
        {"dyn",       {1e-5,    "N"        }},

        // Energy
        {"J",         {1.0,     "N*m"      }},
        {"kJ",        {1e3,     "J"        }},
        {"MJ",        {1e6,     "J"        }},
        {"eV",        {1.602176634e-19, "J"}},
        {"cal",       {4.184,   "J"        }},
        {"kcal",      {4184.0,  "J"        }},
        {"BTU",       {1055.06, "J"        }},
        {"ft*lbf",    {1.35582, "J"        }},

        // Power
        {"W",         {1.0,     "J/s"      }},
        {"kW",        {1e3,     "W"        }},
        {"MW",        {1e6,     "W"        }},
        {"hp",        {745.7,   "W"        }},

        // Pressure
        {"Pa",        {1.0,     "N/m^2"    }},
        {"kPa",       {1e3,     "Pa"       }},
        {"MPa",       {1e6,     "Pa"       }},
        {"GPa",       {1e9,     "Pa"       }},
        {"bar",       {1e5,     "Pa"       }},
        {"atm",       {101325.0,"Pa"       }},
        {"psi",       {6894.76, "Pa"       }},
        {"mmHg",      {133.322, "Pa"       }},
        {"torr",      {133.322, "Pa"       }},

        // Temperature differences (not absolute temperatures)
        {"K",         {1.0,     "K"        }},
        {"degR",      {5.0/9.0, "K"        }},

        // Charge
        {"C",         {1.0,     "A*s"      }},
        {"mAh",       {3.6,     "C"        }},

        // Voltage
        {"V",         {1.0,     "J/C"      }},

        // Capacitance
        {"F",         {1.0,     "C/V"      }},

        // Resistance
        {"ohm",       {1.0,     "V/A"      }},
        {"kOhm",      {1e3,     "ohm"      }},
        {"MOhm",      {1e6,     "ohm"      }},

        // Conductance
        {"S",         {1.0,     "1/ohm"    }},

        // Magnetic
        {"T",         {1.0,     "kg/(A*s^2)"}},
        {"G",         {1e-4,    "T"        }}, // gauss

        // Frequency
        {"Hz",        {1.0,     "1/s"      }},
        {"kHz",       {1e3,     "Hz"       }},
        {"MHz",       {1e6,     "Hz"       }},
        {"GHz",       {1e9,     "Hz"       }},
    };


    //
    // These static factory functions return unit types.
    // They are used to check unit compliance, mainly at the parsing stage.
    // 
    // Primary
    static Unit Less()              { return Unit{1.0, {0,0,0,0,0,0,0}}; } // 1
    static Unit Length()            { return Unit{1.0, {1,0,0,0,0,0,0}}; } // m
    static Unit Time()              { return Unit{1.0, {0,1,0,0,0,0,0}}; } // s
    static Unit Mass()              { return Unit{1.0, {0,0,1,0,0,0,0}}; } // kg
    static Unit Temperature()       { return Unit{1.0, {0,0,0,1,0,0,0}}; } // K
    static Unit Current()           { return Unit{1.0, {0,0,0,0,1,0,0}}; } // A
    static Unit Amount()            { return Unit{1.0, {0,0,0,0,0,1,0}}; } // mol
    static Unit LuminousIntensity() { return Unit{1.0, {0,0,0,0,0,0,1}}; } // cd

    // Geometry and kinematics
    static Unit Area() { return Length()*Length(); }
    static Unit Volume() { return Length()*Length()*Length(); }
    static Unit Velocity() {return Length() / Time();}
    static Unit Acceleration(){ return Velocity() / Time(); }
    // Mechanics
    static Unit Force()       { return Mass() * Acceleration(); } // kg·m/s² = N
    static Unit Momentum()    { return Mass() * Velocity(); }
    static Unit Impulse()     { return Force() * Time(); }
    static Unit Pressure()    { return Force() / Area(); }         // N/m² = Pa
    static Unit Energy()      { return Force() * Length(); }       // N·m = J
    static Unit Power()       { return Energy() / Time(); }        // J/s = W
    static Unit Density()     { return Mass() / Volume(); }
    static Unit SpecificWeight() { return Force() / Volume(); }
    static Unit Work()        { return Energy(); }
    // Thermodynamics
    static Unit SpecificHeatCapacity() {return Energy() / Mass() / Temperature() ;}
    static Unit ThermalConductivity() {return Power()/Length()/Temperature();}
    static Unit ThermalDiffusivity() {return ThermalConductivity() / Density() / SpecificHeatCapacity();}
    static Unit HeatFlux()             { return Power() / Area(); }
    static Unit HeatTransferCoefficient() { return HeatFlux() / Temperature(); }
    // Electromagnetism

    static Unit Charge()      { return Current() * Time(); }                    // A·s = C
    static Unit Voltage()     { return Power() / Current(); }                   // W/A = V
    static Unit Resistance()  { return Voltage() / Current(); }                 // V/A = Ω
    static Unit Capacitance() { return Charge() / Voltage(); }                  // C/V = F
    static Unit Conductance() { return Less() / Resistance(); }                 // S = 1/Ω
    static Unit MagneticFlux(){ return Voltage() * Time(); }                    // V·s = Wb
    static Unit MagneticField(){ return MagneticFlux() / Area(); }              // Wb/m² = T
    static Unit Inductance()  { return MagneticFlux() / Current(); }            // Wb/A = H
    // Chemistry
    static Unit MolarMass()   { return Mass() / Amount(); }                     // kg/mol
    static Unit Concentration(){ return Amount() / Volume(); }                  // mol/m³
    static Unit MolarEnergy() { return Energy() / Amount(); }                   // J/mol
    // Optics and radiation
    static Unit LuminousFlux()      { return LuminousIntensity(); }             // for scalar usage (lm)
    static Unit Illuminance()       { return LuminousFlux() / Area(); }         // lux = lm/m²


    //
    // This is where SYSTEM NORMALIZATION is stored.
    // Should usually be set once (at the beginning of the simulation) and then
    // left constant throughout the simulation.
    // All normalizations are in terms of standard SI units: m,s,kg,K,A,mol,cd.
    //
    inline static std::array<std::pair<std::string,double>, 7> normalization = {
        std::pair{"m",  1.0},
        std::pair{"s",  1.0},
        std::pair{"kg", 1.0},
        std::pair{"K",  1.0},
        std::pair{"A",  1.0},
        std::pair{"mol",1.0},
        std::pair{"cd", 1.0}
    };

    using std::pair<double, std::array<int,7>>::pair;
    Unit(const std::pair<double, std::array<int, 7>>& p) : std::pair<double, std::array<int, 7>>(p) {}
    static std::map<std::string, int> ParseUnitString(const std::string& input);

    
    static Unit Parse(double val, std::string unitstring)
    {
        Unit ret = StringParse(unitstring);
        ret.first *= val;
        return ret;
    }
    static Unit Parse(std::string unit)
    {
        std::vector<std::string> split;
        std::stringstream ss(unit);
        std::string item;
        while (std::getline(ss, item, '_')) {
            split.push_back(item);
        }
        if (split.size() == 2)
        {        
            return Parse(std::stod(split[0]),split[1]);
        }
        else if (split.size() == 1)
        {
            // Attempt to parse the string as if it's a unitless number
            double value = NAN;
            auto first = unit.data();
            auto last = unit.data() + unit.size();
            auto result = std::from_chars(first, last, value);
            
            // If it is a untless number, then return a unitless result
            if (result.ec == std::errc() && result.ptr == last)
            {
                Unit ret = Unit::Less();
                ret.first = value;
                return ret;
            }
            // Otherwise, it must be a pure unit, so return that.
            else
            {
                return StringParse(unit);
            }
        }
        else
        {
            throw std::runtime_error("Invalid unit string "+unit);
        }
    }


    static Unit StringParse(std::string unitstring)
    {
        Unit type = Unit::Less();
        std::map<std::string, int> map = ParseUnitString(unitstring);
        for (auto it = map.begin(); it != map.end(); it++)
        {
            auto [str, exp] = *it;

            if (Unit::base_units.find(str) != base_units.end())
            {
                Unit t = Unit::base_units.at(str);
                type = type * (t^exp);
            }
            else if (compound.find(it->first) != compound.end())
            {
                auto [subfac, substr] = compound.at(str);
                Unit t = StringParse(substr);
                t.first *= subfac;
                type = type * t^exp;
            }
            else
            {
                throw std::runtime_error("Unknown unit: " + unitstring);
            }
        }
        return type;
    }


    static void setLengthUnit(std::string unit)
    {
        Unit parsed = StringParse(unit);
        if (parsed.second != Length().second) throw std::runtime_error("Not a unit of length: " + unit);
        normalization[int(Type::Length)].first = unit;
        normalization[int(Type::Length)].second = 1.0 / parsed.first;
    }
    static void setTimeUnit(std::string unit)
    {
        Unit parsed = StringParse(unit);
        if (parsed.second != Time().second) throw std::runtime_error("Not a unit of time: " + unit);
        normalization[int(Type::Time)].first = unit;
        normalization[int(Type::Time)].second = 1.0 / parsed.first;
    }

    bool isType(const Unit &test) const
    {
        return this->second == test.second;
    }

    Unit operator * (const Unit &rhs)
    {
        Unit ret;
        ret.first = this->first * rhs.first;
        for (unsigned int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] + rhs.second[i];
        return ret;
    }
    Unit operator / (const Unit &rhs)
    {
        Unit ret;
        ret.first = this->first / rhs.first;
        for (unsigned int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] - rhs.second[i];
        return ret;
    }

    Unit operator ^ (const int &rhs)
    {
        Unit ret;
        ret.first = std::pow(this->first,rhs);
        for (unsigned int i = 0; i < this->second.size(); i++) 
        {
            ret.second[i] = (this->second)[i]*rhs;
        }
        return ret;
    }

    friend std::ostream& operator<<(std::ostream &out, const Unit &a) 
    {
        std::cout << a.normalized_value() << " ";
        std::cout << a.normalized_unitstring();
        return out;
    }

    double normalized_value() const
    {
        double fac = 1.0;
        for (unsigned int i = 0 ; i < this->second.size(); i++)
            fac *= std::pow(normalization[i].second, this->second[i]);
        
        return this->first * fac;
    }

    std::string normalized_unitstring() const
    {
        // little utility to join strings with delimiter
        auto join = [] (std::vector<std::string>vec,std::string sep){
            std::ostringstream oss;
            for (size_t i = 0; i < vec.size(); ++i) {
                oss << vec[i];
                if (i < vec.size() - 1) {
                    oss << sep;
                }
            }
            return oss.str();
        };

        std::vector<std::string> num, den;
        for (unsigned int i = 0 ; i < this->second.size(); i++)
        {
            if (this->second[i] > 1)        
                num.push_back(normalization[i].first + "^" + std::to_string(this->second[i]));
            else if (this->second[i] == 1)  
                num.push_back(normalization[i].first);
            else if (this->second[i] == -1) 
                den.push_back(normalization[i].first);
            else if (this->second[i] < -1)  
                den.push_back(normalization[i].first + "^" + std::to_string(std::abs(this->second[i])));
        }
        
        std::string ret;
        if (num.size()) ret = join(num,"*");
        else ret = "1";
        if (den.size()) ret += "/" + join(den,"/");
        
        return ret;
    }
    





};


inline std::map<std::string, int> 
Unit::ParseUnitString(const std::string& input) {
    std::map<std::string, int> result;
    enum class Op { Mul = 1, Div = -1 };
    Op currentOp = Op::Mul;
    std::string token;
    size_t i = 0;
    auto parse_unit = [&](const std::string& t, Op op) {
        if (t.empty()) return;

        size_t caret = t.find('^');
        std::string base = t;
        int exp = 1;

        if (caret != std::string::npos) {
            base = t.substr(0, caret);
            std::string exp_str = t.substr(caret + 1);
            try {
                exp = std::stoi(exp_str);
            } catch (...) {
                throw std::runtime_error("Invalid exponent: " + exp_str);
            }
        }
        int signed_exp = (op == Op::Mul ? +1 : -1) * exp;
        result[base] += signed_exp;
    };

    while (i < input.size()) {
        char c = input[i];

        if (c == '*') 
        {
            parse_unit(token, currentOp);
            token.clear();
            currentOp = Op::Mul;
            ++i;
        } 
        else if (c == '/') 
        {
            parse_unit(token, currentOp);
            token.clear();
            currentOp = Op::Div;
            ++i;
        } 
        else 
        {
            token += c;
            ++i;
        }
    }

    // Parse final token
    parse_unit(token, currentOp);

    return result;
}










#endif
