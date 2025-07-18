#ifndef UTIL_UNIT_H
#define UTIL_UNIT_H


#include <array>
#include <iostream>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <cmath>
#include <regex>

//#include "Set/Base.H"
//#include "String.H"

namespace Util
{
namespace Unit
{

static const std::map<std::string, std::pair<double,std::array<int,7>>> base_units =
{
    {"m", {1.0, {1,0,0,0,0,0,0}}},
    {"s", {1.0, {0,1,0,0,0,0,0}}},
    {"kg",{1.0, {0,0,1,0,0,0,0}}}
};

static const std::map<std::string, std::pair<double,std::string>> compound = {
    
    {"km", {1e+3, "m"}},
    {"cm", {1e-2, "m"}},
    {"mm", {1e-3, "m"}},
    {"um", {1e-6, "m"}},
    {"nm", {1e-9, "m"}},

    {"N",   {  1.0, "kg*m/s^2"}},
    {"kN",  {  1e3, "N"       }},
    {"Pa",  {  1.0, "N/m^2"}},
    {"kPa", {  1e3, "N/m^2"}},
    {"MPa", {  1e6, "N/m^2"}},
    {"GPa", { 1e9, "N/m^2"}},
};

using UnitMap = std::map<std::string, int>;
UnitMap ParseUnitString(const std::string& input);


struct Type : std::pair<double,std::array<int,7>>
{





    using std::pair<double, std::array<int,7>>::pair;
    // Explicit constructor from the base pair
    Type(const std::pair<double, std::array<int, 7>>& p) : std::pair<double, std::array<int, 7>>(p) {}


    static Type Parse(std::string unitstring)
    {
        Type type = Type::Unity();
        double factor = 1.0;
        UnitMap map = ParseUnitString(unitstring);
        for (auto it = map.begin(); it != map.end(); it++)
        {
            auto [str, exp] = *it;

            if (base_units.find(str) != base_units.end())
            {
                Type t = Util::Unit::base_units.at(str);
                type = type * (t^exp);
            }
            else if (compound.find(it->first) != compound.end())
            {
                auto [subfac, substr] = compound.at(str);
                Type t = Parse(substr);
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



    static Type Unity()  { return Type{1.0, {0,0,0,0,0,0,0}}; }
    static Type Length() { return Type{1.0, {1,0,0,0,0,0,0}}; }
    static Type Time()   { return Type{1.0, {0,1,0,0,0,0,0}}; }
    static Type Mass()   { return Type{1.0, {0,0,1,0,0,0,0}}; }

    static Type Area() { return Type{1.0, {2,0,0,0,0,0,0}}; }
    static Type Volume() { return Type{1.0, {3,0,0,0,0,0,0}}; }
    
    static Type Velocity() {return Length() / Time();}
    static Type Force() {return Mass()*Length() / Time()^2;}
    static Type Energy() {return Force()*Length();}
    static Type Pressure() {return Force()/Area();}


    void setLengthUnit(std::string unit)
    {
        Type parsed = Parse(unit);
        if (parsed.second != Length().second)
            std::cout << "not a unit of length" << std::endl;
        normalization[0].first = unit;
        normalization[0].second = 1.0 / parsed.first;
    }


    Type operator * (const Type &rhs)
    {
        Type ret;
        ret.first = this->first * rhs.first;
        for (int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] + rhs.second[i];
        return ret;
    }
    Type operator / (const Type &rhs)
    {
        Type ret;
        ret.first = this->first / rhs.first;
        for (int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] - rhs.second[i];
        return ret;
    }

    Type operator ^ (const int &rhs)
    {
        Type ret;
        ret.first = std::pow(this->first,rhs);
        for (int i = 0; i < this->second.size(); i++) 
        {
            ret.second[i] = (this->second)[i]*rhs;
        }
        return ret;
    }

    void print() 
    {
        std::cout << normalized_value() << " ";
        std::cout << normalized_unitstring() << "  [";
        for (int i = 0; i < this->second.size(); i++) std::cout << (this->second)[i] << " ";
        std::cout << "]";
        std::cout << std::endl;
    }

    double normalized_value()
    {
        double fac = 1.0;
        for (int i = 0 ; i < this->second.size(); i++)
            fac *= std::pow(normalization[i].second, this->second[i]);
        
        return this->first * fac;
    }

    std::string normalized_unitstring()
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
        for (int i = 0 ; i < this->second.size(); i++)
        {
            if (this->second[i] > 1)        
                num.push_back(normalization[i].first + "^" + std::to_string(this->second[i]));
            else if (this->second[i] == 1)  
                num.push_back(normalization[i].first);
            else if (this->second[i] == -1) 
                den.push_back(normalization[i].first);
            else if (this->second[i] < -1)  
                den.push_back(normalization[i].first + "^" + std::to_string(abs(this->second[i])));
        }
        
        std::string ret;
        if (num.size()) ret = join(num,"*");
        else ret = "1";
        if (den.size()) ret += "/" + join(den,"/");
        
        return ret;
    }
    

    inline static std::array<std::pair<std::string,double>, 7> normalization = {
        std::pair{"m",  1.0},
        std::pair{"s",  1.0},
        std::pair{"kg", 1.0},
        std::pair{"K",  1.0},
        std::pair{"A",  1.0},
        std::pair{"mol",1.0},
        std::pair{"cd", 1.0}
    };

};












UnitMap ParseUnitString(const std::string& input) {
    UnitMap result;
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


/*
//enum Type {
//    None,
//    Length, Time, Mass, Temperature, Current,
//    Velocity, Acceleration,
//    Force, Pressure, Energy, Power,
//    Viscosity, Density
//};x

const std::array<std::unordered_map<std::string, double>,3> map =
{
    // None
    std::unordered_map<std::string, double>{
    },
    // Length
    std::unordered_map<std::string, double>{
        {"m",   1.0},
        {"km",  1e3},
        {"cm",  1e-2},
        {"mm",  1e-3},
        {"um",  1e-6},  // micron (µm)
        {"nm",  1e-9}
    },
    // Time
    std::unordered_map<std::string, double>{
        {"s",   1.0},
        {"ms",  1e-3},
        {"us",  1e-6},  // microsecond
        {"ns",  1e-9}
    }
};


using Type = std::pair<double,std::string>;
std::unordered_map<std::string, Type> relations = {
    
    {"km", {1e+3, "m"}},
    {"cm", {1e-2, "m"}},
    {"mm", {1e-3, "m"}},
    {"um", {1e-6, "m"}},
    {"nm", {1e-9, "m"}},

    {"N",   {  1.0, "kg*m/s^2"}},
    {"kN",  {  1e3, "N"       }},
    {"Pa",  {  1.0, "N/m^2"}},
    {"kPa", {  1e3, "N/m^2"}},
    {"MPa", {  1e6, "N/m^2"}},
    {"GPa", { 1e9, "N/m^2"}},
};




std::unordered_map<std::string, double> length = {
    {"m",   1.0},
    {"km",  1e3},
    {"cm",  1e-2},
    {"mm",  1e-3},
    {"um",  1e-6},  // micron (µm)
    {"nm",  1e-9}
};
inline double convert_length(double val, std::string unit)
{
    return 0;
}


std::unordered_map<std::string, std::string> compound_units = {
    {"kph", "km/hour"},
};

inline double convert(UnitType, double val, std::string unit)
{
    if      (length.contains(unit)) return convert_length(val,unit);
    else if (time.contains(unit))   return convert_time(val,unit);
    return 0;
}


UnitType ParseUnitString(const std::string& input) {
    std::string unit_expr = input;

    // Expand compound units recursively
    while (compound_units.count(unit_expr)) {
        unit_expr = compound_units.at(unit_expr);
    }

    UnitType total = {};
    std::regex token_re(R"(([a-zA-Z]+)(\^(-?\d+))?)");
    std::regex div_re(R"(([^/]+)(?:/(.+))?)");

    std::smatch div_match;
    if (std::regex_match(unit_expr, div_match, div_re)) {
        std::string numerator = div_match[1];
        std::string denominator = div_match[2];

        auto parse_side = [&](const std::string& side, int sign) {
            std::sregex_iterator it(side.begin(), side.end(), token_re), end;
            for (; it != end; ++it) {
                std::string sym = (*it)[1];
                int power = (*it)[3].matched ? std::stoi((*it)[3]) : 1;
                UnitType u;

                // Recursive expansion if needed
                if (compound_units.count(sym)) {
                    u = ParseUnitString(compound_units.at(sym));
                } else if (base_units.count(sym)) {
                    u = base_units.at(sym);
                } else {
                    throw std::runtime_error("Unknown unit: " + sym);
                }

                total += u * (sign * power);
            }
        };

        parse_side(numerator, +1);
        if (denominator != "") parse_side(denominator, -1);
    }

    return total;
}


template <int T>
inline double Parse(const std::string& input);
*/


}
}


#endif
