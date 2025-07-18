#ifndef UTIL_UNIT_H
#define UTIL_UNIT_H


#include <array>
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

namespace Util
{

struct Unit : std::pair<double,std::array<int,7>>
{

    //
    // STATIC CONST DEFINITIONS, TRUE EVERYWHERE
    // 
    enum class Type {
        Length=0, Time=1, Mass=2, Temperature=3, Current=4, Amount=5, Intensity=6
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
    // Compound unit relationships and nomenclature
    inline static const std::map<std::string, std::pair<double,std::string>> compound = {
        // Lengths
        {"km", {1e+3,   "m"}},
        {"cm", {1e-2,   "m"}},
        {"mm", {1e-3,   "m"}},
        {"um", {1e-6,   "m"}},
        {"nm", {1e-9,   "m"}},
        {"ft", {0.3048, "m"}},
        {"mi", {5280., "ft"}},
        // Times
        {"ms", {1e-3, "s"}},
        {"us", {1e-6, "s"}},
        {"ns", {1e-9, "s"}},
        {"min", {60.0, "s"}},
        {"hr", {60.0, "min"}},
        // Speeds
        {"mph", {1.0, "mi/hr"}},
        // Force
        {"N",   {  1.0, "kg*m/s^2"}},
        {"kN",  {  1e3, "N"       }},
        // Pressure
        {"Pa",  {  1.0, "N/m^2"}},
        {"kPa", {  1e3, "N/m^2"}},
        {"MPa", {  1e6, "N/m^2"}},
        {"GPa", { 1e9, "N/m^2"}},
    };

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
            return StringParse(unit);
        }
        else
        {
            throw std::runtime_error("Invalid unit string "+unit);
        }
    }


    static Unit StringParse(std::string unitstring)
    {
        Unit type = Unit::Unity();
        double factor = 1.0;
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

    static Unit Unity()  { return Unit{1.0, {0,0,0,0,0,0,0}}; }
    static Unit Length() { return Unit{1.0, {1,0,0,0,0,0,0}}; }
    static Unit Time()   { return Unit{1.0, {0,1,0,0,0,0,0}}; }
    static Unit Mass()   { return Unit{1.0, {0,0,1,0,0,0,0}}; }

    static Unit Area() { return Unit{1.0, {2,0,0,0,0,0,0}}; }
    static Unit Volume() { return Unit{1.0, {3,0,0,0,0,0,0}}; }
    
    static Unit Velocity() {return Length() / Time();}
    static Unit Force() {return Mass()*Length() / Time()^2;}
    static Unit Energy() {return Force()*Length();}
    static Unit Pressure() {return Force()/Area();}

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
        for (int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] + rhs.second[i];
        return ret;
    }
    Unit operator / (const Unit &rhs)
    {
        Unit ret;
        ret.first = this->first / rhs.first;
        for (int i = 0; i < this->second.size(); i++) ret.second[i] = (this->second)[i] - rhs.second[i];
        return ret;
    }

    Unit operator ^ (const int &rhs)
    {
        Unit ret;
        ret.first = std::pow(this->first,rhs);
        for (int i = 0; i < this->second.size(); i++) 
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
        for (int i = 0 ; i < this->second.size(); i++)
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






}



#endif
