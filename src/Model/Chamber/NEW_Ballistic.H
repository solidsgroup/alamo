#ifndef MODEL_CHAMBER_BALLISTIC_H_
#define MODEL_CHAMBER_BALLISTIC_H_


// TODO: Imports
#include "Util/Util.H"
#include "IO/ParmParse.H"
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
    Ballistic(IO::ParmParse &pp, std::string name) { pp_queryclass(name, *this); }
    virtual ~Ballistic() = default;

public:

    // TODO: Write function to compute chamber values


    static void Parse (Ballistic &value, IO::ParmParse &pp)
    { 
        // TODO: Complete out query
        // Sample pp.query - Throat area [m^2]
        pp.query_default("At",value.At,"0.1_m^2", Unit::Area());
        
    }

};
}
}


#endif