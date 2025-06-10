#ifndef MODEL_SURFACECOMBUSTION_HOMOGENIZE_H
#define MODEL_SURFACECOMBUSTION_HOMOGENIZE_H

#include "Model/SurfaceCombustion/SurfaceCombustion.H"
#include "SurfaceCombustion.H"

namespace Model
{
namespace SurfaceCombustion
{
class Homogenize
{
public:
    static constexpr const char* name = "homogenize";

    Homogenize() = default;
    Homogenize(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }


    void set_pressure(const Set::Scalar P) 
    {
        // TODO
    }
    
    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) 
    {
        // TODO
    }

    static void Parse(Homogenize& value, IO::ParmParse& pp)
    {
        // TODO
    }

    // TODO
};

}
}



#endif


