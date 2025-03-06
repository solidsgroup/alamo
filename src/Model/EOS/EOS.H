#ifndef MODEL_EOS_H_
#define MODEL_EOS_H_

//#include <AMReX.H>
//#include <AMReX_REAL.H>

#include "Set/Set.H"

namespace Model
{
namespace EOS
{
class EOS
{
public:
    EOS() {};

    virtual ~EOS() {};
    
    //Inputs should be density, internal energy (per unit mass), species composition, [additional parameters for real gas]
    virtual Set::Scalar pressure(const Set::Scalar density, const Set::Scalar energy) const
    {Util::Abort(INFO, "sdlf"):; return 0.0;};
    virtual Set::Scalar pressure(const Set::Scalar density, const Set::Scalar energy, const std::vector<Set::Scalar> Yfraction) const;
    {Util::Abort(INFO, "sdlf"):; return 0.0;};

    virtual void Advance(Set::Scalar /*dt*/, Set::Scalar /*density*/, Set::Scalar /*energy*/, std::vector<Set::Scalar> /*mass fractions*/) {}
};

}
}
    
#endif
