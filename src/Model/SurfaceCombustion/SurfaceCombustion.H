#ifndef MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H
#define MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H

#include "IO/ParmParse.H"

namespace Model
{
namespace SurfaceCombustion
{
class SurfaceCombustion
{
public:
    virtual void set_pressure(const Set::Scalar){};

    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) = 0;

};
}
}


#endif


