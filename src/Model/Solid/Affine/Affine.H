#ifndef MODEL_SOLID_AFFINE_H_
#define MODEL_SOLID_AFFINE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Affine
{

template<Set::Sym SYM>
class Affine : public Solid<SYM>
{
protected:
    Set::Matrix F0;
    void SetF0(Set::Matrix a_F0) { F0 = a_F0; }
public:
    KinematicVariable kinvar = KinematicVariable::gradu;
};
}
}
}
#endif
