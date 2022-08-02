//
// "Affine" generally means "linear with an offset". Here we use "affine" to
// refer to models that are elastic with an eigenstrain, i.e.
//
// .. math::
//
//    \sigma = \mathbb{C}(\varepsilon - \varepsilon_0)
//
// The quantity :math:`\varepsilon_0` is any kind of eigenstrain - examples 
// include thermal strain, plastic strain, or viscous strain.
// This class can be used directly, where the eigenstrain is read in or
// set by the integrator.
// There are also classes (particularly visco/plastic models) that inherit
// from this type of model.
//
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
public:
    Set::Matrix F0;
    void SetF0(Set::Matrix a_F0) { F0 = a_F0; }
public:
    KinematicVariable kinvar = KinematicVariable::gradu;
};
}
}
}
#endif
