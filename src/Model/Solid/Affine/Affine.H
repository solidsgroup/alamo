#ifndef MODEL_SOLID_AFFINE_H_
#define MODEL_SOLID_AFFINE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Affine
{
public:
    struct
    {
        Set::Matrix epsp;
        Set::Scalar alpha;
        Set::Matrix beta;
    }PlasticState;
};
}
}
}
#endif