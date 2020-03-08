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
//class Affine
//{
//public:
    struct PlasticState{
        Set::Matrix epsp;   // Plastic strain
        Set::Scalar alpha;  // Equivalent plastic strain
        Set::Matrix beta;
    };
    /*struct{
        bool on;        // slip system active or not
        Set::Vector n;  // slip plane normal
        Set::Vector s;  // slip direction
        Set::Scalar tcrss;      // critically resolved shear stress
        Set::Scalar tcrss0;     // pristine material
        Set::Scalar exp;        // power law exponent
        Set::Scalar gam0dot;    // gamma0 dot.
    } slipsystem;*/
//};
}
}
}
#endif