#ifndef OPERATORCELL_ELASTIC_POLYCRYSTAL_CUBIC_H_
#define OPERATORCELL_ELASTIC_POLYCRYSTAL_CUBIC_H_

#include <limits>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>

#include "Operator/Elastic/PolyCrystal/PolyCrystal.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
namespace PolyCrystal
{
class Cubic : public PolyCrystalModel
{
public:

  /// \brief Constructor with randomized orientation
  Cubic (amrex::Real _C11, amrex::Real _C12, amrex::Real _C44);

  /// \brief Base constructor
  Cubic (amrex::Real _C11, amrex::Real _C12, amrex::Real _C44,
	 amrex::Real phi1, amrex::Real Phi, amrex::Real phi2);

protected:
  
  amrex::Real C(const int i, const int j, const int k, const int l) const;

  amrex::Real C11, C12, C44;
  
  Eigen::Matrix3d R;
};
}
}
}
#endif

