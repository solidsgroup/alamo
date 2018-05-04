#ifndef OPERATOR_ELASTIC_CUBIC_ROTATION_H_
#define OPERATOR_ELASTIC_CUBIC_ROTATION_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "eigen3/Eigen/Core"

#include "Operator/Elastic/Elastic.H"

///
/// \class CubicRotation
/// \brief linop class that creates a rotated elastic modulus tensor \f$\mathbb{C}_{ijkl}\f$
///
/// ## Supported constructor inputs ##
///
/// - 3x3 rotation matrix, as eigen type
/// - Bunge Euler angles
///
/// Use examples: `code`
/// ```cpp
/// Operator::Elastic::CubicRotation mlabec(0.9, 0.2, 0.1, C11, C12, C44);  //Bunge Euler angles in radians
/// // OR
/// Eigen::Matrix<amrex::Real, 3, 3> R;
///  R << 0.8755949, -0.3817528,  0.2959702,
///       0.4200312,  0.9043038, -0.0762129,
///       -0.2385525,  0.1910484, 0.9521519;  //3d rotation matrix for 30 deg rotation angle about [1,2,3]
/// Operator::Elastic::CubicRotation mlabec(R, C11, C12, C44);
/// ```

using namespace amrex;

namespace Operator
{
namespace Elastic
{
  class CubicRotation : public Elastic
  {
  public:

    /// \brief Base constructor
    CubicRotation ();
    /// \brief Typical constructor from Cubic
    CubicRotation (const Vector<Geometry>& a_geom,
		   const Vector<BoxArray>& a_grids,
		   const Vector<DistributionMapping>& a_dmap,
		   const LPInfo& a_info);
    /// \brief Constructor, using rotation matrix
    CubicRotation (const Eigen::Matrix<amrex::Real, 3, 3>,
	       amrex::Real, amrex::Real, amrex::Real);
	/// \brief Constructor, using Bunge Euler angles (in radians)
    CubicRotation (const amrex::Real, const amrex::Real, const amrex::Real,
		   amrex::Real, amrex::Real, amrex::Real);
    /// \brief Destructor (does nothing)
    virtual ~CubicRotation () {};

  protected:

    amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		  const int amrlev, const int mglev, const MFIter &mfi) const;

    Eigen::Matrix<amrex::Real, 3, 3> R;
    amrex::Real Cc[3][3][3][3];
  };
}
}
#endif
