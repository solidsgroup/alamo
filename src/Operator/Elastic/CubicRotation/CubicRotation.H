#ifndef OPERATOR_ELASTIC_CUBIC_ROTATION_H_
#define OPERATOR_ELASTIC_CUBIC_ROTATION_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "eigen3/Eigen/Core"

#include "Operator/Elastic/Elastic.H"

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
	/// \brief Constructor, using Bunge Euler angles (in degrees)
    CubicRotation (const amrex::Real, const amrex::Real, const amrex::Real,
		   amrex::Real, amrex::Real, amrex::Real);
    /// \brief Destructor (does nothing)
    virtual ~CubicRotation () {};

    void SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc);

  protected:

    amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		  const int amrlev, const int mglev, const MFIter &mfi) const;

    Eigen::Matrix<amrex::Real, 3, 3> R;
    amrex::Real Cc[3][3][3][3];
  };
}
}
#endif
