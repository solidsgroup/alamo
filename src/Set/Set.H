#ifndef SET_SET_
#define SET_SET_

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_Vector.H"
#include "AMReX_BLassert.H"
//#define EIGEN_MATRIX_PLUGIN "Set/Eigen_Amrex.H"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/SVD>
#include <AMReX_REAL.H>

/// \brief A collection of data types and simple data structures
namespace Set
{
using Scalar   = amrex::Real;
using Vector   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;

namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
}
}

namespace Util
{
    Set::Scalar Random();
}
#endif
