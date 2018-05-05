#ifndef SET_SET_
#define SET_SET_

#include <eigen3/Eigen/Core>
#include <AMReX_REAL.H>

namespace Set
{
using Scalar = amrex::Real;
using Vector = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Matrix = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;
}
#endif
