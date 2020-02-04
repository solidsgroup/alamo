#ifndef SET_BASE_H
#define SET_BASE_H

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_Vector.H"
#include "AMReX_BLassert.H"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/SVD>
#include <AMReX_REAL.H>

namespace Set
{
using Scalar   = amrex::Real;
using Vector   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;
using iMatrix  = Eigen::Matrix<int,AMREX_SPACEDIM,AMREX_SPACEDIM>;

enum Sym {None,Major,Minor,MajorMinor,Diagonal,Full,Isotropic};
template<int dim,int sym> class Matrix4{};
template<int dim, int sym>
std::ostream&
operator<< (std::ostream& os, const Matrix4<dim,sym>& b)
{
	Matrix4<dim,sym> bcopy = b;
   	bcopy.Print(os);
	return os;
}
}
#endif