#ifndef SET_BASE_H
#define SET_BASE_H

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_Vector.H"
#include "AMReX_BLassert.H"
#include "AMReX_Geometry.H"
#include "AMReX_IndexType.H"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/SVD>
#include <AMReX_REAL.H>
#include "Util/Util.H"

namespace Set
{
using Scalar   = amrex::Real;
using Vector   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;
using iMatrix  = Eigen::Matrix<int,AMREX_SPACEDIM,AMREX_SPACEDIM>;

AMREX_FORCE_INLINE
Vector Position(const int &i, const int &j, const int &k, const amrex::Geometry &geom, const amrex::IndexType &ixType)
{
    (void)k;
    if (ixType == amrex::IndexType::TheNodeType())
    {
        return Vector(AMREX_D_DECL(
                geom.ProbLo()[0] + ((amrex::Real)(i)) * geom.CellSize()[0],
                geom.ProbLo()[1] + ((amrex::Real)(j)) * geom.CellSize()[1],
                geom.ProbLo()[2] + ((amrex::Real)(k)) * geom.CellSize()[2]));
    }
    else if (ixType == amrex::IndexType::TheCellType())
    {
        return Vector(AMREX_D_DECL(
                geom.ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom.CellSize()[0],
                geom.ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom.CellSize()[1],
                geom.ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom.CellSize()[2]));
    }
    Util::Abort(INFO,"Unsupported index type");
    return Set::Vector::Zero(); 
}

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
