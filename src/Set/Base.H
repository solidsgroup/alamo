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
using Vector3d = Eigen::Vector3d;
using MultiVector = Eigen::Matrix<amrex::Real, Eigen::Dynamic, 1>; //Vectors with N-components
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;
using Matrix3d = Eigen::Matrix3d;
using MultiMatrix = Eigen::Matrix<amrex::Real, Eigen::Dynamic, Eigen::Dynamic>; //Vectors with N-components
using iMatrix  = Eigen::Matrix<int,AMREX_SPACEDIM,AMREX_SPACEDIM>;

AMREX_FORCE_INLINE
Set::Matrix reduce(const Set::Matrix3d &in)
{
    return in.block<AMREX_SPACEDIM,AMREX_SPACEDIM>(0,0);
}

AMREX_FORCE_INLINE
Set::Matrix3d expand(const Set::Matrix &in)
{
    Set::Matrix3d ret = Set::Matrix3d::Zero();
    ret.block<AMREX_SPACEDIM,AMREX_SPACEDIM>(0,0) = in;
    return ret;
}

//using Quaternion = Eigen::Quaterniond;

class Quaternion : public Eigen::Quaterniond
{
public:
    Quaternion()
        : Eigen::Quaterniond() {}
    Quaternion(Set::Scalar w, Set::Scalar x, Set::Scalar y, Set::Scalar z)
        : Eigen::Quaterniond(w,x,y,z) {}
    Quaternion(const Eigen::Matrix3d &R)
        : Eigen::Quaterniond(R) {}
    Quaternion(Eigen::Quaterniond &q)
        : Eigen::Quaterniond(q.w(),q.x(),q.y(),q.z()) {}

    AMREX_FORCE_INLINE
    void operator = (const Eigen::Quaterniond &rhs)
    {
        w() = rhs.w(); x() = rhs.x(); y() = rhs.y(); z() = rhs.z();
    }

    AMREX_FORCE_INLINE
    void operator += (const Quaternion &rhs)
    {
        w() += rhs.w(); x() += rhs.x(); y() += rhs.y(); z() += rhs.z();
    }
    // AMREX_FORCE_INLINE
    // void operator * (const Set::Scalar alpha)
    // {
    //     w += alpha; x *= alpha; y *= alpha; z *= alpha;
    // }
    friend Quaternion operator * (const Set::Scalar alpha, const Quaternion b);
    friend Quaternion operator * (const Quaternion b,const Set::Scalar alpha);
    friend Quaternion operator + (const Quaternion a, const Quaternion b);
    friend Quaternion operator - (const Quaternion a, const Quaternion b);
    friend bool operator == (const Quaternion a, const Quaternion b);
};
AMREX_FORCE_INLINE
Quaternion operator * (const Set::Scalar alpha, const Quaternion b)
{
    Quaternion ret;
    ret.w() = b.w()*alpha; ret.x() = b.x()*alpha; ret.y() = b.y()*alpha; ret.z() = b.z()*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
Quaternion operator * (const Quaternion b, const Set::Scalar alpha)
{
    Quaternion ret;
    ret.w() = b.w()*alpha; ret.x() = b.x()*alpha; ret.y() = b.y()*alpha; ret.z() = b.z()*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
Quaternion operator + (const Quaternion a, const Quaternion b)
{
    Quaternion ret;
    ret.w() = a.w()+b.w(); ret.x() = a.x()+b.x(); ret.y() = a.y()+b.y(); ret.z() = a.z()+b.z();
    return ret;
}    

AMREX_FORCE_INLINE
Quaternion operator - (const Quaternion a, const Quaternion b)
{
    Quaternion ret;
    ret.w() = a.w()-b.w(); ret.x() = a.x()-b.x(); ret.y() = a.y()-b.y(); ret.z() = a.z()-b.z();
    return ret;
}    

AMREX_FORCE_INLINE
bool operator == (const Quaternion a, const Quaternion b)
{
    if (fabs((a.w() - b.w())/(a.w()+a.w())) > 1E-8) return false;
    if (fabs((a.x() - b.x())/(a.x()+a.x())) > 1E-8) return false;
    if (fabs((a.y() - b.y())/(a.y()+a.y())) > 1E-8) return false;
    if (fabs((a.z() - b.z())/(a.z()+a.z())) > 1E-8) return false;
    return true;
}


static Scalar Garbage = NAN;

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

AMREX_FORCE_INLINE
Vector Normal(  AMREX_D_DECL(bool xmin,bool ymin, bool zmin),
                AMREX_D_DECL(bool xmax,bool ymax, bool zmax))
{
    Set::Vector N =  Set::Vector::Zero();
    if (xmin)      N(0) = -1.0;
    else if (xmax) N(0) =  1.0;
#if AMREX_SPACEDIM>1
    if (ymin)      N(1) = -1.0;
    else if (ymax) N(1) =  1.0;
#endif
#if AMREX_SPACEDIM>2
    if (zmin)      N(2) = -1.0;
    else if (zmax) N(2) =  1.0;
#endif 
    return N;
}

AMREX_FORCE_INLINE
Vector Size(const amrex::Geometry &geom)
{
    Set::Vector size;
    size(0) = geom.ProbHi()[0] - geom.ProbLo()[0];
#if AMREX_SPACEDIM>1
    size(1) = geom.ProbHi()[1] - geom.ProbLo()[1];
#endif
#if AMREX_SPACEDIM>2
    size(2) = geom.ProbHi()[2] - geom.ProbLo()[2];
#endif
    return size;
}


AMREX_FORCE_INLINE
Vector Volume(const int &i, const int &j, const int &k, const amrex::Geometry &geom, const amrex::IndexType &ixType)
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
