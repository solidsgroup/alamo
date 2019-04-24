#ifndef NUMERIC_STENCIL_H_
#define NUMERIC_STENCIL_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"

namespace Numeric
{
namespace Stencil
{

enum Type {Lo, Hi, Central};

template<int x, int y, int z>
AMREX_FORCE_INLINE
Set::Scalar D(const amrex::Array4<const amrex::Real> &field,
	      const int &i, const int &j, const int &k, const int &m,
	      const Set::Scalar &dx,const Set::Scalar &dy,const Set::Scalar &dz,
	      Type tx = Type::Central, Type ty = Type::Central, Type tz = Type::Central);

template<>
AMREX_FORCE_INLINE
Set::Scalar D<4,0,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &, const Set::Scalar &,
		     Type, Type, Type)
{
	return ((field(i+2,j,k,m)) - 4.*(field(i+1,j,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i-1,j,k,m)) + (field(i-2,j,k,m))) / (dx*dx*dx*dx);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,4,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &, const Set::Scalar &dy, const Set::Scalar &,
		     Type, Type, Type)
{
	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) / (dy*dy*dy*dy);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<3,1,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &,
		     Type, Type, Type)
{
	// (2nd order in x, 4th order in y)
	return	((   -f(i+2,j+2,k,m) + 8.0*f(i+2,j+1,k,m) - 8.0*f(i+2,j-1,k,m) + f(i+2,j-2,k,m))
		 -2*(-f(i+1,j+2,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i+1,j-1,k,m) + f(i+1,j-2,k,m))
		 +2*(-f(i-1,j+2,k,m) + 8.0*f(i-1,j+1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-1,j-2,k,m))
		 -(  -f(i-2,j+2,k,m) + 8.0*f(i-2,j+1,k,m) - 8.0*f(i-2,j-1,k,m) + f(i-2,j-2,k,m))) /
		(24.0*dx*dx*dx*dy);
	// (2nd order in x, 2nd order in y)
	// return 
	// 	0.5*( - (- 0.5*f(i-2,j-1,k,m) + f(i-1,j-1,k,m) - f(i+1,j-1,k,m) + 0.5*f(i+2,j-1,k,m))
	// 	      + (- 0.5*f(i-2,j+1,k,m) + f(i-1,j+1,k,m) - f(i+1,j+1,k,m) + 0.5*f(i+2,j+1,k,m)) ) /
	// 	(dx*dx*dx*dy);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<1,3,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &,
		     Type, Type, Type)
{
	return ((   -f(i+2,j+2,k,m) + 8.0*f(i+1,j+2,k,m) - 8.0*f(i-1,j+2,k,m) + f(i-2,j+2,k,m))
		-2*(-f(i+2,j+1,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i-1,j+1,k,m) + f(i-2,j+1,k,m))
		+2*(-f(i+2,j-1,k,m) + 8.0*f(i+1,j-1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-2,j-1,k,m))
		-(  -f(i+2,j-2,k,m) + 8.0*f(i+1,j-2,k,m) - 8.0*f(i-1,j-2,k,m) + f(i-2,j-2,k,m))) /
		(24.0*dx*dy*dy*dy) ;
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<2,2,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &,
		     Type, Type, Type)
{
	return 	(-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
		 +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
		 -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
		 +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
		 -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
		(144.0*dx*dx*dy*dy) ;
}


// template<>
// AMREX_FORCE_INLINE
// Set::Scalar D<0,4,0>(const amrex::Array4<const amrex::Real> &field,
// 		     const int &i, const int &j, const int &k, const int &m,
// 		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &dz,
// 		     Type tx = Type::Central, Type ty = Type::Central, Type tz = Type::Central)
// {
// 	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) / (dy*dy*dy*dy);
// }



}
}
#endif
