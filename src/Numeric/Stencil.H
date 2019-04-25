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

static Type DefaultType[AMREX_SPACEDIM] = {AMREX_D_DECL(Type::Central, Type::Central, Type::Central)};

template<int x, int y, int z>
AMREX_FORCE_INLINE
Set::Scalar D(const amrex::Array4<const amrex::Real> &field,
	      const int &i, const int &j, const int &k, const int &m,
	      const Set::Scalar dx[AMREX_SPACEDIM],
	      Type stencil[AMREX_SPACEDIM] = DefaultType);

template<>
AMREX_FORCE_INLINE
Set::Scalar D<1,0,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM],
		     Type stencil[AMREX_SPACEDIM])
{
	if (stencil[0] == Type::Lo)
		return (f(i,j,k,m) - f(i-1,j,k,m)) / dx[0];
	else if (stencil[0] == Type::Hi)
		return (f(i+1,j,k,m) - f(i,j,k,m)) / dx[0];
	else
		return 0.5 * (f(i+1,j,k,m) - f(i-1,j,k,m)) / dx[0];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,1,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM],
		     Type stencil[AMREX_SPACEDIM])
{
	if (stencil[0] == Type::Lo)
		return (f(i,j,k,m) - f(i,j-1,k,m)) / dx[1];
	else if (stencil[0] == Type::Hi)
		return (f(i,j+1,k,m) - f(i,j,k,m)) / dx[1];
	else
		return 0.5 * (f(i,j+1,k,m) - f(i,j-1,k,m)) / dx[1];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,0,1>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM],
		     Type stencil[AMREX_SPACEDIM])
{
	if (stencil[0] == Type::Lo)
		return (f(i,j,k,m) - f(i,j,k-1,m)) / dx[2];
	else if (stencil[0] == Type::Hi)
		return (f(i,j,k+1,m) - f(i,j,k,m)) / dx[2];
	else
		return 0.5 * (f(i,j,k+1,m) - f(i,j,k-1,m)) / dx[2];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<2,0,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return (f(i+1,j,k,m) - 2.0*f(i,j,k,m) + f(i-1,j,k,m))/dx[0]/dx[0];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,2,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return (f(i,j+1,k,m) - 2.0*f(i,j,k,m) + f(i,j-1,k,m))/dx[1]/dx[1];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,0,2>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return (f(i,j,k+1,m) - 2.0*f(i,j,k,m) + f(i,j,k-1,m))/dx[2]/dx[2];
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<1,1,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return 0.25*(f(i+1,j+1,k,m) + f(i-1,j-1,k,m) - f(i+1,j-1,k,m) - f(i-1,j+1,k,m))/(dx[0]*dx[1]);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<1,0,1>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return 0.25*(f(i+1,j,k+1,m) + f(i-1,j,k-1,m) - f(i+1,j,k-1,m) - f(i-1,j,k+1,m))/(dx[0]*dx[2]);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,1,1>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return 0.25*(f(i,j+1,k+1,m) + f(i,j-1,k-1,m) - f(i,j+1,k-1,m) - f(i,j-1,k+1,m))/(dx[1]*dx[2]);
}


template<>
AMREX_FORCE_INLINE
Set::Scalar D<4,0,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return (+    (field(i+2,j,k,m))
		- 4.*(field(i+1,j,k,m))
		+ 6.*(field(i  ,j,k,m))
		- 4.*(field(i-1,j,k,m))
		+    (field(i-2,j,k,m))) /
		(dx[0]*dx[0]*dx[0]*dx[0]);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,4,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) /
		(dx[1]*dx[1]*dx[1]*dx[1]);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<3,1,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	// (2nd order in x, 4th order in y)
	return	((   -f(i+2,j+2,k,m) + 8.0*f(i+2,j+1,k,m) - 8.0*f(i+2,j-1,k,m) + f(i+2,j-2,k,m))
		 -2*(-f(i+1,j+2,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i+1,j-1,k,m) + f(i+1,j-2,k,m))
		 +2*(-f(i-1,j+2,k,m) + 8.0*f(i-1,j+1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-1,j-2,k,m))
		 -(  -f(i-2,j+2,k,m) + 8.0*f(i-2,j+1,k,m) - 8.0*f(i-2,j-1,k,m) + f(i-2,j-2,k,m))) /
		(24.0*dx[0]*dx[0]*dx[0]*dx[1]);
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
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return ((   -f(i+2,j+2,k,m) + 8.0*f(i+1,j+2,k,m) - 8.0*f(i-1,j+2,k,m) + f(i-2,j+2,k,m))
		-2*(-f(i+2,j+1,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i-1,j+1,k,m) + f(i-2,j+1,k,m))
		+2*(-f(i+2,j-1,k,m) + 8.0*f(i+1,j-1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-2,j-1,k,m))
		-(  -f(i+2,j-2,k,m) + 8.0*f(i+1,j-2,k,m) - 8.0*f(i-1,j-2,k,m) + f(i-2,j-2,k,m))) /
		(24.0*dx[0]*dx[1]*dx[1]*dx[1]) ;
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<2,2,0>(const amrex::Array4<const amrex::Real> &f,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar dx[AMREX_SPACEDIM], Type[])
{
	return 	(-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
		 +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
		 -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
		 +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
		 -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
		(144.0*dx[0]*dx[0]*dx[1]*dx[1]) ;
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
