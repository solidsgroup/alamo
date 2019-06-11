#ifndef NUMERIC_STENCIL_H_
#define NUMERIC_STENCIL_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"

/// \brief This namespace contains some numerical tools
namespace Numeric
{

enum StencilType {Lo, Hi, Central};
static std::array<StencilType,AMREX_SPACEDIM>
DefaultType = {AMREX_D_DECL(StencilType::Central, StencilType::Central, StencilType::Central)};

static
AMREX_FORCE_INLINE
std::array<StencilType,AMREX_SPACEDIM>
GetStencil(const int i, const int j, const int k, const amrex::Box domain)
{
	(void)i; (void)j; (void)k; // Suppress "unused variable" warnings.
	std::array<StencilType,AMREX_SPACEDIM> sten;
	const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
	AMREX_D_TERM(sten[0] = (i == lo.x ? Numeric::StencilType::Hi :
				i == hi.x ? Numeric::StencilType::Lo :
				Numeric::StencilType::Central);,
		     sten[1] = (j == lo.y ? Numeric::StencilType::Hi :
				j == hi.y ? Numeric::StencilType::Lo :
				Numeric::StencilType::Central);,
		     sten[2] = (k == lo.z ? Numeric::StencilType::Hi :
				k == hi.z ? Numeric::StencilType::Lo :
				Numeric::StencilType::Central););
	return sten;
}

template<class T,int x, int y, int z>
struct Stencil
{};

template<class T>
struct Stencil<T,1,0,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
	{
		if (stencil[0] == StencilType::Lo)
			//return (f(i,j,k,m) - f(i-1,j,k,m)) / dx[0];
			return (f(i,j,k,m)*(1.5) + f(i-1,j,k,m)*(-2.0) + f(i-2,j,k,m)*(0.5)) / dx[0]; // 2nd order stencil
		else if (stencil[0] == StencilType::Hi)
			//return (f(i+1,j,k,m) - f(i,j,k,m)) / dx[0]; // 1st order stencil
			return (f(i,j,k,m)*(-1.5) + f(i+1,j,k,m)*(2.0) + f(i+2,j,k,m)*(-0.5)) / dx[0]; // 2nd order stencil
		else
			return (f(i+1,j,k,m) - f(i-1,j,k,m))*0.5 / dx[0];
	};
};

template<class T>
struct Stencil<T,0,1,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
	{
		if (stencil[1] == StencilType::Lo)
			//return (f(i,j,k,m) - f(i,j-1,k,m)) / dx[1];
			return (f(i,j,k,m)*( 1.5) + f(i,j-1,k,m)*(-2.0) + f(i,j-2,k,m)*(0.5)) / dx[1]; // 2nd order stencil
		else if (stencil[1] == StencilType::Hi)
			//return (f(i,j+1,k,m) - f(i,j,k,m)) / dx[1];
			return (f(i,j,k,m)*(-1.5) + f(i,j+1,k,m)*(2.0) + f(i,j+2,k,m)*(-0.5)) / dx[1]; // 2nd order stencil
		else
			return (f(i,j+1,k,m) - f(i,j-1,k,m))*0.5 / dx[1];
	};
};

template<class T>
struct Stencil<T,0,0,1>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
	{
		if (stencil[2] == StencilType::Lo)
			//return (f(i,j,k,m) - f(i,j,k-1,m)) / dx[2];
			return (f(i,j,k,m)*( 1.5) + f(i,j,k-1,m)*(-2.0) + f(i,j,k-2,m)*(0.5)) / dx[2]; // 2nd order stencil
		else if (stencil[2] == StencilType::Hi)
			//return (f(i,j,k+1,m) - f(i,j,k,m)) / dx[2];
			return (f(i,j,k,m)*(-1.5) + f(i,j,k+1,m)*( 2.0) + f(i,j,k+2,m)*(-0.5)) / dx[2]; // 2nd order stencil
		else
			return (f(i,j,k+1,m) - f(i,j,k-1,m))*0.5 / dx[2];
	};
};

template<class T>
struct Stencil<T,2,0,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return (f(i+1,j,k,m) - 2.0*f(i,j,k,m) + f(i-1,j,k,m))/dx[0]/dx[0];
	};
};

template<class T>
struct Stencil<T,0,2,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return (f(i,j+1,k,m) - 2.0*f(i,j,k,m) + f(i,j-1,k,m))/dx[1]/dx[1];
	};
};

template<class T>
struct Stencil<T,0,0,2>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return (f(i,j,k+1,m) - 2.0*f(i,j,k,m) + f(i,j,k-1,m))/dx[2]/dx[2];
	};
};

template<class T>
struct Stencil<T,1,1,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return 0.25*(f(i+1,j+1,k,m) + f(i-1,j-1,k,m) - f(i+1,j-1,k,m) - f(i-1,j+1,k,m))/(dx[0]*dx[1]);
	};
};
template<class T>
struct Stencil<T,1,0,1>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return 0.25*(f(i+1,j,k+1,m) + f(i-1,j,k-1,m) - f(i+1,j,k-1,m) - f(i-1,j,k+1,m))/(dx[0]*dx[2]);
	};
};
template<class T>
struct Stencil<T,0,1,1>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return 0.25*(f(i,j+1,k+1,m) + f(i,j-1,k-1,m) - f(i,j+1,k-1,m) - f(i,j-1,k+1,m))/(dx[1]*dx[2]);
	};
};
template<class T>
struct Stencil<T,4,0,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return (+    (f(i+2,j,k,m))
			- 4.*(f(i+1,j,k,m))
			+ 6.*(f(i  ,j,k,m))
			- 4.*(f(i-1,j,k,m))
			+    (f(i-2,j,k,m))) /
			(dx[0]*dx[0]*dx[0]*dx[0]);
	};
};
template<class T>
struct Stencil<T,0,4,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return ((f(i,j+2,k,m)) - 4.*(f(i,j+1,k,m)) + 6.*(f(i,j,k,m)) - 4.*(f(i,j-1,k,m)) + (f(i,j-2,k,m))) /
			(dx[1]*dx[1]*dx[1]*dx[1]);
	};
};
template<class T>
struct Stencil<T,3,1,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
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
	};
};
template<class T>
struct Stencil<T,1,3,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return ((   -f(i+2,j+2,k,m) + 8.0*f(i+1,j+2,k,m) - 8.0*f(i-1,j+2,k,m) + f(i-2,j+2,k,m))
			-2*(-f(i+2,j+1,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i-1,j+1,k,m) + f(i-2,j+1,k,m))
			+2*(-f(i+2,j-1,k,m) + 8.0*f(i+1,j-1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-2,j-1,k,m))
			-(  -f(i+2,j-2,k,m) + 8.0*f(i+1,j-2,k,m) - 8.0*f(i-1,j-2,k,m) + f(i-2,j-2,k,m))) /
			(24.0*dx[0]*dx[1]*dx[1]*dx[1]) ;
	};
};
template<class T>
struct Stencil<T,2,2,0>
{
	AMREX_FORCE_INLINE
	static T D(const amrex::Array4<const T> &f,
		   const int &i, const int &j, const int &k, const int &m,
		   const Set::Scalar dx[AMREX_SPACEDIM])
	{
		return 	(-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
			 +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
			 -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
			 +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
			 -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
			(144.0*dx[0]*dx[0]*dx[1]*dx[1]) ;
	};
};
}
#endif
