#ifndef NUMERIC_STENCIL_H_
#define NUMERIC_STENCIL_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"


namespace Numeric
{

enum StencilType {Lo, Hi, Central};
static StencilType DefaultType[AMREX_SPACEDIM] = {AMREX_D_DECL(StencilType::Central, StencilType::Central, StencilType::Central)};

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
		   StencilType stencil[AMREX_SPACEDIM] = DefaultType)
	{
		if (stencil[0] == StencilType::Lo)
			return (f(i,j,k,m) - f(i-1,j,k,m)) / dx[0];
		else if (stencil[0] == StencilType::Hi)
			return (f(i+1,j,k,m) - f(i,j,k,m)) / dx[0];
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
		   StencilType stencil[AMREX_SPACEDIM] = DefaultType)
	{
		if (stencil[1] == StencilType::Lo)
			return (f(i,j,k,m) - f(i,j-1,k,m)) / dx[1];
		else if (stencil[1] == StencilType::Hi)
			return (f(i,j+1,k,m) - f(i,j,k,m)) / dx[1];
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
		   StencilType stencil[AMREX_SPACEDIM] = DefaultType)
	{
		if (stencil[2] == StencilType::Lo)
			return (f(i,j,k,m) - f(i,j,k-1,m)) / dx[2];
		else if (stencil[2] == StencilType::Hi)
			return (f(i,j,k+1,m) - f(i,j,k,m)) / dx[2];
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
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
		   const Set::Scalar dx[AMREX_SPACEDIM],
		   StencilType [AMREX_SPACEDIM] = DefaultType)
	{
		return 	(-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
			 +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
			 -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
			 +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
			 -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
			(144.0*dx[0]*dx[0]*dx[1]*dx[1]) ;
	};
};

//
// OLD
//


// template<class T>
// AMREX_FORCE_INLINE
// T D100(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM],
// 		 Type stencil[AMREX_SPACEDIM] = DefaultType)
// {
// 	if (stencil[0] == Type::Lo)
// 		return (f(i,j,k,m) - f(i-1,j,k,m)) / dx[0];
// 	else if (stencil[0] == Type::Hi)
// 		return (f(i+1,j,k,m) - f(i,j,k,m)) / dx[0];
// 	else
// 		return (f(i+1,j,k,m) - f(i-1,j,k,m))*0.5 / dx[0];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D010(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM],
// 		 Type stencil[AMREX_SPACEDIM] = DefaultType)
// {
// 	if (stencil[0] == Type::Lo)
// 		return (f(i,j,k,m) - f(i,j-1,k,m)) / dx[1];
// 	else if (stencil[0] == Type::Hi)
// 		return (f(i,j+1,k,m) - f(i,j,k,m)) / dx[1];
// 	else
// 		return (f(i,j+1,k,m) - f(i,j-1,k,m))*0.5 / dx[1];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D001(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM],
// 		 Type stencil[AMREX_SPACEDIM] = DefaultType)
// {
// 	if (stencil[0] == Type::Lo)
// 		return (f(i,j,k,m) - f(i,j,k-1,m)) / dx[2];
// 	else if (stencil[0] == Type::Hi)
// 		return (f(i,j,k+1,m) - f(i,j,k,m)) / dx[2];
// 	else
// 		return (f(i,j,k+1,m) - f(i,j,k-1,m))*0.5 / dx[2];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D200(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return (f(i+1,j,k,m) - 2.0*f(i,j,k,m) + f(i-1,j,k,m))/dx[0]/dx[0];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D020(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return (f(i,j+1,k,m) - 2.0*f(i,j,k,m) + f(i,j-1,k,m))/dx[1]/dx[1];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D002(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return (f(i,j,k+1,m) - 2.0*f(i,j,k,m) + f(i,j,k-1,m))/dx[2]/dx[2];
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D110(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return 0.25*(f(i+1,j+1,k,m) + f(i-1,j-1,k,m) - f(i+1,j-1,k,m) - f(i-1,j+1,k,m))/(dx[0]*dx[1]);
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D101(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return 0.25*(f(i+1,j,k+1,m) + f(i-1,j,k-1,m) - f(i+1,j,k-1,m) - f(i-1,j,k+1,m))/(dx[0]*dx[2]);
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D011(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return 0.25*(f(i,j+1,k+1,m) + f(i,j-1,k-1,m) - f(i,j+1,k-1,m) - f(i,j-1,k+1,m))/(dx[1]*dx[2]);
// }


// template<class T>
// AMREX_FORCE_INLINE
// T D400(const amrex::Array4<const T> &field,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return (+    (field(i+2,j,k,m))
// 		- 4.*(field(i+1,j,k,m))
// 		+ 6.*(field(i  ,j,k,m))
// 		- 4.*(field(i-1,j,k,m))
// 		+    (field(i-2,j,k,m))) /
// 		(dx[0]*dx[0]*dx[0]*dx[0]);
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D040(const amrex::Array4<const T> &field,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) /
// 		(dx[1]*dx[1]*dx[1]*dx[1]);
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D310(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	// (2nd order in x, 4th order in y)
// 	return	((   -f(i+2,j+2,k,m) + 8.0*f(i+2,j+1,k,m) - 8.0*f(i+2,j-1,k,m) + f(i+2,j-2,k,m))
// 		 -2*(-f(i+1,j+2,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i+1,j-1,k,m) + f(i+1,j-2,k,m))
// 		 +2*(-f(i-1,j+2,k,m) + 8.0*f(i-1,j+1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-1,j-2,k,m))
// 		 -(  -f(i-2,j+2,k,m) + 8.0*f(i-2,j+1,k,m) - 8.0*f(i-2,j-1,k,m) + f(i-2,j-2,k,m))) /
// 		(24.0*dx[0]*dx[0]*dx[0]*dx[1]);
// 	// (2nd order in x, 2nd order in y)
// 	// return 
// 	// 	0.5*( - (- 0.5*f(i-2,j-1,k,m) + f(i-1,j-1,k,m) - f(i+1,j-1,k,m) + 0.5*f(i+2,j-1,k,m))
// 	// 	      + (- 0.5*f(i-2,j+1,k,m) + f(i-1,j+1,k,m) - f(i+1,j+1,k,m) + 0.5*f(i+2,j+1,k,m)) ) /
// 	// 	(dx*dx*dx*dy);
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D130(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return ((   -f(i+2,j+2,k,m) + 8.0*f(i+1,j+2,k,m) - 8.0*f(i-1,j+2,k,m) + f(i-2,j+2,k,m))
// 		-2*(-f(i+2,j+1,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i-1,j+1,k,m) + f(i-2,j+1,k,m))
// 		+2*(-f(i+2,j-1,k,m) + 8.0*f(i+1,j-1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-2,j-1,k,m))
// 		-(  -f(i+2,j-2,k,m) + 8.0*f(i+1,j-2,k,m) - 8.0*f(i-1,j-2,k,m) + f(i-2,j-2,k,m))) /
// 		(24.0*dx[0]*dx[1]*dx[1]*dx[1]) ;
// }

// template<class T>
// AMREX_FORCE_INLINE
// T D220(const amrex::Array4<const T> &f,
// 		 const int &i, const int &j, const int &k, const int &m,
// 		 const Set::Scalar dx[AMREX_SPACEDIM], Type[] = DefaultType)
// {
// 	return 	(-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
// 		 +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
// 		 -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
// 		 +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
// 		 -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
// 		(144.0*dx[0]*dx[0]*dx[1]*dx[1]) ;
// }


// template<class T>
// AMREX_FORCE_INLINE
// Set::Scalar D<0,4,0>(const amrex::Array4<const T> &field,
// 		     const int &i, const int &j, const int &k, const int &m,
// 		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &dz,
// 		     Type tx = Type::Central, Type ty = Type::Central, Type tz = Type::Central)
// {
// 	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) / (dy*dy*dy*dy);
// }



// }
}
#endif
