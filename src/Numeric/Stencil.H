#ifndef NUMERIC_STENCIL_H_
#define NUMERIC_STENCIL_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"

namespace Numeric
{
namespace Stencil
{

template<int x, int y, int z>
AMREX_FORCE_INLINE
Set::Scalar D(const amrex::Array4<const amrex::Real> &field,
	      const int &i, const int &j, const int &k, const int &m,
	      const Set::Scalar &dx,const Set::Scalar &dy,const Set::Scalar &dz);

template<>
AMREX_FORCE_INLINE
Set::Scalar D<4,0,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &dz)
{
	return ((field(i+2,j,k,m)) - 4.*(field(i+1,j,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i-1,j,k,m)) + (field(i-2,j,k,m))) / (dx*dx*dx*dx);
}

template<>
AMREX_FORCE_INLINE
Set::Scalar D<0,4,0>(const amrex::Array4<const amrex::Real> &field,
		     const int &i, const int &j, const int &k, const int &m,
		     const Set::Scalar &dx, const Set::Scalar &dy, const Set::Scalar &dz)
{
	return ((field(i,j+2,k,m)) - 4.*(field(i,j+1,k,m)) + 6.*(field(i,j,k,m)) - 4.*(field(i,j-1,k,m)) + (field(i,j-2,k,m))) / (dy*dy*dy*dy);
}


}
}
#endif
