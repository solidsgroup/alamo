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

//
// FIRST order derivatives
//

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
            return (f(i,j,k,m) - f(i-1,j,k,m)) / dx[0]; // 1st order stencil
            //return (f(i,j,k,m)*(1.5) + f(i-1,j,k,m)*(-2.0) + f(i-2,j,k,m)*(0.5)) / dx[0]; // 2nd order stencil
        else if (stencil[0] == StencilType::Hi)
            return (f(i+1,j,k,m) - f(i,j,k,m)) / dx[0]; // 1st order stencil
            //return (f(i,j,k,m)*(-1.5) + f(i+1,j,k,m)*(2.0) + f(i+2,j,k,m)*(-0.5)) / dx[0]; // 2nd order stencil
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
            return (f(i,j,k,m) - f(i,j-1,k,m)) / dx[1];
            //return (f(i,j,k,m)*( 1.5) + f(i,j-1,k,m)*(-2.0) + f(i,j-2,k,m)*(0.5)) / dx[1]; // 2nd order stencil
        else if (stencil[1] == StencilType::Hi)
            return (f(i,j+1,k,m) - f(i,j,k,m)) / dx[1];
            //return (f(i,j,k,m)*(-1.5) + f(i,j+1,k,m)*(2.0) + f(i,j+2,k,m)*(-0.5)) / dx[1]; // 2nd order stencil
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
            return (f(i,j,k,m) - f(i,j,k-1,m)) / dx[2];
            //return (f(i,j,k,m)*( 1.5) + f(i,j,k-1,m)*(-2.0) + f(i,j,k-2,m)*(0.5)) / dx[2]; // 2nd order stencil
        else if (stencil[2] == StencilType::Hi)
            return (f(i,j,k+1,m) - f(i,j,k,m)) / dx[2];
            //return (f(i,j,k,m)*(-1.5) + f(i,j,k+1,m)*( 2.0) + f(i,j,k+2,m)*(-0.5)) / dx[2]; // 2nd order stencil
        else
            return (f(i,j,k+1,m) - f(i,j,k-1,m))*0.5 / dx[2];
    };
};

//
// SECOND order derivatives
//

template<class T>
struct Stencil<T,2,0,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        if (stencil[0] == StencilType::Central)
            return (f(i+1,j,k,m) - 2.0*f(i,j,k,m) + f(i-1,j,k,m))/dx[0]/dx[0];
        else
            return 0.0 * f(i,j,k,m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T,0,2,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        if (stencil[1] == StencilType::Central)
            return (f(i,j+1,k,m) - 2.0*f(i,j,k,m) + f(i,j-1,k,m))/dx[1]/dx[1];
        else
            return 0.0 * f(i,j,k,m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T,0,0,2>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        if (stencil[2] == StencilType::Central)
            return (f(i,j,k+1,m) - 2.0*f(i,j,k,m) + f(i,j,k-1,m))/dx[2]/dx[2];
        else
            return 0.0 * f(i,j,k,m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T,1,1,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        int ihi = 1, ilo = 1, jhi = 1, jlo = 1;
        Set::Scalar ifac = 0.5, jfac = 0.5;
        if (stencil[0] == StencilType::Hi) {ilo = 0; ifac = 1.0;}
        if (stencil[0] == StencilType::Lo) {ihi = 0; ifac = 1.0;}
        if (stencil[1] == StencilType::Hi) {jlo = 0; jfac = 1.0;}
        if (stencil[1] == StencilType::Lo) {jhi = 0; jfac = 1.0;}

        return ifac*jfac*(f(i+ihi,j+jhi,k,m) + f(i-ilo,j-jlo,k,m) - f(i+ihi,j-jlo,k,m) - f(i-ilo,j+jhi,k,m))/(dx[0]*dx[1]);
    };
};
template<class T>
struct Stencil<T,1,0,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        int khi = 1, klo = 1, ihi = 1, ilo = 1;
        Set::Scalar kfac = 0.5, ifac = 0.5;
        if (stencil[0] == StencilType::Hi) {ilo = 0; ifac = 1.0;}
        if (stencil[0] == StencilType::Lo) {ihi = 0; ifac = 1.0;}
        if (stencil[2] == StencilType::Hi) {klo = 0; kfac = 1.0;}
        if (stencil[2] == StencilType::Lo) {khi = 0; kfac = 1.0;}

        return kfac*ifac*(f(i+ihi,j,k+khi,m) + f(i-ilo,j,k-klo,m) - f(i+ihi,j,k-klo,m) - f(i-ilo,j,k+khi,m))/(dx[0]*dx[2]);
    };
};
template<class T>
struct Stencil<T,0,1,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
    {
        int jhi = 1, jlo = 1, khi = 1, klo = 1;
        Set::Scalar jfac = 0.5, kfac = 0.5;
        if (stencil[1] == StencilType::Hi) {jlo = 0; jfac = 1.0;}
        if (stencil[1] == StencilType::Lo) {jhi = 0; jfac = 1.0;}
        if (stencil[2] == StencilType::Hi) {klo = 0; kfac = 1.0;}
        if (stencil[2] == StencilType::Lo) {khi = 0; kfac = 1.0;}

        return jfac*kfac*(f(i,j+jhi,k+khi,m) + f(i,j-jlo,k-klo,m) - f(i,j+jhi,k-klo,m) - f(i,j-jlo,k+khi,m))/(dx[1]*dx[2]);
    };
};

//
// FOURTH order derivatives
//

template<class T>
struct Stencil<T,4,0,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((f(i+2,j,k,m)) - 4.*(f(i+1,j,k,m)) + 6.*(f(i,j,k,m)) - 4.*(f(i-1,j,k,m)) + (f(i-2,j,k,m))) /
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
struct Stencil<T,0,0,4>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((f(i,j,k+2,m)) - 4.*(f(i,j,k+1,m)) + 6.*(f(i,j,k,m)) - 4.*(f(i,j,k-1,m)) + (f(i,j,k-2,m))) /
                (dx[2]*dx[2]*dx[2]*dx[2]);
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_1\partial x_2\f]
template<class T>
struct Stencil<T,3,1,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((   -f(i+2,j+2,k,m) + 8.0*f(i+2,j+1,k,m) - 8.0*f(i+2,j-1,k,m) + f(i+2,j-2,k,m))
                -2*(-f(i+1,j+2,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i+1,j-1,k,m) + f(i+1,j-2,k,m))
                +2*(-f(i-1,j+2,k,m) + 8.0*f(i-1,j+1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-1,j-2,k,m))
                    -(-f(i-2,j+2,k,m) + 8.0*f(i-2,j+1,k,m) - 8.0*f(i-2,j-1,k,m) + f(i-2,j-2,k,m))) /
            (24.0*dx[0]*dx[0]*dx[0]*dx[1]);
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_2\partial x_1\f]
template<class T>
struct Stencil<T,1,3,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((-f(i+2,j+2,k,m) + 8.0*f(i+1,j+2,k,m) - 8.0*f(i-1,j+2,k,m) + f(i-2,j+2,k,m))
            -2*(-f(i+2,j+1,k,m) + 8.0*f(i+1,j+1,k,m) - 8.0*f(i-1,j+1,k,m) + f(i-2,j+1,k,m))
            +2*(-f(i+2,j-1,k,m) + 8.0*f(i+1,j-1,k,m) - 8.0*f(i-1,j-1,k,m) + f(i-2,j-1,k,m))
            -(  -f(i+2,j-2,k,m) + 8.0*f(i+1,j-2,k,m) - 8.0*f(i-1,j-2,k,m) + f(i-2,j-2,k,m))) /
            (24.0*dx[0]*dx[1]*dx[1]*dx[1]) ;
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_2\partial x_3\f]
template<class T>
struct Stencil<T,0,3,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((   -f(i,j+2,k+2,m) + 8.0*f(i,j+2,k+1,m) - 8.0*f(i,j+2,k-1,m) + f(i,j+2,k-2,m))
                -2*(-f(i,j+1,k+2,m) + 8.0*f(i,j+1,k+1,m) - 8.0*f(i,j+1,k-1,m) + f(i,j+1,k-2,m))
                +2*(-f(i,j-1,k+2,m) + 8.0*f(i,j-1,k+1,m) - 8.0*f(i,j-1,k-1,m) + f(i,j-1,k-2,m))
                    -(-f(i,j-2,k+2,m) + 8.0*f(i,j-2,k+1,m) - 8.0*f(i,j-2,k-1,m) + f(i,j-2,k-2,m))) /
            (24.0*dx[1]*dx[1]*dx[1]*dx[2]);
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_3\partial x_2}\f]
template<class T>
struct Stencil<T,0,1,3>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((   -f(i,j+2,k+2,m) + 8.0*f(i,j+1,k+2,m) - 8.0*f(i,j-1,k+2,m) + f(i,j-2,k+2,m))
                -2*(-f(i,j+2,k+1,m) + 8.0*f(i,j+1,k+1,m) - 8.0*f(i,j-1,k+1,m) + f(i,j-2,k+1,m))
                +2*(-f(i,j+2,k-1,m) + 8.0*f(i,j+1,k-1,m) - 8.0*f(i,j-1,k-1,m) + f(i,j-2,k-1,m))
                    -(-f(i,j+2,k-2,m) + 8.0*f(i,j+1,k-2,m) - 8.0*f(i,j-1,k-2,m) + f(i,j-2,k-2,m))) /
            (24.0*dx[1]*dx[2]*dx[2]*dx[2]);
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_3\partial x_1}\f]
template<class T>
struct Stencil<T,1,0,3>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((   -f(i+2,j,k+2,m) + 8.0*f(i+1,j,k+2,m) - 8.0*f(i-1,j,k+2,m) + f(i-2,j,k+2,m))
                -2*(-f(i+2,j,k+1,m) + 8.0*f(i+1,j,k+1,m) - 8.0*f(i-1,j,k+1,m) + f(i-2,j,k+1,m))
                +2*(-f(i+2,j,k-1,m) + 8.0*f(i+1,j,k-1,m) - 8.0*f(i-1,j,k-1,m) + f(i-2,j,k-1,m))
                    -(-f(i+2,j,k-2,m) + 8.0*f(i+1,j,k-2,m) - 8.0*f(i-1,j,k-2,m) + f(i-2,j,k-2,m))) /
            (24.0*dx[0]*dx[2]*dx[2]*dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_1\partial x_3}\f]
template<class T>
struct Stencil<T,3,0,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((   -f(i+2,j,k+2,m) + 8.0*f(i+2,j,k+1,m) - 8.0*f(i+2,j,k-1,m) + f(i+2,j,k-2,m))
                -2*(-f(i+1,j,k+2,m) + 8.0*f(i+1,j,k+1,m) - 8.0*f(i+1,j,k-1,m) + f(i+1,j,k-2,m))
                +2*(-f(i-1,j,k+2,m) + 8.0*f(i-1,j,k+1,m) - 8.0*f(i-1,j,k-1,m) + f(i-1,j,k-2,m))
                    -(-f(i-2,j,k+2,m) + 8.0*f(i-2,j,k+1,m) - 8.0*f(i-2,j,k-1,m) + f(i-2,j,k-2,m))) /
            (24.0*dx[0]*dx[0]*dx[0]*dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial^2 x_2}\f]
template<class T>
struct Stencil<T,2,2,0>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(   -f(i+2,j+2,k,m)+16.0*f(i+1,j+2,k,m)-30.0*f(i  ,j+2,k,m) +16.0*f(i-1,j+2,k,m)-f(i-2,j+2,k,m))
                +16*(-f(i+2,j+1,k,m)+16.0*f(i+1,j+1,k,m)-30.0*f(i  ,j+1,k,m) +16.0*f(i-1,j+1,k,m)-f(i-2,j+1,k,m))
                -30*(-f(i+2,j  ,k,m)+16.0*f(i+1,j  ,k,m)-30.0*f(i  ,j  ,k,m) +16.0*f(i-1,j  ,k,m)-f(i-2,j  ,k,m))
                +16*(-f(i+2,j-1,k,m)+16.0*f(i+1,j-1,k,m)-30.0*f(i  ,j-1,k,m) +16.0*f(i-1,j-1,k,m)-f(i-2,j-1,k,m))
                -(   -f(i+2,j-2,k,m)+16.0*f(i+1,j-2,k,m)-30.0*f(i  ,j-2,k,m) +16.0*f(i-1,j-2,k,m)-f(i-2,j-2,k,m))) /
            (144.0*dx[0]*dx[0]*dx[1]*dx[1]) ;
    };
};

/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_2\partial^2 x_3}\f]
template<class T>
struct Stencil<T,0,2,2>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(   -f(i,j+2,k+2,m)+16.0*f(i,j+2,k+1,m)-30.0*f(i,j+2,k,m) +16.0*f(i,j+2,k-1,m)-f(i,j+2,k-2,m))
                +16*(-f(i,j+1,k+2,m)+16.0*f(i,j+1,k+1,m)-30.0*f(i,j+1,k,m) +16.0*f(i,j+1,k-1,m)-f(i,j+1,k-2,m))
                -30*(-f(i,j  ,k+2,m)+16.0*f(i,j  ,k+1,m)-30.0*f(i,j  ,k,m) +16.0*f(i,j  ,k-1,m)-f(i,j  ,k-2,m))
                +16*(-f(i,j-1,k+2,m)+16.0*f(i,j-1,k+1,m)-30.0*f(i,j-1,k,m) +16.0*f(i,j-1,k-1,m)-f(i,j-1,k-2,m))
                -(   -f(i,j-2,k+2,m)+16.0*f(i,j-2,k+1,m)-30.0*f(i,j-2,k,m) +16.0*f(i,j-2,k-1,m)-f(i,j-2,k-2,m))) /
            (144.0*dx[0]*dx[0]*dx[2]*dx[2]) ;
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial^2 x_3}\f]
template<class T>
struct Stencil<T,2,0,2>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
                const int &i, const int &j, const int &k, const int &m,
                const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(   -f(i+2,j,k+2,m)+16.0*f(i+2,j,k+1,m)-30.0*f(i+2,j,k,m) +16.0*f(i+2,j,k-1,m)-f(i+2,j,k-2,m))
                +16*(-f(i+1,j,k+2,m)+16.0*f(i+1,j,k+1,m)-30.0*f(i+1,j,k,m) +16.0*f(i+1,j,k-1,m)-f(i+1,j,k-2,m))
                -30*(-f(i  ,j,k+2,m)+16.0*f(i  ,j,k+1,m)-30.0*f(i  ,j,k,m) +16.0*f(i  ,j,k-1,m)-f(i  ,j,k-2,m))
                +16*(-f(i-1,j,k+2,m)+16.0*f(i-1,j,k+1,m)-30.0*f(i-1,j,k,m) +16.0*f(i-1,j,k-1,m)-f(i-1,j,k-2,m))
                -(   -f(i-2,j,k+2,m)+16.0*f(i-2,j,k+1,m)-30.0*f(i-2,j,k,m) +16.0*f(i-2,j,k-1,m)-f(i-2,j,k-2,m))) /
            (144.0*dx[0]*dx[0]*dx[2]*dx[2]) ;
    };
};


/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial x_2\partial x_3}\f]
template<class T>
struct Stencil<T,2,1,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ( + (f(i+1,j+1,k+1,m) - 2.0*f(i,j+1,k+1,m) + f(i-1,j+1,k+1,m))
                - (f(i+1,j-1,k+1,m) - 2.0*f(i,j-1,k+1,m) + f(i-1,j-1,k+1,m))
                - (f(i+1,j+1,k-1,m) - 2.0*f(i,j+1,k-1,m) + f(i-1,j+1,k-1,m))
                + (f(i+1,j-1,k-1,m) - 2.0*f(i,j-1,k-1,m) + f(i-1,j-1,k-1,m)))
                / (4.0*dx[0]*dx[0]*dx[1]*dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial x_1\partial^2 x_2\partial x_3}\f]
template<class T>
struct Stencil<T,1,2,1>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ( + (f(i+1,j+1,k+1,m) - 2.0*f(i+1,j,k+1,m) + f(i+1,j-1,k+1,m))
                - (f(i-1,j+1,k+1,m) - 2.0*f(i-1,j,k+1,m) + f(i-1,j-1,k+1,m))
                - (f(i+1,j+1,k-1,m) - 2.0*f(i+1,j,k-1,m) + f(i+1,j-1,k-1,m))
                + (f(i-1,j+1,k-1,m) - 2.0*f(i-1,j,k-1,m) + f(i-1,j-1,k-1,m)))
                / (4.0*dx[0]*dx[1]*dx[1]*dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial x_1\partial x_2\partial^2 x_3}\f]
template<class T>
struct Stencil<T,1,1,2>
{
    AMREX_FORCE_INLINE
    static T D(const amrex::Array4<const T> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ( + (f(i+1,j+1,k+1,m) - 2.0*f(i+1,j+1,k,m) + f(i+1,j+1,k-1,m))
                - (f(i-1,j+1,k+1,m) - 2.0*f(i-1,j+1,k,m) + f(i-1,j+1,k-1,m))
                - (f(i+1,j-1,k+1,m) - 2.0*f(i+1,j-1,k,m) + f(i+1,j-1,k-1,m))
                + (f(i-1,j-1,k+1,m) - 2.0*f(i-1,j-1,k,m) + f(i-1,j-1,k-1,m)))
                / (4.0*dx[0]*dx[1]*dx[1]*dx[2]);

    };
};

AMREX_FORCE_INLINE
Set::Scalar
Laplacian(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k, const int &m,
        const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Scalar ret = 0.0;
    ret += (Numeric::Stencil<Set::Scalar,2,0,0>::D(f,i,j,k,m,dx));
#if AMREX_SPACEDIM > 1
    ret += (Numeric::Stencil<Set::Scalar,0,2,0>::D(f,i,j,k,m,dx));
#if AMREX_SPACEDIM > 2
    ret += (Numeric::Stencil<Set::Scalar,0,0,2>::D(f,i,j,k,m,dx));
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
Set::Vector
Divergence(const amrex::Array4<const Set::Matrix> &dw,
        const int &i, const int &j, const int &k,
        const Set::Scalar DX[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> &stencil = DefaultType)
{
    Set::Vector ret = Set::Vector::Zero();
    
    if (stencil[0] == StencilType::Central)
    {
        AMREX_D_TERM(ret(0) += (dw(i + 1, j, k)(0, 0) - dw(i - 1, j, k)(0, 0)) / 2. / DX[0];,
                    ret(1) += (dw(i + 1, j, k)(1, 0) - dw(i - 1, j, k)(1, 0)) / 2. / DX[0];, 
                    ret(2) += (dw(i + 1, j, k)(2, 0) - dw(i - 1, j, k)(2, 0)) / 2. / DX[0];)
    }
    else if (stencil[0] == StencilType::Lo)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j, k)(0, 0) - dw(i - 1, j, k)(0, 0)) / DX[0];,
                    ret(1) += (dw(i, j, k)(1, 0) - dw(i - 1, j, k)(1, 0)) / DX[0];,
                    ret(2) += (dw(i, j, k)(2, 0) - dw(i - 1, j, k)(2, 0)) / DX[0];)
    }
    else if (stencil[0] == StencilType::Hi)
    {
        AMREX_D_TERM(ret(0) += (dw(i + 1, j, k)(0, 0) - dw(i, j, k)(0, 0)) / DX[0];,
                    ret(1) += (dw(i + 1, j, k)(1, 0) - dw(i, j, k)(1, 0)) / DX[0];,
                    ret(2) += (dw(i + 1, j, k)(2, 0) - dw(i, j, k)(2, 0)) / DX[0];)
    }

#if AMREX_SPACEDIM > 1
    if (stencil[1] == StencilType::Central)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j + 1, k)(0, 1) - dw(i, j - 1, k)(0, 1)) / 2. / DX[1];,
                    ret(1) += (dw(i, j + 1, k)(1, 1) - dw(i, j - 1, k)(1, 1)) / 2. / DX[1];,
                    ret(2) += (dw(i, j + 1, k)(2, 1) - dw(i, j - 1, k)(2, 1)) / 2. / DX[1];)
    }
    else if (stencil[1] == StencilType::Lo)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j, k)(0, 1) - dw(i, j - 1, k)(0, 1)) / DX[1];,
                    ret(1) += (dw(i, j, k)(1, 1) - dw(i, j - 1, k)(1, 1)) / DX[1];,
                    ret(2) += (dw(i, j, k)(2, 1) - dw(i, j - 1, k)(2, 1)) / DX[1];)
    }
    else if (stencil[1] == StencilType::Hi)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j + 1, k)(0, 1) - dw(i, j, k)(0, 1)) / DX[1];,
                    ret(1) += (dw(i, j + 1, k)(1, 1) - dw(i, j, k)(1, 1)) / DX[1];,
                    ret(2) += (dw(i, j + 1, k)(2, 1) - dw(i, j, k)(2, 1)) / DX[1];)
    }
#endif         
#if AMREX_SPACEDIM > 2
    if (stencil[2] == StencilType::Central)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j, k + 1)(0, 2) - dw(i, j, k - 1)(0, 2)) / 2. / DX[2];,
                    ret(1) += (dw(i, j, k + 1)(1, 2) - dw(i, j, k - 1)(1, 2)) / 2. / DX[2];,
                    ret(2) += (dw(i, j, k + 1)(2, 2) - dw(i, j, k - 1)(2, 2)) / 2. / DX[2];)
    }
    else if (stencil[2] == StencilType::Lo)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j, k)(0, 2) - dw(i, j, k - 1)(0, 2)) / DX[2];,
                    ret(1) += (dw(i, j, k)(1, 2) - dw(i, j, k - 1)(1, 2)) / DX[2];,
                    ret(2) += (dw(i, j, k)(2, 2) - dw(i, j, k - 1)(2, 2)) / DX[2];)
    }
    else if (stencil[2] == StencilType::Hi)
    {
        AMREX_D_TERM(ret(0) += (dw(i, j, k + 1)(0, 2) - dw(i, j, k)(0, 2)) / DX[2];,
                    ret(1) += (dw(i, j, k + 1)(1, 2) - dw(i, j, k)(1, 2)) / DX[2];,
                    ret(2) += (dw(i, j, k + 1)(2, 2) - dw(i, j, k)(2, 2)) / DX[2];)
    }
#endif
    return ret;
}


AMREX_FORCE_INLINE
Set::Vector
Gradient(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k, const int &m,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Vector ret;
    ret(0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,m,dx,stencil));
#if AMREX_SPACEDIM > 1
    ret(1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,m,dx,stencil));
#if AMREX_SPACEDIM > 2
    ret(2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,m,dx,stencil));
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
Set::Matrix
Gradient(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Matrix ret;
    ret(0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,0,dx,stencil));
#if AMREX_SPACEDIM > 1
    ret(0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,0,dx,stencil));
    ret(1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,1,dx,stencil));
    ret(1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,1,dx,stencil));
#if AMREX_SPACEDIM > 2
    ret(0,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,0,dx,stencil));
    ret(2,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,2,dx,stencil));
    ret(1,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,1,dx,stencil));
    ret(2,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,2,dx,stencil));
    ret(2,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,2,dx,stencil));
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
Set::Matrix
Gradient(const amrex::Array4<const Set::Vector> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Matrix ret;

#if AMREX_SPACEDIM > 0
    ret.col(0) = Numeric::Stencil<Set::Vector,1,0,0>::D(f,i,j,k,0,dx,stencil);
#endif
#if AMREX_SPACEDIM > 1
    ret.col(1) = Numeric::Stencil<Set::Vector,0,1,0>::D(f,i,j,k,0,dx,stencil);
#endif
#if AMREX_SPACEDIM > 2
    ret.col(2) = Numeric::Stencil<Set::Vector,0,0,1>::D(f,i,j,k,0,dx,stencil);
#endif

    return ret;
}

AMREX_FORCE_INLINE
Set::Matrix3
Gradient(const amrex::Array4<const Set::Matrix> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Matrix3 ret;

#if AMREX_SPACEDIM > 0
    ret[0] = Numeric::Stencil<Set::Matrix,1,0,0>::D(f,i,j,k,0,dx,stencil);
#endif
#if AMREX_SPACEDIM > 1
    ret[1] = Numeric::Stencil<Set::Matrix,0,1,0>::D(f,i,j,k,0,dx,stencil);
#endif
#if AMREX_SPACEDIM > 2
    ret[2] = Numeric::Stencil<Set::Matrix,0,0,1>::D(f,i,j,k,0,dx,stencil);
#endif

    return ret;
}


AMREX_FORCE_INLINE
Set::Matrix3
MatrixGradient(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Matrix3 ret;
#if AMREX_SPACEDIM == 1
    ret[0](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,0,dx,stencil));
#elif AMREX_SPACEDIM == 2
    ret[0](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,0,dx,stencil));
    ret[0](0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,0,dx,stencil));
    ret[0](1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,1,dx,stencil));
    ret[0](1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,1,dx,stencil));
    ret[1](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,2,dx,stencil));
    ret[1](0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,2,dx,stencil));
    ret[1](1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,3,dx,stencil));
    ret[1](1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,3,dx,stencil));
#elif AMREX_SPACEDIM == 3
    ret[0](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,0,dx,stencil));
    ret[0](1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,1,dx,stencil));
    ret[0](2,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,2,dx,stencil));
    ret[1](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,3,dx,stencil));
    ret[1](1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,4,dx,stencil));
    ret[1](2,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,5,dx,stencil));
    ret[2](0,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,6,dx,stencil));
    ret[2](1,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,7,dx,stencil));
    ret[2](2,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(f,i,j,k,8,dx,stencil));
    ret[0](0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,0,dx,stencil));
    ret[0](1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,1,dx,stencil));
    ret[0](2,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,2,dx,stencil));
    ret[1](0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,3,dx,stencil));
    ret[1](1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,4,dx,stencil));
    ret[1](2,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,5,dx,stencil));
    ret[2](0,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,6,dx,stencil));
    ret[2](1,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,7,dx,stencil));
    ret[2](2,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(f,i,j,k,8,dx,stencil));
    ret[0](0,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,0,dx,stencil));
    ret[0](1,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,1,dx,stencil));
    ret[0](2,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,2,dx,stencil));
    ret[1](0,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,3,dx,stencil));
    ret[1](1,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,4,dx,stencil));
    ret[1](2,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,5,dx,stencil));
    ret[2](0,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,6,dx,stencil));
    ret[2](1,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,7,dx,stencil));
    ret[2](2,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(f,i,j,k,8,dx,stencil));
#endif
    return ret;
}


AMREX_FORCE_INLINE
Set::Matrix
Hessian(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k, const int &m,
        const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix ret;
    ret(0,0) = (Numeric::Stencil<Set::Scalar,2,0,0>::D(f,i,j,k,m,dx));
#if AMREX_SPACEDIM > 1
    ret(1,1) = (Numeric::Stencil<Set::Scalar,0,2,0>::D(f,i,j,k,m,dx));
    ret(0,1) = (Numeric::Stencil<Set::Scalar,1,1,0>::D(f,i,j,k,m,dx));
    ret(1,0) = ret(0,1);
#if AMREX_SPACEDIM > 2
    ret(2,2) = (Numeric::Stencil<Set::Scalar,0,0,2>::D(f,i,j,k,m,dx));
    ret(1,2) = (Numeric::Stencil<Set::Scalar,0,1,1>::D(f,i,j,k,m,dx));
    ret(2,0) = (Numeric::Stencil<Set::Scalar,1,0,1>::D(f,i,j,k,m,dx));
    ret(2,1) = ret(1,2);
    ret(0,2) = ret(2,0);
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
Set::Matrix3
Hessian(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar DX[AMREX_SPACEDIM],
        std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
{
    Set::Matrix3 ret;
    // 1D 
    #if AMREX_SPACEDIM>0
    ret(0, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 0, DX,stencil));
    #endif
    // 2D
    #if AMREX_SPACEDIM>1
    ret(0, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 0, DX,stencil));
    ret(0, 1, 0) = ret(0, 0, 1);
    ret(0, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 0, DX,stencil));
    ret(1, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 1, DX,stencil));
    ret(1, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 1, DX,stencil));
    ret(1, 1, 0) = ret(1, 0, 1);
    ret(1, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 1, DX,stencil));
    #endif 
    // 3D
    #if AMREX_SPACEDIM>2
    ret(0, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 0, DX,stencil));
    ret(0, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 0, DX,stencil));
    ret(0, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 0, DX,stencil));
    ret(0, 2, 0) = ret(0, 0, 2);
    ret(0, 2, 1) = ret(0, 1, 2);;
    ret(1, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 1, DX,stencil));
    ret(1, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 1, DX,stencil));
    ret(1, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 1, DX,stencil));
    ret(1, 2, 0) = ret(1, 0, 2);
    ret(1, 2, 1) = ret(1, 1, 2);;
    ret(2, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 1, 0) = ret(2, 0, 1);
    ret(2, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 2, DX,stencil));
    ret(2, 2, 0) = ret(2, 0, 2);
    ret(2, 2, 1) = ret(2, 1, 2);;
    #endif
    return ret;
}


// Returns Hessian of a vector field.
// Return value: ret[i](j,k) = ret_{i,jk}
AMREX_FORCE_INLINE
Set::Matrix3
Hessian(const amrex::Array4<const Set::Vector> &f,
        const int &i, const int &j, const int &k,
        const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix3 ret;

    #if AMREX_SPACEDIM>0    
    Set::Vector f_11 = Numeric::Stencil<Set::Vector,2,0,0>::D(f,i,j,k,0,dx);
    ret[0](0,0) = f_11(0);
    #endif
    #if AMREX_SPACEDIM>1
    ret[1](0,0) = f_11(1);
    Set::Vector f_12 = Numeric::Stencil<Set::Vector,1,1,0>::D(f,i,j,k,0,dx);
    ret[0](1,0) = ret[0](0,1) = f_12(0);
    ret[1](1,0) = ret[1](0,1) = f_12(1);
    Set::Vector f_22 = Numeric::Stencil<Set::Vector,0,2,0>::D(f,i,j,k,0,dx);
    ret[0](1,1) = f_22(0);
    ret[1](1,1) = f_22(1);
    #endif
    #if AMREX_SPACEDIM>2
    ret[2](0,0) = f_11(2);
    ret[2](1,0) = ret[2](0,1) = f_12(2);
    Set::Vector f_13 = Numeric::Stencil<Set::Vector,1,0,1>::D(f,i,j,k,0,dx);
    ret[0](2,0) = ret[0](0,2) = f_13(0);
    ret[1](2,0) = ret[1](0,2) = f_13(1);
    ret[2](2,0) = ret[2](0,2) = f_13(2);
    ret[2](1,1) = f_22(2);
    Set::Vector f_23 = Numeric::Stencil<Set::Vector,0,1,1>::D(f,i,j,k,0,dx);
    ret[0](1,2) = ret[0](2,1) = f_23(0);
    ret[1](1,2) = ret[1](2,1) = f_23(1);
    ret[2](1,2) = ret[2](2,1) = f_23(2);
    Set::Vector f_33 = Numeric::Stencil<Set::Vector,0,0,2>::D(f,i,j,k,0,dx);
    ret[0](2,2) = f_33(0);
    ret[1](2,2) = f_33(1);
    ret[2](2,2) = f_33(2);
    #endif
    return ret;
}


AMREX_FORCE_INLINE
Set::Matrix
FieldToMatrix(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k)
{
    Set::Matrix ret;
#if AMREX_SPACEDIM == 1
    ret(0,0) = f(i,j,k,0);

#elif AMREX_SPACEDIM == 2
    ret(0,0) = f(i,j,k,0); ret(0,1) = f(i,j,k,1);
    ret(1,0) = f(i,j,k,2); ret(1,1) = f(i,j,k,3);

#elif AMREX_SPACEDIM == 3
    ret(0,0) = f(i,j,k,0); ret(0,1) = f(i,j,k,1); ret(0,2) = f(i,j,k,2);
    ret(1,0) = f(i,j,k,3); ret(1,1) = f(i,j,k,4); ret(1,2) = f(i,j,k,5);
    ret(2,0) = f(i,j,k,6); ret(2,1) = f(i,j,k,7); ret(2,2) = f(i,j,k,8);
#endif

    return ret;
}

AMREX_FORCE_INLINE
Set::Matrix
FieldToMatrix(const amrex::Array4<Set::Scalar> &f,
        const int &i, const int &j, const int &k)
{
    Set::Matrix ret;
#if AMREX_SPACEDIM == 1
    ret(0,0) = f(i,j,k,0);

#elif AMREX_SPACEDIM == 2
    ret(0,0) = f(i,j,k,0); ret(0,1) = f(i,j,k,1);
    ret(1,0) = f(i,j,k,2); ret(1,1) = f(i,j,k,3);

#elif AMREX_SPACEDIM == 3
    ret(0,0) = f(i,j,k,0); ret(0,1) = f(i,j,k,1); ret(0,2) = f(i,j,k,2);
    ret(1,0) = f(i,j,k,3); ret(1,1) = f(i,j,k,4); ret(1,2) = f(i,j,k,5);
    ret(2,0) = f(i,j,k,6); ret(2,1) = f(i,j,k,7); ret(2,2) = f(i,j,k,8);
#endif

    return ret;
}

AMREX_FORCE_INLINE
Set::Vector
FieldToVector(const amrex::Array4<const Set::Scalar> &f,
        const int &i, const int &j, const int &k)
{
    Set::Vector ret;
    ret(0) = f(i,j,k,0);
#if AMREX_SPACEDIM > 1
    ret(1) = f(i,j,k,1);
#if AMREX_SPACEDIM > 2
    ret(2) = f(i,j,k,2);
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
Set::Vector
FieldToVector(const amrex::Array4<Set::Scalar> &f,
        const int &i, const int &j, const int &k)
{
    Set::Vector ret;
    ret(0) = f(i,j,k,0);
#if AMREX_SPACEDIM > 1
    ret(1) = f(i,j,k,1);
#if AMREX_SPACEDIM > 2
    ret(2) = f(i,j,k,2);
#endif
#endif
    return ret;
}

AMREX_FORCE_INLINE
void
MatrixToField(const amrex::Array4<Set::Scalar> &f,
        const int &i, const int &j, const int &k,
        Set::Matrix matrix)
{
#if AMREX_SPACEDIM == 1
    f(i,j,k,0) = matrix(0,0);
#elif AMREX_SPACEDIM == 2
    f(i,j,k,0) = matrix(0,0); f(i,j,k,1) = matrix(0,1);
    f(i,j,k,2) = matrix(1,0); f(i,j,k,3) = matrix(1,1);
#elif AMREX_SPACEDIM == 3
    f(i,j,k,0) = matrix(0,0); f(i,j,k,1) = matrix(0,1); f(i,j,k,2) = matrix(0,2);
    f(i,j,k,3) = matrix(1,0); f(i,j,k,4) = matrix(1,1); f(i,j,k,5) = matrix(1,2);
    f(i,j,k,6) = matrix(2,0); f(i,j,k,7) = matrix(2,1); f(i,j,k,8) = matrix(2,2);
#endif
}

AMREX_FORCE_INLINE
void
VectorToField(const amrex::Array4<Set::Scalar> &f,
        const int &i, const int &j, const int &k,
        Set::Vector vector)
{
    f(i,j,k,0) = vector(0);
#if AMREX_SPACEDIM > 1
    f(i,j,k,1) = vector(1);
#if AMREX_SPACEDIM > 2
    f(i,j,k,2) = vector(2);
#endif
#endif
}

template<int index,int SYM>
Set::Matrix3
Divergence(const amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM,SYM>> &,
            const int , const int , const int , const Set::Scalar[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType)
            {Util::Abort(INFO,"Not implemented yet"); return Set::Matrix3::Zero();}

template<>
AMREX_FORCE_INLINE
Set::Matrix3 Divergence<2,Set::Sym::Isotropic>(const amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>> &C,
            const int i, const int j, const int k, const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil)
{
    Set::Matrix3 ret;
    
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> gradCx =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>,1,0,0>::D(C,i,j,k,0,dx,stencil);
    #if AMREX_SPACEDIM>1
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> gradCy =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>,0,1,0>::D(C,i,j,k,0,dx,stencil);
    #endif
    #if AMREX_SPACEDIM>2
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> gradCz =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>,0,0,1>::D(C,i,j,k,0,dx,stencil);
    #endif

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int k = 0; k < AMREX_SPACEDIM; k++)
            for (int l = 0; l < AMREX_SPACEDIM; l++)
                {
                    ret(i,k,l) = 0.0;
                    ret(i,k,l) = gradCx(i,0,k,l);
                    #if AMREX_SPACEDIM>1
                    ret(i,k,l) = gradCy(i,1,k,l);
                    #endif
                    #if AMREX_SPACEDIM>2
                    ret(i,k,l) = gradCz(i,2,k,l);
                    #endif
                }
    return ret;
}




template<int dim>
AMREX_FORCE_INLINE
Set::Matrix4<dim,Set::Sym::Full>
DoubleHessian(const amrex::Array4<const Set::Scalar> &f,
                const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM]);

template<>
AMREX_FORCE_INLINE
Set::Matrix4<2,Set::Sym::Full>
DoubleHessian<2>(const amrex::Array4<const Set::Scalar> &f,
                const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix4<2,Set::Sym::Full> ret;
    // [0,0,0,0]
    ret(0, 0, 0, 0) = Stencil<Set::Scalar,4,0,0>::D(f,i,j,k,m,dx);
    // [0, 0, 0, 1]
    ret(0, 0, 0, 1) = Stencil<Set::Scalar,3,1,0>::D(f,i,j,k,m,dx);
    // [0, 0, 1, 1]
    ret(0, 0, 1, 1) = Stencil<Set::Scalar,2,2,0>::D(f,i,j,k,m,dx);
    // [0, 1, 1, 1]
    ret(0, 1, 1, 1) = Stencil<Set::Scalar,1,3,0>::D(f,i,j,k,m,dx);
    // [1, 1, 1, 1]
    ret(1, 1, 1, 1) = Stencil<Set::Scalar,0,4,0>::D(f,i,j,k,m,dx);
    return ret;
}

template<>
AMREX_FORCE_INLINE
Set::Matrix4<3,Set::Sym::Full>
DoubleHessian<3>(const amrex::Array4<const Set::Scalar> &f,
                const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix4<3,Set::Sym::Full> ret;
    // [0,0,0,0]
    ret(0, 0, 0, 0) = Stencil<Set::Scalar,4,0,0>::D(f,i,j,k,m,dx);
    // [0, 0, 0, 1]
    ret(0, 0, 0, 1) = Stencil<Set::Scalar,3,1,0>::D(f,i,j,k,m,dx);
    // [0, 0, 0, 2]
    ret(0, 0, 0, 2) = Stencil<Set::Scalar,3,0,1>::D(f,i,j,k,m,dx);
    // [0, 0, 1, 1]
    ret(0, 0, 1, 1) = Stencil<Set::Scalar,2,2,0>::D(f,i,j,k,m,dx);
    // [0, 0, 1, 2]
    ret(0, 0, 1, 2) = Stencil<Set::Scalar,2,1,1>::D(f,i,j,k,m,dx);
    // [0, 0, 2, 2]
    ret(0, 0, 2, 2) = Stencil<Set::Scalar,2,0,2>::D(f,i,j,k,m,dx);
    // [0, 1, 1, 1]
    ret(0, 1, 1, 1) = Stencil<Set::Scalar,1,3,0>::D(f,i,j,k,m,dx);
    // [0, 1, 1, 2]
    ret(0, 1, 1, 2) = Stencil<Set::Scalar,1,2,1>::D(f,i,j,k,m,dx);
    // [0, 1, 2, 2]
    ret(0, 1, 2, 2) = Stencil<Set::Scalar,1,1,2>::D(f,i,j,k,m,dx);
    // [0, 2, 2, 2]
    ret(0, 2, 2, 2) = Stencil<Set::Scalar,1,0,3>::D(f,i,j,k,m,dx);
    // [1, 1, 1, 1]
    ret(1, 1, 1, 1) = Stencil<Set::Scalar,0,4,0>::D(f,i,j,k,m,dx);
    // [1, 1, 1, 2]
    ret(1, 1, 1, 2) = Stencil<Set::Scalar,0,3,1>::D(f,i,j,k,m,dx);
    // [1, 1, 2, 2]
    ret(1, 1, 2, 2) = Stencil<Set::Scalar,0,2,2>::D(f,i,j,k,m,dx);
    // [1, 2, 2, 2]
    ret(1, 2, 2, 2) = Stencil<Set::Scalar,0,1,3>::D(f,i,j,k,m,dx);
    // [2, 2, 2, 2]
    ret(2, 2, 2, 2) = Stencil<Set::Scalar,0,0,4>::D(f,i,j,k,m,dx);
    return ret;
}

struct Interpolate
{
public:
    template<class T>
    AMREX_FORCE_INLINE
    static T CellToNodeAverage(const amrex::Array4<const T> &f,
                                const int &i, const int &j, const int &k, const int &m)
    {
        return ( AMREX_D_TERM(      f(i  ,j  ,k  ,m) + f(i-1,j  ,k  ,m)
                                , 
                                + f(i  ,j-1,k  ,m) + f(i-1,j-1,k  ,m)
                                , 
                                + f(i  ,j  ,k-1,m) + f(i-1,j  ,k-1,m)
                                + f(i  ,j-1,k-1,m) + f(i-1,j-1,k-1,m)
                                )) * fac;
    }
    template<class T>
    AMREX_FORCE_INLINE
    static T NodeToCellAverage(const amrex::Array4<const T> &f,
                                const int &i, const int &j, const int &k, const int &m)
    {
        return ( AMREX_D_TERM(      f(i  ,j  ,k  ,m) + f(i+1,j  ,k  ,m)
                                , 
                                + f(i  ,j+1,k  ,m) + f(i+1,j+1,k  ,m)
                                , 
                                + f(i  ,j  ,k+1,m) + f(i+1,j  ,k+1,m)
                                + f(i  ,j+1,k+1,m) + f(i+1,j+1,k+1,m)
                                )) * fac;
    }
    constexpr static Set::Scalar fac = AMREX_D_PICK(0.5,0.25,0.125);
};

}
#endif
