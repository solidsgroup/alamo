#ifndef NUMERIC_STENCIL_H_
#define NUMERIC_STENCIL_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Set/Base.H"
#include "Set/Matrix3.H"
#include "Set/Matrix4.H"

#include "Set/Matrix4_Major.H"
#include "Set/Matrix4_MajorMinor.H"
#include "Set/Matrix4_Isotropic.H"
#include "Set/Matrix4_Diagonal.H"
#include "Set/Matrix4_Full.H"

/// \brief This namespace contains some numerical tools
namespace Numeric
{

enum StencilType { Lo, Hi, Central };
AMREX_GPU_HOST_DEVICE
static std::array<StencilType, AMREX_SPACEDIM>
DefaultType()
{
    return {AMREX_D_DECL(StencilType::Central, StencilType::Central, StencilType::Central)};
};
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
XLo         = { AMREX_D_DECL(StencilType::Lo,      StencilType::Central, StencilType::Central) };
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
XHi         = { AMREX_D_DECL(StencilType::Hi,      StencilType::Central, StencilType::Central) };
#if AMREX_SPACEDIM>1
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
YLo         = { AMREX_D_DECL(StencilType::Central, StencilType::Lo,      StencilType::Central) };
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
YHi         = { AMREX_D_DECL(StencilType::Central, StencilType::Hi,      StencilType::Central) };
#endif
#if AMREX_SPACEDIM>2
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
ZLo         = { AMREX_D_DECL(StencilType::Central, StencilType::Central, StencilType::Lo     ) };
[[maybe_unused]] static std::array<StencilType, AMREX_SPACEDIM>
ZHi         = { AMREX_D_DECL(StencilType::Central, StencilType::Central, StencilType::Hi     ) };
#endif

[[nodiscard]]
static
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
std::array<StencilType, AMREX_SPACEDIM>
GetStencil(const int i, const int j, const int k, const amrex::Box domain)
{
    (void)i; (void)j; (void)k; // Suppress "unused variable" warnings.
    std::array<StencilType, AMREX_SPACEDIM> sten;
    const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
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

template<class T, int x, int y, int z>
struct Stencil
{};

//
// FIRST order derivatives
//

template<class T>
struct Stencil<T, 1, 0, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE 
        static T D( const amrex::Array4<const T>& f,
                    const int& i, const int& j, const int& k, const int& m,
                    const Set::Scalar dx[AMREX_SPACEDIM],
                    std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[0] == StencilType::Lo)
            return (f(i, j, k, m) - f(i - 1, j, k, m)) / dx[0]; // 1st order stencil
        else if (stencil[0] == StencilType::Hi)
            return (f(i + 1, j, k, m) - f(i, j, k, m)) / dx[0]; // 1st order stencil
        else
            return (f(i + 1, j, k, m) - f(i - 1, j, k, m)) * 0.5 / dx[0];
    };
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static Set::Jet<Set::Scalar,T>
    Dsplit( const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[0] == StencilType::Lo)
            return Set::Jet<Set::Scalar,T>(
                1.0 / dx[0],
                -f(i - 1, j, k, m) / dx[0]
                ); // 1st order
        else if (stencil[0] == StencilType::Hi)
            return Set::Jet<Set::Scalar,T>(
                -1.0 / dx[0],
                f(i + 1, j, k, m) / dx[0]
                ); // 1st order
        else
            return Set::Jet<Set::Scalar,T>(
                0.0,
                (f(i + 1, j, k, m) - f(i - 1, j, k, m)) * 0.5 / dx[0]
                );
    };
};

template<class T>
struct Stencil<T, 0, 1, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[1] == StencilType::Lo)
            return (f(i, j, k, m) - f(i, j - 1, k, m)) / dx[1];
        else if (stencil[1] == StencilType::Hi)
            return (f(i, j + 1, k, m) - f(i, j, k, m)) / dx[1];
        else
            return (f(i, j + 1, k, m) - f(i, j - 1, k, m)) * 0.5 / dx[1];
    };
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static Set::Jet<Set::Scalar,T>
    Dsplit( const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[1] == StencilType::Lo)
            return Set::Jet<Set::Scalar,T>(
                1.0 / dx[1],
                - f(i, j - 1, k, m) / dx[1]
                );
        else if (stencil[1] == StencilType::Hi)
            return Set::Jet<Set::Scalar,T>(
                - 1.0 / dx[1],
                f(i, j + 1, k, m) / dx[1]
                );
        else
            return Set::Jet<Set::Scalar,T>(
                0.0,
                (f(i, j + 1, k, m) - f(i, j - 1, k, m)) * 0.5 / dx[1]
                );
    };
};

template<class T>
struct Stencil<T, 0, 0, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[2] == StencilType::Lo)
            return (f(i, j, k, m) - f(i, j, k - 1, m)) / dx[2];
        else if (stencil[2] == StencilType::Hi)
            return (f(i, j, k + 1, m) - f(i, j, k, m)) / dx[2];
        else
            return (f(i, j, k + 1, m) - f(i, j, k - 1, m)) * 0.5 / dx[2];
    };
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static Set::Jet<Set::Scalar, T>
    Dsplit( const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[2] == StencilType::Lo)
            return Set::Jet<Set::Scalar,T>(
                1.0 / dx[2],
                -f(i, j, k - 1, m) / dx[2]
                );
        else if (stencil[2] == StencilType::Hi)
            return Set::Jet<Set::Scalar,T>(
                -1.0 / dx[2],
                f(i, j, k + 1, m) / dx[2]
                );
        else
            return Set::Jet<Set::Scalar,T>(
                0.0,
                (f(i, j, k + 1, m) - f(i, j, k - 1, m)) * 0.5 / dx[2]
                );
    };
};

//
// SECOND order derivatives
//

template<class T>
struct Stencil<T, 2, 0, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[0] == StencilType::Central)
            return (f(i + 1, j, k, m) - 2.0 * f(i, j, k, m) + f(i - 1, j, k, m)) / dx[0] / dx[0];
        else
            return 0.0 * f(i, j, k, m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T, 0, 2, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[1] == StencilType::Central)
            return (f(i, j + 1, k, m) - 2.0 * f(i, j, k, m) + f(i, j - 1, k, m)) / dx[1] / dx[1];
        else
            return 0.0 * f(i, j, k, m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T, 0, 0, 2>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        if (stencil[2] == StencilType::Central)
            return (f(i, j, k + 1, m) - 2.0 * f(i, j, k, m) + f(i, j, k - 1, m)) / dx[2] / dx[2];
        else
            return 0.0 * f(i, j, k, m); // TODO this is not a great way to do this
    };
};

template<class T>
struct Stencil<T, 1, 1, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        int ihi = 1, ilo = 1, jhi = 1, jlo = 1;
        Set::Scalar ifac = 0.5, jfac = 0.5;
        if (stencil[0] == StencilType::Hi) { ilo = 0; ifac = 1.0; }
        if (stencil[0] == StencilType::Lo) { ihi = 0; ifac = 1.0; }
        if (stencil[1] == StencilType::Hi) { jlo = 0; jfac = 1.0; }
        if (stencil[1] == StencilType::Lo) { jhi = 0; jfac = 1.0; }

        return ifac * jfac * (f(i + ihi, j + jhi, k, m) + f(i - ilo, j - jlo, k, m) - f(i + ihi, j - jlo, k, m) - f(i - ilo, j + jhi, k, m)) / (dx[0] * dx[1]);
    };
};
template<class T>
struct Stencil<T, 1, 0, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        int khi = 1, klo = 1, ihi = 1, ilo = 1;
        Set::Scalar kfac = 0.5, ifac = 0.5;
        if (stencil[0] == StencilType::Hi) { ilo = 0; ifac = 1.0; }
        if (stencil[0] == StencilType::Lo) { ihi = 0; ifac = 1.0; }
        if (stencil[2] == StencilType::Hi) { klo = 0; kfac = 1.0; }
        if (stencil[2] == StencilType::Lo) { khi = 0; kfac = 1.0; }

        return kfac * ifac * (f(i + ihi, j, k + khi, m) + f(i - ilo, j, k - klo, m) - f(i + ihi, j, k - klo, m) - f(i - ilo, j, k + khi, m)) / (dx[0] * dx[2]);
    };
};
template<class T>
struct Stencil<T, 0, 1, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T D( const amrex::Array4<const T>& f,
                const int& i, const int& j, const int& k, const int& m,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        int jhi = 1, jlo = 1, khi = 1, klo = 1;
        Set::Scalar jfac = 0.5, kfac = 0.5;
        if (stencil[1] == StencilType::Hi) { jlo = 0; jfac = 1.0; }
        if (stencil[1] == StencilType::Lo) { jhi = 0; jfac = 1.0; }
        if (stencil[2] == StencilType::Hi) { klo = 0; kfac = 1.0; }
        if (stencil[2] == StencilType::Lo) { khi = 0; kfac = 1.0; }

        return jfac * kfac * (f(i, j + jhi, k + khi, m) + f(i, j - jlo, k - klo, m) - f(i, j + jhi, k - klo, m) - f(i, j - jlo, k + khi, m)) / (dx[1] * dx[2]);
    };
};

//
// FOURTH order derivatives
//

template<class T>
struct Stencil<T, 4, 0, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((f(i + 2, j, k, m)) - 4. * (f(i + 1, j, k, m)) + 6. * (f(i, j, k, m)) - 4. * (f(i - 1, j, k, m)) + (f(i - 2, j, k, m))) /
            (dx[0] * dx[0] * dx[0] * dx[0]);
    };
};
template<class T>
struct Stencil<T, 0, 4, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((f(i, j + 2, k, m)) - 4. * (f(i, j + 1, k, m)) + 6. * (f(i, j, k, m)) - 4. * (f(i, j - 1, k, m)) + (f(i, j - 2, k, m))) /
            (dx[1] * dx[1] * dx[1] * dx[1]);
    };
};
template<class T>
struct Stencil<T, 0, 0, 4>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((f(i, j, k + 2, m)) - 4. * (f(i, j, k + 1, m)) + 6. * (f(i, j, k, m)) - 4. * (f(i, j, k - 1, m)) + (f(i, j, k - 2, m))) /
            (dx[2] * dx[2] * dx[2] * dx[2]);
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_1\partial x_2\f]
template<class T>
struct Stencil<T, 3, 1, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((-f(i + 2, j + 2, k, m) + 8.0 * f(i + 2, j + 1, k, m) - 8.0 * f(i + 2, j - 1, k, m) + f(i + 2, j - 2, k, m))
            - 2 * (-f(i + 1, j + 2, k, m) + 8.0 * f(i + 1, j + 1, k, m) - 8.0 * f(i + 1, j - 1, k, m) + f(i + 1, j - 2, k, m))
            + 2 * (-f(i - 1, j + 2, k, m) + 8.0 * f(i - 1, j + 1, k, m) - 8.0 * f(i - 1, j - 1, k, m) + f(i - 1, j - 2, k, m))
            - (-f(i - 2, j + 2, k, m) + 8.0 * f(i - 2, j + 1, k, m) - 8.0 * f(i - 2, j - 1, k, m) + f(i - 2, j - 2, k, m))) /
            (24.0 * dx[0] * dx[0] * dx[0] * dx[1]);
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_2\partial x_1\f]
template<class T>
struct Stencil<T, 1, 3, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return ((-f(i + 2, j + 2, k, m) + 8.0 * f(i + 1, j + 2, k, m) - 8.0 * f(i - 1, j + 2, k, m) + f(i - 2, j + 2, k, m))
            - 2 * (-f(i + 2, j + 1, k, m) + 8.0 * f(i + 1, j + 1, k, m) - 8.0 * f(i - 1, j + 1, k, m) + f(i - 2, j + 1, k, m))
            + 2 * (-f(i + 2, j - 1, k, m) + 8.0 * f(i + 1, j - 1, k, m) - 8.0 * f(i - 1, j - 1, k, m) + f(i - 2, j - 1, k, m))
            - (-f(i + 2, j - 2, k, m) + 8.0 * f(i + 1, j - 2, k, m) - 8.0 * f(i - 1, j - 2, k, m) + f(i - 2, j - 2, k, m))) /
            (24.0 * dx[0] * dx[1] * dx[1] * dx[1]);
    };
};

/// Compute
/// \f[ \frac{\partial^4}{\partial^3 x_2\partial x_3\f]
template<class T>
struct Stencil<T, 0, 3, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((-f(i, j + 2, k + 2, m) + 8.0 * f(i, j + 2, k + 1, m) - 8.0 * f(i, j + 2, k - 1, m) + f(i, j + 2, k - 2, m))
            - 2 * (-f(i, j + 1, k + 2, m) + 8.0 * f(i, j + 1, k + 1, m) - 8.0 * f(i, j + 1, k - 1, m) + f(i, j + 1, k - 2, m))
            + 2 * (-f(i, j - 1, k + 2, m) + 8.0 * f(i, j - 1, k + 1, m) - 8.0 * f(i, j - 1, k - 1, m) + f(i, j - 1, k - 2, m))
            - (-f(i, j - 2, k + 2, m) + 8.0 * f(i, j - 2, k + 1, m) - 8.0 * f(i, j - 2, k - 1, m) + f(i, j - 2, k - 2, m))) /
            (24.0 * dx[1] * dx[1] * dx[1] * dx[2]);
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_3\partial x_2}\f]
template<class T>
struct Stencil<T, 0, 1, 3>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((-f(i, j + 2, k + 2, m) + 8.0 * f(i, j + 1, k + 2, m) - 8.0 * f(i, j - 1, k + 2, m) + f(i, j - 2, k + 2, m))
            - 2 * (-f(i, j + 2, k + 1, m) + 8.0 * f(i, j + 1, k + 1, m) - 8.0 * f(i, j - 1, k + 1, m) + f(i, j - 2, k + 1, m))
            + 2 * (-f(i, j + 2, k - 1, m) + 8.0 * f(i, j + 1, k - 1, m) - 8.0 * f(i, j - 1, k - 1, m) + f(i, j - 2, k - 1, m))
            - (-f(i, j + 2, k - 2, m) + 8.0 * f(i, j + 1, k - 2, m) - 8.0 * f(i, j - 1, k - 2, m) + f(i, j - 2, k - 2, m))) /
            (24.0 * dx[1] * dx[2] * dx[2] * dx[2]);
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_3\partial x_1}\f]
template<class T>
struct Stencil<T, 1, 0, 3>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((-f(i + 2, j, k + 2, m) + 8.0 * f(i + 1, j, k + 2, m) - 8.0 * f(i - 1, j, k + 2, m) + f(i - 2, j, k + 2, m))
            - 2 * (-f(i + 2, j, k + 1, m) + 8.0 * f(i + 1, j, k + 1, m) - 8.0 * f(i - 1, j, k + 1, m) + f(i - 2, j, k + 1, m))
            + 2 * (-f(i + 2, j, k - 1, m) + 8.0 * f(i + 1, j, k - 1, m) - 8.0 * f(i - 1, j, k - 1, m) + f(i - 2, j, k - 1, m))
            - (-f(i + 2, j, k - 2, m) + 8.0 * f(i + 1, j, k - 2, m) - 8.0 * f(i - 1, j, k - 2, m) + f(i - 2, j, k - 2, m))) /
            (24.0 * dx[0] * dx[2] * dx[2] * dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^3 x_1\partial x_3}\f]
template<class T>
struct Stencil<T, 3, 0, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    ((-f(i + 2, j, k + 2, m) + 8.0 * f(i + 2, j, k + 1, m) - 8.0 * f(i + 2, j, k - 1, m) + f(i + 2, j, k - 2, m))
            - 2 * (-f(i + 1, j, k + 2, m) + 8.0 * f(i + 1, j, k + 1, m) - 8.0 * f(i + 1, j, k - 1, m) + f(i + 1, j, k - 2, m))
            + 2 * (-f(i - 1, j, k + 2, m) + 8.0 * f(i - 1, j, k + 1, m) - 8.0 * f(i - 1, j, k - 1, m) + f(i - 1, j, k - 2, m))
            - (-f(i - 2, j, k + 2, m) + 8.0 * f(i - 2, j, k + 1, m) - 8.0 * f(i - 2, j, k - 1, m) + f(i - 2, j, k - 2, m))) /
            (24.0 * dx[0] * dx[0] * dx[0] * dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial^2 x_2}\f]
template<class T>
struct Stencil<T, 2, 2, 0>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(-f(i + 2, j + 2, k, m) + 16.0 * f(i + 1, j + 2, k, m) - 30.0 * f(i, j + 2, k, m) + 16.0 * f(i - 1, j + 2, k, m) - f(i - 2, j + 2, k, m))
            + 16 * (-f(i + 2, j + 1, k, m) + 16.0 * f(i + 1, j + 1, k, m) - 30.0 * f(i, j + 1, k, m) + 16.0 * f(i - 1, j + 1, k, m) - f(i - 2, j + 1, k, m))
            - 30 * (-f(i + 2, j, k, m) + 16.0 * f(i + 1, j, k, m) - 30.0 * f(i, j, k, m) + 16.0 * f(i - 1, j, k, m) - f(i - 2, j, k, m))
            + 16 * (-f(i + 2, j - 1, k, m) + 16.0 * f(i + 1, j - 1, k, m) - 30.0 * f(i, j - 1, k, m) + 16.0 * f(i - 1, j - 1, k, m) - f(i - 2, j - 1, k, m))
            - (-f(i + 2, j - 2, k, m) + 16.0 * f(i + 1, j - 2, k, m) - 30.0 * f(i, j - 2, k, m) + 16.0 * f(i - 1, j - 2, k, m) - f(i - 2, j - 2, k, m))) /
            (144.0 * dx[0] * dx[0] * dx[1] * dx[1]);
    };
};

/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_2\partial^2 x_3}\f]
template<class T>
struct Stencil<T, 0, 2, 2>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(-f(i, j + 2, k + 2, m) + 16.0 * f(i, j + 2, k + 1, m) - 30.0 * f(i, j + 2, k, m) + 16.0 * f(i, j + 2, k - 1, m) - f(i, j + 2, k - 2, m))
            + 16 * (-f(i, j + 1, k + 2, m) + 16.0 * f(i, j + 1, k + 1, m) - 30.0 * f(i, j + 1, k, m) + 16.0 * f(i, j + 1, k - 1, m) - f(i, j + 1, k - 2, m))
            - 30 * (-f(i, j, k + 2, m) + 16.0 * f(i, j, k + 1, m) - 30.0 * f(i, j, k, m) + 16.0 * f(i, j, k - 1, m) - f(i, j, k - 2, m))
            + 16 * (-f(i, j - 1, k + 2, m) + 16.0 * f(i, j - 1, k + 1, m) - 30.0 * f(i, j - 1, k, m) + 16.0 * f(i, j - 1, k - 1, m) - f(i, j - 1, k - 2, m))
            - (-f(i, j - 2, k + 2, m) + 16.0 * f(i, j - 2, k + 1, m) - 30.0 * f(i, j - 2, k, m) + 16.0 * f(i, j - 2, k - 1, m) - f(i, j - 2, k - 2, m))) /
            (144.0 * dx[0] * dx[0] * dx[2] * dx[2]);
    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial^2 x_3}\f]
template<class T>
struct Stencil<T, 2, 0, 2>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return     (-(-f(i + 2, j, k + 2, m) + 16.0 * f(i + 2, j, k + 1, m) - 30.0 * f(i + 2, j, k, m) + 16.0 * f(i + 2, j, k - 1, m) - f(i + 2, j, k - 2, m))
            + 16 * (-f(i + 1, j, k + 2, m) + 16.0 * f(i + 1, j, k + 1, m) - 30.0 * f(i + 1, j, k, m) + 16.0 * f(i + 1, j, k - 1, m) - f(i + 1, j, k - 2, m))
            - 30 * (-f(i, j, k + 2, m) + 16.0 * f(i, j, k + 1, m) - 30.0 * f(i, j, k, m) + 16.0 * f(i, j, k - 1, m) - f(i, j, k - 2, m))
            + 16 * (-f(i - 1, j, k + 2, m) + 16.0 * f(i - 1, j, k + 1, m) - 30.0 * f(i - 1, j, k, m) + 16.0 * f(i - 1, j, k - 1, m) - f(i - 1, j, k - 2, m))
            - (-f(i - 2, j, k + 2, m) + 16.0 * f(i - 2, j, k + 1, m) - 30.0 * f(i - 2, j, k, m) + 16.0 * f(i - 2, j, k - 1, m) - f(i - 2, j, k - 2, m))) /
            (144.0 * dx[0] * dx[0] * dx[2] * dx[2]);
    };
};


/// \brief Compute \f[ \frac{\partial^4}{\partial^2 x_1\partial x_2\partial x_3}\f]
template<class T>
struct Stencil<T, 2, 1, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    (+(f(i + 1, j + 1, k + 1, m) - 2.0 * f(i, j + 1, k + 1, m) + f(i - 1, j + 1, k + 1, m))
            - (f(i + 1, j - 1, k + 1, m) - 2.0 * f(i, j - 1, k + 1, m) + f(i - 1, j - 1, k + 1, m))
            - (f(i + 1, j + 1, k - 1, m) - 2.0 * f(i, j + 1, k - 1, m) + f(i - 1, j + 1, k - 1, m))
            + (f(i + 1, j - 1, k - 1, m) - 2.0 * f(i, j - 1, k - 1, m) + f(i - 1, j - 1, k - 1, m)))
            / (4.0 * dx[0] * dx[0] * dx[1] * dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial x_1\partial^2 x_2\partial x_3}\f]
template<class T>
struct Stencil<T, 1, 2, 1>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    (+(f(i + 1, j + 1, k + 1, m) - 2.0 * f(i + 1, j, k + 1, m) + f(i + 1, j - 1, k + 1, m))
            - (f(i - 1, j + 1, k + 1, m) - 2.0 * f(i - 1, j, k + 1, m) + f(i - 1, j - 1, k + 1, m))
            - (f(i + 1, j + 1, k - 1, m) - 2.0 * f(i + 1, j, k - 1, m) + f(i + 1, j - 1, k - 1, m))
            + (f(i - 1, j + 1, k - 1, m) - 2.0 * f(i - 1, j, k - 1, m) + f(i - 1, j - 1, k - 1, m)))
            / (4.0 * dx[0] * dx[1] * dx[1] * dx[2]);

    };
};
/// \brief Compute \f[ \frac{\partial^4}{\partial x_1\partial x_2\partial^2 x_3}\f]
template<class T>
struct Stencil<T, 1, 1, 2>
{
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        static T D(const amrex::Array4<const T>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM])
    {
        return    (+(f(i + 1, j + 1, k + 1, m) - 2.0 * f(i + 1, j + 1, k, m) + f(i + 1, j + 1, k - 1, m))
            - (f(i - 1, j + 1, k + 1, m) - 2.0 * f(i - 1, j + 1, k, m) + f(i - 1, j + 1, k - 1, m))
            - (f(i + 1, j - 1, k + 1, m) - 2.0 * f(i + 1, j - 1, k, m) + f(i + 1, j - 1, k - 1, m))
            + (f(i - 1, j - 1, k + 1, m) - 2.0 * f(i - 1, j - 1, k, m) + f(i - 1, j - 1, k - 1, m)))
            / (4.0 * dx[0] * dx[1] * dx[1] * dx[2]);

    };
};

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Scalar
Laplacian(  const amrex::Array4<const Set::Scalar>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Scalar ret = 0.0;
    ret += (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 1
    ret += (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 2
    ret += (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, m, dx, stencil));
#endif
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
Laplacian(  const amrex::Array4<const Set::Vector>& f,
            const int& i, const int& j, const int& k,
            const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret = Set::Vector::Zero();
    ret += (Numeric::Stencil<Set::Vector, 2, 0, 0>::D(f, i, j, k, 0, dx));
#if AMREX_SPACEDIM > 1
    ret += (Numeric::Stencil<Set::Vector, 0, 2, 0>::D(f, i, j, k, 0, dx));
#if AMREX_SPACEDIM > 2
    ret += (Numeric::Stencil<Set::Vector, 0, 0, 2>::D(f, i, j, k, 0, dx));
#endif
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
Divergence( const amrex::Array4<const Set::Matrix>& dw,
            const int& i, const int& j, const int& k,
            const Set::Scalar DX[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
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

/// \brief Backward stress divergence paired with `ForwardGradient`.
///
/// This is intentionally a different operator from `Divergence`: composing
/// the two gives the ordinary one-cell second-difference rather than the
/// checkerboard-null central/central composition.
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
BackwardDivergence(const amrex::Array4<const Set::Matrix>& dw,
                   const int& i, const int& j, const int& k,
                   const Set::Scalar DX[AMREX_SPACEDIM])
{
    Set::Vector ret = Set::Vector::Zero();

    AMREX_D_TERM(
        ret += (dw(i, j, k).col(0) - dw(i - 1, j, k).col(0)) / DX[0];,
        ret += (dw(i, j, k).col(1) - dw(i, j - 1, k).col(1)) / DX[1];,
        ret += (dw(i, j, k).col(2) - dw(i, j, k - 1).col(2)) / DX[2];);

    return ret;
}


[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Scalar
Divergence( const amrex::Array4<const Set::Scalar> &f,
            const int &i, const int &j, const int &k, const int &m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType,AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Scalar ret;
    ret = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 1
    ret += (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, m, dx, stencil));
#endif
#if AMREX_SPACEDIM > 2
    ret += (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, m, dx, stencil));
#endif
    return ret;
}


[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
Gradient(   const amrex::Array4<const Set::Scalar>& f,
            const int& i, const int& j, const int& k, const int& m,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Vector ret;
    ret(0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 1
    ret(1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 2
    ret(2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, m, dx, stencil));
#endif
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
CellGradientOnNode(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret;
#if AMREX_SPACEDIM == 1
    ret(0) = (f(i, j, k, m) - f(i - 1, j, k, m)) / dx[0];
#elif AMREX_SPACEDIM == 2
    ret(0) = 0.5 * (f(i, j, k, m) - f(i - 1, j, k, m) + f(i, j - 1, k, m) - f(i - 1, j - 1, k, m)) / dx[0];
    ret(1) = 0.5 * (f(i, j, k, m) - f(i, j - 1, k, m) + f(i - 1, j, k, m) - f(i - 1, j - 1, k, m)) / dx[1];
#elif AMREX_SPACEDIM == 3
    ret(0) = 0.25 * (f(i, j, k, m) - f(i - 1, j, k, m) + f(i, j - 1, k, m) - f(i - 1, j - 1, k, m) + f(i, j, k - 1, m) - f(i - 1, j, k - 1, m) + f(i, j - 1, k - 1, m) - f(i - 1, j - 1, k - 1, m)) / dx[0];
    ret(1) = 0.25 * (f(i, j, k, m) - f(i, j - 1, k, m) + f(i - 1, j, k, m) - f(i - 1, j - 1, k, m) + f(i, j, k - 1, m) - f(i, j - 1, k - 1, m) + f(i - 1, j, k - 1, m) - f(i - 1, j - 1, k - 1, m)) / dx[1];
    ret(2) = 0.25 * (f(i, j, k, m) - f(i, j, k - 1, m) + f(i - 1, j, k, m) - f(i - 1, j, k - 1, m) + f(i, j - 1, k, m) - f(i, j - 1, k - 1, m) + f(i - 1, j - 1, k, m) - f(i - 1, j - 1, k - 1, m)) / dx[2];
#endif
    return ret;
}


template <class T>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
std::array<T, AMREX_SPACEDIM>
CellGradientOnNode(const amrex::Array4<const T>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    std::array<T, AMREX_SPACEDIM> ret;
    Util::Abort(INFO);
#if AMREX_SPACEDIM == 1
    ret[0] = (f(i, j, k, m) - f(i - 1, j, k, m)) / dx[0];
#elif AMREX_SPACEDIM == 2
    ret[0] = (f(i, j, k, m) - f(i - 1, j, k, m) + f(i, j - 1, k, m) - f(i - 1, j - 1, k, m)) * 0.5 / dx[0];
    ret[1] = (f(i, j, k, m) - f(i, j - 1, k, m) + f(i - 1, j, k, m) - f(i - 1, j - 1, k, m)) * 0.5 / dx[1];
#elif AMREX_SPACEDIM == 3
    ret[0] = (f(i, j, k, m) - f(i - 1, j, k, m) + f(i, j - 1, k, m) - f(i - 1, j - 1, k, m) + f(i, j, k - 1, m) - f(i - 1, j, k - 1, m) + f(i, j - 1, k - 1, m) - f(i - 1, j - 1, k - 1, m)) * 0.25 / dx[0];
    ret[1] = (f(i, j, k, m) - f(i, j - 1, k, m) + f(i - 1, j, k, m) - f(i - 1, j - 1, k, m) + f(i, j, k - 1, m) - f(i, j - 1, k - 1, m) + f(i - 1, j, k - 1, m) - f(i - 1, j - 1, k - 1, m)) * 0.25 / dx[1];
    ret[2] = (f(i, j, k, m) - f(i, j, k - 1, m) + f(i - 1, j, k, m) - f(i - 1, j, k - 1, m) + f(i, j - 1, k, m) - f(i, j - 1, k - 1, m) + f(i - 1, j - 1, k, m) - f(i - 1, j - 1, k - 1, m)) * 0.25 / dx[2];
#endif
    return ret;
}



[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
Gradient(   const amrex::Array4<const Set::Scalar>& f,
            const int& i, const int& j, const int& k,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix ret;
    ret(0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil));
#if AMREX_SPACEDIM > 1
    ret(0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 0, dx, stencil));
    ret(1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 1, dx, stencil));
    ret(1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 1, dx, stencil));
#if AMREX_SPACEDIM > 2
    ret(0, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 0, dx, stencil));
    ret(2, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 2, dx, stencil));
    ret(1, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 1, dx, stencil));
    ret(2, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 2, dx, stencil));
    ret(2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 2, dx, stencil));
#endif
#endif
    return ret;
}

/// \brief Forward displacement gradient used by the paired elastic stencil.
///
/// At a high physical boundary (`StencilType::Lo`) this falls back to the
/// backward one-sided difference so callers can safely use it in boundary
/// paths.  Interior and low-boundary directions use the forward difference.
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
ForwardGradient(const amrex::Array4<const Set::Scalar>& f,
                const int& i, const int& j, const int& k,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix ret;

    AMREX_D_TERM(
        ret(0, 0) = (stencil[0] == StencilType::Lo ? f(i, j, k, 0) - f(i - 1, j, k, 0) : f(i + 1, j, k, 0) - f(i, j, k, 0)) / dx[0];,
        ret(0, 1) = (stencil[1] == StencilType::Lo ? f(i, j, k, 0) - f(i, j - 1, k, 0) : f(i, j + 1, k, 0) - f(i, j, k, 0)) / dx[1];
        ret(1, 0) = (stencil[0] == StencilType::Lo ? f(i, j, k, 1) - f(i - 1, j, k, 1) : f(i + 1, j, k, 1) - f(i, j, k, 1)) / dx[0];
        ret(1, 1) = (stencil[1] == StencilType::Lo ? f(i, j, k, 1) - f(i, j - 1, k, 1) : f(i, j + 1, k, 1) - f(i, j, k, 1)) / dx[1];,
        ret(0, 2) = (stencil[2] == StencilType::Lo ? f(i, j, k, 0) - f(i, j, k - 1, 0) : f(i, j, k + 1, 0) - f(i, j, k, 0)) / dx[2];
        ret(1, 2) = (stencil[2] == StencilType::Lo ? f(i, j, k, 1) - f(i, j, k - 1, 1) : f(i, j, k + 1, 1) - f(i, j, k, 1)) / dx[2];
        ret(2, 0) = (stencil[0] == StencilType::Lo ? f(i, j, k, 2) - f(i - 1, j, k, 2) : f(i + 1, j, k, 2) - f(i, j, k, 2)) / dx[0];
        ret(2, 1) = (stencil[1] == StencilType::Lo ? f(i, j, k, 2) - f(i, j - 1, k, 2) : f(i, j + 1, k, 2) - f(i, j, k, 2)) / dx[1];
        ret(2, 2) = (stencil[2] == StencilType::Lo ? f(i, j, k, 2) - f(i, j, k - 1, 2) : f(i, j, k + 1, 2) - f(i, j, k, 2)) / dx[2];);

    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
Gradient(   const amrex::Array4<const Set::Vector>& f,
            const int& i, const int& j, const int& k,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix ret;

#if AMREX_SPACEDIM > 0
    ret.col(0) = Numeric::Stencil<Set::Vector, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil);
#endif
#if AMREX_SPACEDIM > 1
    ret.col(1) = Numeric::Stencil<Set::Vector, 0, 1, 0>::D(f, i, j, k, 0, dx, stencil);
#endif
#if AMREX_SPACEDIM > 2
    ret.col(2) = Numeric::Stencil<Set::Vector, 0, 0, 1>::D(f, i, j, k, 0, dx, stencil);
#endif

    return ret;
}

/// \brief Forward displacement gradient used by the paired elastic stencil.
///
/// See the scalar-multifab overload above for the high-boundary fallback.
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
ForwardGradient(const amrex::Array4<const Set::Vector>& f,
                const int& i, const int& j, const int& k,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix ret;

#if AMREX_SPACEDIM > 0
    ret.col(0) = (stencil[0] == StencilType::Lo ? f(i, j, k, 0) - f(i - 1, j, k, 0) : f(i + 1, j, k, 0) - f(i, j, k, 0)) / dx[0];
#endif
#if AMREX_SPACEDIM > 1
    ret.col(1) = (stencil[1] == StencilType::Lo ? f(i, j, k, 0) - f(i, j - 1, k, 0) : f(i, j + 1, k, 0) - f(i, j, k, 0)) / dx[1];
#endif
#if AMREX_SPACEDIM > 2
    ret.col(2) = (stencil[2] == StencilType::Lo ? f(i, j, k, 0) - f(i, j, k - 1, 0) : f(i, j, k + 1, 0) - f(i, j, k, 0)) / dx[2];
#endif

    return ret;
}

/// \brief Elastic gradient with a conservative paired bulk stencil and a
/// traction-compatible physical-boundary closure.
///
/// In the interior this is the forward member of the `D- G+` pair.  On a
/// physical face it retains the one-sided normal derivative but uses centered
/// tangential derivatives.  That is the natural strong-traction closure: the
/// same closure is used for the nonlinear residual, tangent, and the boundary
/// flux consumed by the adjacent interior divergence row.
/// Periodic ghost-shell stencil edges are deliberately *not* physical faces;
/// their forward/backward fallback is retained so the periodic operator stays
/// translation invariant.
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
PairedGradient(const amrex::Array4<const Set::Scalar>& f,
               const int& i, const int& j, const int& k,
               const Set::Scalar dx[AMREX_SPACEDIM],
               std::array<StencilType, AMREX_SPACEDIM> stencil,
               const std::array<bool, AMREX_SPACEDIM>& periodic)
{
    auto physical_stencil = stencil;
    AMREX_D_TERM(
        if (periodic[0]) physical_stencil[0] = StencilType::Central;,
        if (periodic[1]) physical_stencil[1] = StencilType::Central;,
        if (periodic[2]) physical_stencil[2] = StencilType::Central;);

    const bool on_physical_boundary = AMREX_D_TERM(
        physical_stencil[0] != StencilType::Central,
        || physical_stencil[1] != StencilType::Central,
        || physical_stencil[2] != StencilType::Central);
    return on_physical_boundary ? Gradient(f, i, j, k, dx, physical_stencil)
                                : ForwardGradient(f, i, j, k, dx, stencil);
}

/// \brief `Set::Vector` overload of `PairedGradient`.
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
PairedGradient(const amrex::Array4<const Set::Vector>& f,
               const int& i, const int& j, const int& k,
               const Set::Scalar dx[AMREX_SPACEDIM],
               std::array<StencilType, AMREX_SPACEDIM> stencil,
               const std::array<bool, AMREX_SPACEDIM>& periodic)
{
    auto physical_stencil = stencil;
    AMREX_D_TERM(
        if (periodic[0]) physical_stencil[0] = StencilType::Central;,
        if (periodic[1]) physical_stencil[1] = StencilType::Central;,
        if (periodic[2]) physical_stencil[2] = StencilType::Central;);

    const bool on_physical_boundary = AMREX_D_TERM(
        physical_stencil[0] != StencilType::Central,
        || physical_stencil[1] != StencilType::Central,
        || physical_stencil[2] != StencilType::Central);
    return on_physical_boundary ? Gradient(f, i, j, k, dx, physical_stencil)
                                : ForwardGradient(f, i, j, k, dx, stencil);
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Jet<Set::Vector,Set::Matrix>
GradientSplit(  const amrex::Array4<const Set::Vector>& f,
                const int& i, const int& j, const int& k,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Vector diag;
    Set::Matrix offdiag;

    Set::Jet<Set::Scalar, Set::Vector> ret;
#if AMREX_SPACEDIM > 0
    ret = Numeric::Stencil<Set::Vector, 1, 0, 0>::Dsplit(f, i, j, k, 0, dx, stencil);
    diag(0)        = ret.zero;
    offdiag.col(0) = ret.first;
#endif
#if AMREX_SPACEDIM > 1
    ret = Numeric::Stencil<Set::Vector, 0, 1, 0>::Dsplit(f, i, j, k, 0, dx, stencil);
    diag(1)        = ret.zero;
    offdiag.col(1) = ret.first;
#endif
#if AMREX_SPACEDIM > 2
    ret = Numeric::Stencil<Set::Vector, 0, 0, 1>::Dsplit(f, i, j, k, 0, dx, stencil);
    diag(2)        = ret.zero;
    offdiag.col(2) = ret.first;
#endif
    return Set::Jet(diag, offdiag);
}


[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Vector
NodeGradientOnCell(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret;
#if AMREX_SPACEDIM == 1
    ret(0) = (f(i + 1, j, k) - f(i, j, k)) / dx[0];
#elif AMREX_SPACEDIM == 2
    ret(0) = (f(i + 1, j, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i, j + 1, k)) * 0.5 / dx[0];
    ret(1) = (f(i, j + 1, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i + 1, j, k)) * 0.5 / dx[1];
#elif AMREX_SPACEDIM == 3
    ret(0) = (f(i + 1, j, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i, j + 1, k) + f(i + 1, j, k + 1) - f(i, j, k + 1) + f(i + 1, j + 1, k + 1) - f(i, j + 1, k + 1)) * 0.25 / dx[0];
    ret(1) = (f(i, j + 1, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i + 1, j, k) + f(i, j + 1, k + 1) - f(i, j, k + 1) + f(i + 1, j + 1, k + 1) - f(i + 1, j, k + 1)) * 0.25 / dx[1];
    ret(2) = (f(i, j, k + 1) - f(i, j, k) + f(i, j + 1, k + 1) - f(i, j + 1, k) + f(i + 1, j, k + 1) - f(i + 1, j, k) + f(i + 1, j + 1, k + 1) - f(i + 1, j + 1, k)) * 0.25 / dx[2];
#endif
    return ret;
}

[[nodiscard]] AMREX_FORCE_INLINE
Set::Matrix
NodeGradientOnCell(const amrex::Array4<const Set::Vector>& f,
    const int& i, const int& j, const int& k,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix ret;
#if AMREX_SPACEDIM == 1
    ret.col(0) = (f(i + 1, j, k) - f(i, j, k)) / dx[0];
#elif AMREX_SPACEDIM == 2
    ret.col(0) = (f(i + 1, j, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i, j + 1, k)) * 0.5 / dx[0];
    ret.col(1) = (f(i, j + 1, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i + 1, j, k)) * 0.5 / dx[1];
#elif AMREX_SPACEDIM == 3
    ret.col(0) = (f(i + 1, j, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i, j + 1, k) + f(i + 1, j, k + 1) - f(i, j, k + 1) + f(i + 1, j + 1, k + 1) - f(i, j + 1, k + 1)) * 0.25 / dx[0];
    ret.col(1) = (f(i, j + 1, k) - f(i, j, k) + f(i + 1, j + 1, k) - f(i + 1, j, k) + f(i, j + 1, k + 1) - f(i, j, k + 1) + f(i + 1, j + 1, k + 1) - f(i + 1, j, k + 1)) * 0.25 / dx[1];
    ret.col(2) = (f(i, j, k + 1) - f(i, j, k) + f(i, j + 1, k + 1) - f(i, j + 1, k) + f(i + 1, j, k + 1) - f(i + 1, j, k) + f(i + 1, j + 1, k + 1) - f(i + 1, j + 1, k)) * 0.25 / dx[2];
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
NodeGradientOnCell(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret;
#if AMREX_SPACEDIM == 1
    ret(0) = (f(i + 1, j, k, m) - f(i, j, k, m)) / dx[0];
#elif AMREX_SPACEDIM == 2
    ret(0) = 0.5 * (f(i + 1, j + 1, k, m) - f(i, j + 1, k, m) + f(i + 1, j, k, m) - f(i, j, k, m)) / dx[0];
    ret(1) = 0.5 * (f(i + 1, j + 1, k, m) - f(i + 1, j, k, m) + f(i, j + 1, k, m) - f(i, j, k, m)) / dx[1];
#elif AMREX_SPACEDIM == 3
    ret(0) = 0.25 * (f(i + 1, j + 1, k + 1, m) - f(i, j + 1, k + 1, m) + f(i + 1, j, k + 1, m) - f(i, j, k + 1, m) + f(i + 1, j + 1, k, m) - f(i, j + 1, k, m) + f(i + 1, j, k, m) - f(i, j, k, m)) / dx[0];
    ret(1) = 0.25 * (f(i + 1, j + 1, k + 1, m) - f(i + 1, j, k + 1, m) + f(i, j + 1, k + 1, m) - f(i, j, k + 1, m) + f(i + 1, j + 1, k, m) - f(i + 1, j, k, m) + f(i, j + 1, k, m) - f(i, j, k, m)) / dx[1];
    ret(2) = 0.25 * (f(i + 1, j + 1, k + 1, m) - f(i + 1, j + 1, k, m) + f(i, j + 1, k + 1, m) - f(i, j + 1, k, m) + f(i + 1, j, k + 1, m) - f(i + 1, j, k, m) + f(i, j, k + 1, m) - f(i, j, k, m)) / dx[2];
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix3
Gradient(   const amrex::Array4<const Set::Matrix>& f,
            const int& i, const int& j, const int& k,
            const Set::Scalar dx[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix3 ret;

#if AMREX_SPACEDIM > 0
    ret[0] = Numeric::Stencil<Set::Matrix, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil);
#endif
#if AMREX_SPACEDIM > 1
    ret[1] = Numeric::Stencil<Set::Matrix, 0, 1, 0>::D(f, i, j, k, 0, dx, stencil);
#endif
#if AMREX_SPACEDIM > 2
    ret[2] = Numeric::Stencil<Set::Matrix, 0, 0, 1>::D(f, i, j, k, 0, dx, stencil);
#endif

    return ret;
}

template<class T>
[[nodiscard]] AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
std::array<T,AMREX_SPACEDIM>
Gradient_Diagonal(  const Set::Scalar dx[AMREX_SPACEDIM],
                    std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType());

template<>
[[nodiscard]] AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
std::array<Set::Matrix,AMREX_SPACEDIM>
Gradient_Diagonal(  const Set::Scalar dx[AMREX_SPACEDIM],
                    std::array<StencilType, AMREX_SPACEDIM> stencil)
{

    std::array<Set::Matrix,AMREX_SPACEDIM> gradu;
    AMREX_D_TERM(   bool xmin = stencil[0] == StencilType::Hi;
                    bool xmax = stencil[0] == StencilType::Lo;, 
                    bool ymin = stencil[1] == StencilType::Hi;
                    bool ymax = stencil[1] == StencilType::Lo;,
                    bool zmin = stencil[2] == StencilType::Hi;
                    bool zmax = stencil[2] == StencilType::Lo; );
        
    for (int p = 0; p < AMREX_SPACEDIM; p++)
        for (int q = 0; q < AMREX_SPACEDIM; q++)
        {
            AMREX_D_TERM(
                gradu[p](q, 0) = ((!xmax ? 0.0 : (p == q ? 1.0 : 0.0)) - (!xmin ? 0.0 : (p == q ? 1.0 : 0.0))) / ((xmin || xmax ? 1.0 : 2.0) * dx[0]);,
                gradu[p](q, 1) = ((!ymax ? 0.0 : (p == q ? 1.0 : 0.0)) - (!ymin ? 0.0 : (p == q ? 1.0 : 0.0))) / ((ymin || ymax ? 1.0 : 2.0) * dx[1]);,
                gradu[p](q, 2) = ((!zmax ? 0.0 : (p == q ? 1.0 : 0.0)) - (!zmin ? 0.0 : (p == q ? 1.0 : 0.0))) / ((zmin || zmax ? 1.0 : 2.0) * dx[2]););
        }
    return gradu;
}

template<>
[[nodiscard]] AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
std::array<Set::Matrix3,AMREX_SPACEDIM>
Gradient_Diagonal(  const Set::Scalar dx[AMREX_SPACEDIM],
                    std::array<StencilType, AMREX_SPACEDIM> stencil)
{
    AMREX_D_TERM(   Util::Assert(INFO,TEST(stencil[0]==StencilType::Central));,
                    Util::Assert(INFO,TEST(stencil[1]==StencilType::Central));,
                    Util::Assert(INFO,TEST(stencil[2]==StencilType::Central)););

    std::array<Set::Matrix3,AMREX_SPACEDIM> gradgradu;

    for (int p = 0; p < AMREX_SPACEDIM; p++)
        for (int q = 0; q < AMREX_SPACEDIM; q++)
        {

            AMREX_D_TERM(
                gradgradu[p](q, 0, 0) = (p == q ? -2.0 : 0.0) / dx[0] / dx[0];
                ,// 2D
                gradgradu[p](q, 0, 1) = 0.0;
                gradgradu[p](q, 1, 0) = 0.0;
                gradgradu[p](q, 1, 1) = (p == q ? -2.0 : 0.0) / dx[1] / dx[1];
                ,// 3D
                gradgradu[p](q, 0, 2) = 0.0;
                gradgradu[p](q, 1, 2) = 0.0;
                gradgradu[p](q, 2, 0) = 0.0;
                gradgradu[p](q, 2, 1) = 0.0;
                gradgradu[p](q, 2, 2) = (p == q ? -2.0 : 0.0) / dx[2] / dx[2]);
        }

    return gradgradu;
}



[[nodiscard]] AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Set::Matrix3
NodeGradientOnCell( const amrex::Array4<const Set::Matrix>& f,
                    const int& i, const int& j, const int& k,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix3 ret;

#if AMREX_SPACEDIM > 0
    ret[0] = (f(i+1,j,k) - f(i,j,k)) / dx[0];
#endif
#if AMREX_SPACEDIM > 1
    ret[1] = (f(i,j+1,k) - f(i,j,k)) / dx[1];
#endif
#if AMREX_SPACEDIM > 2
    ret[2] = (f(i,j,k+1) - f(i,j,k)) / dx[2];
#endif

    return ret;
}


[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix3
MatrixGradient( const amrex::Array4<const Set::Scalar>& f,
                const int& i, const int& j, const int& k,
                const Set::Scalar dx[AMREX_SPACEDIM],
                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix3 ret;
#if AMREX_SPACEDIM == 1
    ret[0](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil));
#elif AMREX_SPACEDIM == 2
    ret[0](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil));
    ret[0](0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 0, dx, stencil));
    ret[0](1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 1, dx, stencil));
    ret[0](1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 1, dx, stencil));
    ret[1](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 2, dx, stencil));
    ret[1](0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 2, dx, stencil));
    ret[1](1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 3, dx, stencil));
    ret[1](1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 3, dx, stencil));
#elif AMREX_SPACEDIM == 3
    ret[0](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 0, dx, stencil));
    ret[0](1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 1, dx, stencil));
    ret[0](2, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 2, dx, stencil));
    ret[1](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 3, dx, stencil));
    ret[1](1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 4, dx, stencil));
    ret[1](2, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 5, dx, stencil));
    ret[2](0, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 6, dx, stencil));
    ret[2](1, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 7, dx, stencil));
    ret[2](2, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(f, i, j, k, 8, dx, stencil));
    ret[0](0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 0, dx, stencil));
    ret[0](1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 1, dx, stencil));
    ret[0](2, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 2, dx, stencil));
    ret[1](0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 3, dx, stencil));
    ret[1](1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 4, dx, stencil));
    ret[1](2, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 5, dx, stencil));
    ret[2](0, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 6, dx, stencil));
    ret[2](1, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 7, dx, stencil));
    ret[2](2, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(f, i, j, k, 8, dx, stencil));
    ret[0](0, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 0, dx, stencil));
    ret[0](1, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 1, dx, stencil));
    ret[0](2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 2, dx, stencil));
    ret[1](0, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 3, dx, stencil));
    ret[1](1, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 4, dx, stencil));
    ret[1](2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 5, dx, stencil));
    ret[2](0, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 6, dx, stencil));
    ret[2](1, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 7, dx, stencil));
    ret[2](2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(f, i, j, k, 8, dx, stencil));
#endif
    return ret;
}


[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
Hessian(const amrex::Array4<const Set::Scalar>& f,
        const int& i, const int& j, const int& k, const int& m,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix ret;
    ret(0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, m, dx, stencil));
#if AMREX_SPACEDIM > 1
    ret(1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, m, dx, stencil));
    ret(0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, m, dx, stencil));
    ret(1, 0) = ret(0, 1);
#if AMREX_SPACEDIM > 2
    ret(2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, m, dx, stencil));
    ret(1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, m, dx, stencil));
    ret(2, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, m, dx, stencil));
    ret(2, 1) = ret(1, 2);
    ret(0, 2) = ret(2, 0);
#endif
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix3
Hessian(const amrex::Array4<const Set::Scalar>& f,
        const int& i, const int& j, const int& k,
        const Set::Scalar DX[AMREX_SPACEDIM],
        std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix3 ret;
    // 1D 
#if AMREX_SPACEDIM>0
    ret(0, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 0, DX, stencil));
#endif
    // 2D
#if AMREX_SPACEDIM>1
    ret(0, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 0, DX, stencil));
    ret(0, 1, 0) = ret(0, 0, 1);
    ret(0, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 0, DX, stencil));
    ret(1, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 1, DX, stencil));
    ret(1, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 1, DX, stencil));
    ret(1, 1, 0) = ret(1, 0, 1);
    ret(1, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 1, DX, stencil));
#endif 
    // 3D
#if AMREX_SPACEDIM>2
    ret(0, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 0, DX, stencil));
    ret(0, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 0, DX, stencil));
    ret(0, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 0, DX, stencil));
    ret(0, 2, 0) = ret(0, 0, 2);
    ret(0, 2, 1) = ret(0, 1, 2);;
    ret(1, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 1, DX, stencil));
    ret(1, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 1, DX, stencil));
    ret(1, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 1, DX, stencil));
    ret(1, 2, 0) = ret(1, 0, 2);
    ret(1, 2, 1) = ret(1, 1, 2);;
    ret(2, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 1, 0) = ret(2, 0, 1);
    ret(2, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(f, i, j, k, 2, DX, stencil));
    ret(2, 2, 0) = ret(2, 0, 2);
    ret(2, 2, 1) = ret(2, 1, 2);;
#endif
    return ret;
}


// Returns Hessian of a vector field.
// Return value: ret[i](j,k) = ret_{i,jk}
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix3
Hessian(const amrex::Array4<const Set::Vector>& f,
        const int& i, const int& j, const int& k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
{
    Set::Matrix3 ret;

#if AMREX_SPACEDIM>0    
    Set::Vector f_11 = Numeric::Stencil<Set::Vector, 2, 0, 0>::D(f, i, j, k, 0, dx, stencil);
    ret[0](0, 0) = f_11(0);
#endif
#if AMREX_SPACEDIM>1
    ret[1](0, 0) = f_11(1);
    Set::Vector f_12 = Numeric::Stencil<Set::Vector, 1, 1, 0>::D(f, i, j, k, 0, dx, stencil);
    ret[0](1, 0) = ret[0](0, 1) = f_12(0);
    ret[1](1, 0) = ret[1](0, 1) = f_12(1);
    Set::Vector f_22 = Numeric::Stencil<Set::Vector, 0, 2, 0>::D(f, i, j, k, 0, dx, stencil);
    ret[0](1, 1) = f_22(0);
    ret[1](1, 1) = f_22(1);
#endif
#if AMREX_SPACEDIM>2
    ret[2](0, 0) = f_11(2);
    ret[2](1, 0) = ret[2](0, 1) = f_12(2);
    Set::Vector f_13 = Numeric::Stencil<Set::Vector, 1, 0, 1>::D(f, i, j, k, 0, dx, stencil);
    ret[0](2, 0) = ret[0](0, 2) = f_13(0);
    ret[1](2, 0) = ret[1](0, 2) = f_13(1);
    ret[2](2, 0) = ret[2](0, 2) = f_13(2);
    ret[2](1, 1) = f_22(2);
    Set::Vector f_23 = Numeric::Stencil<Set::Vector, 0, 1, 1>::D(f, i, j, k, 0, dx, stencil);
    ret[0](1, 2) = ret[0](2, 1) = f_23(0);
    ret[1](1, 2) = ret[1](2, 1) = f_23(1);
    ret[2](1, 2) = ret[2](2, 1) = f_23(2);
    Set::Vector f_33 = Numeric::Stencil<Set::Vector, 0, 0, 2>::D(f, i, j, k, 0, dx, stencil);
    ret[0](2, 2) = f_33(0);
    ret[1](2, 2) = f_33(1);
    ret[2](2, 2) = f_33(2);
#endif
    return ret;
}


[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
FieldToMatrix(  const amrex::Array4<const Set::Scalar>& f,
                const int& i, const int& j, const int& k)
{
    Set::Matrix ret;
#if AMREX_SPACEDIM == 1
    ret(0, 0) = f(i, j, k, 0);

#elif AMREX_SPACEDIM == 2
    ret(0, 0) = f(i, j, k, 0); ret(0, 1) = f(i, j, k, 1);
    ret(1, 0) = f(i, j, k, 2); ret(1, 1) = f(i, j, k, 3);

#elif AMREX_SPACEDIM == 3
    ret(0, 0) = f(i, j, k, 0); ret(0, 1) = f(i, j, k, 1); ret(0, 2) = f(i, j, k, 2);
    ret(1, 0) = f(i, j, k, 3); ret(1, 1) = f(i, j, k, 4); ret(1, 2) = f(i, j, k, 5);
    ret(2, 0) = f(i, j, k, 6); ret(2, 1) = f(i, j, k, 7); ret(2, 2) = f(i, j, k, 8);
#endif

    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
FieldToMatrix(  const amrex::Array4<Set::Scalar>& f,
                const int& i, const int& j, const int& k)
{
    Set::Matrix ret;
#if AMREX_SPACEDIM == 1
    ret(0, 0) = f(i, j, k, 0);

#elif AMREX_SPACEDIM == 2
    ret(0, 0) = f(i, j, k, 0); ret(0, 1) = f(i, j, k, 1);
    ret(1, 0) = f(i, j, k, 2); ret(1, 1) = f(i, j, k, 3);

#elif AMREX_SPACEDIM == 3
    ret(0, 0) = f(i, j, k, 0); ret(0, 1) = f(i, j, k, 1); ret(0, 2) = f(i, j, k, 2);
    ret(1, 0) = f(i, j, k, 3); ret(1, 1) = f(i, j, k, 4); ret(1, 2) = f(i, j, k, 5);
    ret(2, 0) = f(i, j, k, 6); ret(2, 1) = f(i, j, k, 7); ret(2, 2) = f(i, j, k, 8);
#endif

    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
FieldToVector(  const amrex::Array4<const Set::Scalar>& f,
                const int& i, const int& j, const int& k)
{
    Set::Vector ret;
    ret(0) = f(i, j, k, 0);
#if AMREX_SPACEDIM > 1
    ret(1) = f(i, j, k, 1);
#if AMREX_SPACEDIM > 2
    ret(2) = f(i, j, k, 2);
#endif
#endif
    return ret;
}

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
FieldToVector(  const amrex::Array4<Set::Scalar>& f,
                const int& i, const int& j, const int& k)
{
    Set::Vector ret;
    ret(0) = f(i, j, k, 0);
#if AMREX_SPACEDIM > 1
    ret(1) = f(i, j, k, 1);
#if AMREX_SPACEDIM > 2
    ret(2) = f(i, j, k, 2);
#endif
#endif
    return ret;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
MatrixToField(  const amrex::Array4<Set::Scalar>& f,
                const int& i, const int& j, const int& k,
                Set::Matrix matrix)
{
#if AMREX_SPACEDIM == 1
    f(i, j, k, 0) = matrix(0, 0);
#elif AMREX_SPACEDIM == 2
    f(i, j, k, 0) = matrix(0, 0); f(i, j, k, 1) = matrix(0, 1);
    f(i, j, k, 2) = matrix(1, 0); f(i, j, k, 3) = matrix(1, 1);
#elif AMREX_SPACEDIM == 3
    f(i, j, k, 0) = matrix(0, 0); f(i, j, k, 1) = matrix(0, 1); f(i, j, k, 2) = matrix(0, 2);
    f(i, j, k, 3) = matrix(1, 0); f(i, j, k, 4) = matrix(1, 1); f(i, j, k, 5) = matrix(1, 2);
    f(i, j, k, 6) = matrix(2, 0); f(i, j, k, 7) = matrix(2, 1); f(i, j, k, 8) = matrix(2, 2);
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
VectorToField(  const amrex::Array4<Set::Scalar>& f,
                const int& i, const int& j, const int& k,
                Set::Vector vector)
{
    f(i, j, k, 0) = vector(0);
#if AMREX_SPACEDIM > 1
    f(i, j, k, 1) = vector(1);
#if AMREX_SPACEDIM > 2
    f(i, j, k, 2) = vector(2);
#endif
#endif
}

template<int index, int SYM>
[[nodiscard]] 
AMREX_GPU_HOST_DEVICE
Set::Matrix3 
Divergence( const amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM, SYM>>&,
            const int, const int, const int, const Set::Scalar[AMREX_SPACEDIM],
            std::array<StencilType, AMREX_SPACEDIM> /*stecil*/ = DefaultType())
{
    Util::Abort(INFO, "Not implemented yet"); return Set::Matrix3::Zero();
}

template <>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix3 Divergence<2, Set::Sym::Isotropic>(
    const amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>>& C,
    const int i, const int j, const int k, const Set::Scalar dx[AMREX_SPACEDIM],
    std::array<StencilType, AMREX_SPACEDIM> stencil)
{
    Set::Matrix3 ret;

    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic> gradCx =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>, 1, 0, 0>::D(C, i, j, k, 0, dx, stencil);
#if AMREX_SPACEDIM>1
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic> gradCy =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>, 0, 1, 0>::D(C, i, j, k, 0, dx, stencil);
#endif
#if AMREX_SPACEDIM>2
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic> gradCz =
        Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>, 0, 0, 1>::D(C, i, j, k, 0, dx, stencil);
#endif

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int k = 0; k < AMREX_SPACEDIM; k++)
            for (int l = 0; l < AMREX_SPACEDIM; l++)
            {
                ret(i, k, l) = 0.0;
                ret(i, k, l) = gradCx(i, 0, k, l);
#if AMREX_SPACEDIM>1
                ret(i, k, l) = gradCy(i, 1, k, l);
#endif
#if AMREX_SPACEDIM>2
                ret(i, k, l) = gradCz(i, 2, k, l);
#endif
            }
    return ret;
}




template<int dim>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix4<dim, Set::Sym::Full>
DoubleHessian(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM]);

template<>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix4<2, Set::Sym::Full>
DoubleHessian<2>(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix4<2, Set::Sym::Full> ret;
    // [0,0,0,0]
    ret(0, 0, 0, 0) = Stencil<Set::Scalar, 4, 0, 0>::D(f, i, j, k, m, dx);
    // [0, 0, 0, 1]
    ret(0, 0, 0, 1) = Stencil<Set::Scalar, 3, 1, 0>::D(f, i, j, k, m, dx);
    // [0, 0, 1, 1]
    ret(0, 0, 1, 1) = Stencil<Set::Scalar, 2, 2, 0>::D(f, i, j, k, m, dx);
    // [0, 1, 1, 1]
    ret(0, 1, 1, 1) = Stencil<Set::Scalar, 1, 3, 0>::D(f, i, j, k, m, dx);
    // [1, 1, 1, 1]
    ret(1, 1, 1, 1) = Stencil<Set::Scalar, 0, 4, 0>::D(f, i, j, k, m, dx);
    return ret;
}

template<>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix4<3, Set::Sym::Full>
DoubleHessian<3>(const amrex::Array4<const Set::Scalar>& f,
    const int& i, const int& j, const int& k, const int& m,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix4<3, Set::Sym::Full> ret;
    // [0,0,0,0]
    ret(0, 0, 0, 0) = Stencil<Set::Scalar, 4, 0, 0>::D(f, i, j, k, m, dx);
    // [0, 0, 0, 1]
    ret(0, 0, 0, 1) = Stencil<Set::Scalar, 3, 1, 0>::D(f, i, j, k, m, dx);
    // [0, 0, 0, 2]
    ret(0, 0, 0, 2) = Stencil<Set::Scalar, 3, 0, 1>::D(f, i, j, k, m, dx);
    // [0, 0, 1, 1]
    ret(0, 0, 1, 1) = Stencil<Set::Scalar, 2, 2, 0>::D(f, i, j, k, m, dx);
    // [0, 0, 1, 2]
    ret(0, 0, 1, 2) = Stencil<Set::Scalar, 2, 1, 1>::D(f, i, j, k, m, dx);
    // [0, 0, 2, 2]
    ret(0, 0, 2, 2) = Stencil<Set::Scalar, 2, 0, 2>::D(f, i, j, k, m, dx);
    // [0, 1, 1, 1]
    ret(0, 1, 1, 1) = Stencil<Set::Scalar, 1, 3, 0>::D(f, i, j, k, m, dx);
    // [0, 1, 1, 2]
    ret(0, 1, 1, 2) = Stencil<Set::Scalar, 1, 2, 1>::D(f, i, j, k, m, dx);
    // [0, 1, 2, 2]
    ret(0, 1, 2, 2) = Stencil<Set::Scalar, 1, 1, 2>::D(f, i, j, k, m, dx);
    // [0, 2, 2, 2]
    ret(0, 2, 2, 2) = Stencil<Set::Scalar, 1, 0, 3>::D(f, i, j, k, m, dx);
    // [1, 1, 1, 1]
    ret(1, 1, 1, 1) = Stencil<Set::Scalar, 0, 4, 0>::D(f, i, j, k, m, dx);
    // [1, 1, 1, 2]
    ret(1, 1, 1, 2) = Stencil<Set::Scalar, 0, 3, 1>::D(f, i, j, k, m, dx);
    // [1, 1, 2, 2]
    ret(1, 1, 2, 2) = Stencil<Set::Scalar, 0, 2, 2>::D(f, i, j, k, m, dx);
    // [1, 2, 2, 2]
    ret(1, 2, 2, 2) = Stencil<Set::Scalar, 0, 1, 3>::D(f, i, j, k, m, dx);
    // [2, 2, 2, 2]
    ret(2, 2, 2, 2) = Stencil<Set::Scalar, 0, 0, 4>::D(f, i, j, k, m, dx);
    return ret;
}

struct Interpolate
{
public:
    template <class T>
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T CellToNodeAverage( const amrex::Array4<const T>& f,
                                const int& i, const int& j, const int& k, const int& m,
                                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        int AMREX_D_DECL(
            ilo = (stencil[0] == Numeric::StencilType::Lo ? 0 : 1),
            jlo = (stencil[1] == Numeric::StencilType::Lo ? 0 : 1),
            klo = (stencil[2] == Numeric::StencilType::Lo ? 0 : 1));

        return (AMREX_D_TERM(f(i, j, k, m) + f(i - ilo, j, k, m)
            ,
            +f(i, j - jlo, k, m) + f(i - ilo, j - jlo, k, m)
            ,
            +f(i, j, k - klo, m) + f(i - ilo, j, k - klo, m)
            + f(i, j - jlo, k - klo, m) + f(i - ilo, j - jlo, k - klo, m)
        )) * AMREX_D_PICK(0.5, 0.25, 0.125);
    }

    template <class T>
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T CellToNodeAverage( const amrex::Array4<T>& f,
                                const int& i, const int& j, const int& k, const int& m,
                                std::array<StencilType, AMREX_SPACEDIM> stencil = DefaultType())
    {
        int AMREX_D_DECL(
            ilo = (stencil[0] == Numeric::StencilType::Lo ? 0 : 1),
            jlo = (stencil[1] == Numeric::StencilType::Lo ? 0 : 1),
            klo = (stencil[2] == Numeric::StencilType::Lo ? 0 : 1));

        return (AMREX_D_TERM(f(i, j, k, m) + f(i - ilo, j, k, m)
            ,
            +f(i, j - jlo, k, m) + f(i - ilo, j - jlo, k, m)
            ,
            +f(i, j, k - klo, m) + f(i - ilo, j, k - klo, m)
            + f(i, j - jlo, k - klo, m) + f(i - ilo, j - jlo, k - klo, m)
        )) * AMREX_D_PICK(0.5, 0.25, 0.125);
    }

    /// \brief Average the distinct physical cells adjacent to a node.
    ///
    /// Unlike `CellToNodeAverage`, this never reads a non-periodic physical
    /// ghost cell.  At a physical face/edge/corner its clamped duplicate cell
    /// is counted once, while periodic directions use the already-filled
    /// neighboring ghost cell.  This gives the paired elastic residual,
    /// tangent, and diagonal one decomposition-independent raw-psi policy.
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static Set::Scalar CellToNodeAverageBounded(const amrex::Array4<const Set::Scalar>& f,
                                                const int& i, const int& j, const int& k, const int& m,
                                                const amrex::Dim3& cell_lo, const amrex::Dim3& cell_hi,
                                                const std::array<bool, AMREX_SPACEDIM>& periodic)
    {
        int i0 = i - 1, i1 = i;
#if AMREX_SPACEDIM > 1
        int j0 = j - 1, j1 = j;
#endif
#if AMREX_SPACEDIM > 2
        int k0 = k - 1, k1 = k;
#endif

        if (!periodic[0])
        {
            i0 = i0 < cell_lo.x ? cell_lo.x : (i0 > cell_hi.x ? cell_hi.x : i0);
            i1 = i1 < cell_lo.x ? cell_lo.x : (i1 > cell_hi.x ? cell_hi.x : i1);
        }
#if AMREX_SPACEDIM > 1
        if (!periodic[1])
        {
            j0 = j0 < cell_lo.y ? cell_lo.y : (j0 > cell_hi.y ? cell_hi.y : j0);
            j1 = j1 < cell_lo.y ? cell_lo.y : (j1 > cell_hi.y ? cell_hi.y : j1);
        }
#endif
#if AMREX_SPACEDIM > 2
        if (!periodic[2])
        {
            k0 = k0 < cell_lo.z ? cell_lo.z : (k0 > cell_hi.z ? cell_hi.z : k0);
            k1 = k1 < cell_lo.z ? cell_lo.z : (k1 > cell_hi.z ? cell_hi.z : k1);
        }
#endif

        const int nx = i0 == i1 ? 1 : 2;
#if AMREX_SPACEDIM > 1
        const int ny = j0 == j1 ? 1 : 2;
#endif
#if AMREX_SPACEDIM > 2
        const int nz = k0 == k1 ? 1 : 2;
#endif

        Set::Scalar sum = 0.0;
        for (int ii = 0; ii < nx; ++ii)
        {
            const int ix = ii == 0 ? i0 : i1;
#if AMREX_SPACEDIM == 1
            sum += f(ix, j, k, m);
#elif AMREX_SPACEDIM == 2
            for (int jj = 0; jj < ny; ++jj)
            {
                const int iy = jj == 0 ? j0 : j1;
                sum += f(ix, iy, k, m);
            }
#else
            for (int jj = 0; jj < ny; ++jj)
            {
                const int iy = jj == 0 ? j0 : j1;
                for (int kk = 0; kk < nz; ++kk)
                {
                    const int iz = kk == 0 ? k0 : k1;
                    sum += f(ix, iy, iz, m);
                }
            }
#endif
        }

        const int count = nx
#if AMREX_SPACEDIM > 1
            * ny
#endif
#if AMREX_SPACEDIM > 2
            * nz
#endif
            ;
        return sum / static_cast<Set::Scalar>(count);
    }

    template <class T>
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T NodeToCellAverage( const amrex::Array4<const T>& f,
                                const int& i, const int& j, const int& k, const int& m)
    {
        return (AMREX_D_TERM(f(i, j, k, m) + f(i + 1, j, k, m)
            ,
            +f(i, j + 1, k, m) + f(i + 1, j + 1, k, m)
            ,
            +f(i, j, k + 1, m) + f(i + 1, j, k + 1, m)
            + f(i, j + 1, k + 1, m) + f(i + 1, j + 1, k + 1, m)
        )) * AMREX_D_PICK(0.5, 0.25, 0.125);
    }
    template <class T>
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static T NodeToCellAverage( const amrex::Array4<T>& f,
                                const int& i, const int& j, const int& k, const int& m)
    {
        return (AMREX_D_TERM(f(i, j, k, m) + f(i + 1, j, k, m)
            ,
            +f(i, j + 1, k, m) + f(i + 1, j + 1, k, m)
            ,
            +f(i, j, k + 1, m) + f(i + 1, j, k + 1, m)
            + f(i, j + 1, k + 1, m) + f(i + 1, j + 1, k + 1, m)
        )) * AMREX_D_PICK(0.5, 0.25, 0.125);
    }
};

}
#endif
