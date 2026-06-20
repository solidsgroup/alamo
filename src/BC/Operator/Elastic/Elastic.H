//
// Class of BC operators that work with :ref:`Operator::Elastic`.
//

#ifndef BC_OPERATOR_ELASTIC_H
#define BC_OPERATOR_ELASTIC_H

// #include "Operator/Elastic.H"
#include "IO/ParmParse.H"
#include "Numeric/Interpolator/Linear.H"

namespace BC
{
namespace Operator
{
namespace Elastic
{
class Elastic
{
public:
    static const int min = 1, max = 2;

    virtual ~Elastic() = default;

public:
    enum Type {Displacement, Traction, Periodic, Neumann, None};

    #if AMREX_SPACEDIM==2
    enum Face{
        XLO, YLO, XHI, YHI,
        XLO_YLO, XLO_YHI, XHI_YLO, XHI_YHI,
        INT
    };
    #elif AMREX_SPACEDIM==3
    enum Face{
        XLO, YLO, ZLO, XHI, YHI, ZHI, // 6
        YLO_ZLO, YLO_ZHI, YHI_ZLO, YHI_ZHI, // 12
        ZLO_XLO, ZLO_XHI, ZHI_XLO, ZHI_XHI,
        XLO_YLO, XLO_YHI, XHI_YLO, XHI_YHI,
        XLO_YLO_ZLO, XLO_YLO_ZHI, XLO_YHI_ZLO, XLO_YHI_ZHI, // 8
        XHI_YLO_ZLO, XHI_YLO_ZHI, XHI_YHI_ZLO, XHI_YHI_ZHI,
        INT
    };
    #endif

    enum Direction {AMREX_D_DECL(X=0,Y=1,Z=2)};

    #if AMREX_SPACEDIM==2
    static const int m_nfaces = 8;
    #elif AMREX_SPACEDIM==3
    static const int m_nfaces = 26;
    #endif

    // The boundary-condition *type* per (face,direction). This is the only state
    // the in-kernel boundary operator needs, and it is plain-old-data (an array of
    // Type enums), so it can be captured by value into a device kernel. The
    // value data (Constant::m_bc_val interpolators, Expression::m_bc_func parsers)
    // lives in the derived classes and is only consumed by Init() on the host.
    using BcTypeArray = std::array<std::array<Type,AMREX_SPACEDIM>, m_nfaces>;

    void
    SetTime(const Set::Scalar a_time) {m_time = a_time;}

    virtual void
    Init(amrex::MultiFab * a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const = 0;

    virtual void
    Init(amrex::FabArray<amrex::BaseFab<Set::Vector>> * a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const = 0;

    void
    Init(Set::Field<Set::Scalar> &a_rhs,
        const amrex::Vector<amrex::Geometry> &a_geom,
        bool a_homogeneous = false) const
    {
        for (int ilev = 0; ilev <= a_rhs.finest_level; ilev++)
            Init(a_rhs[ilev].get(),a_geom[ilev],a_homogeneous);
    }

    void
    Init(Set::Field<Set::Vector> &a_rhs,
        const amrex::Vector<amrex::Geometry> &a_geom,
        bool a_homogeneous = false) const
    {
        for (int ilev = 0; ilev <= a_rhs.finest_level; ilev++)
            Init(a_rhs[ilev].get(),a_geom[ilev],a_homogeneous);
    }

#define SQRT3INV 0.57735026919
#define SQRT2INV 0.70710678118

    // Device-copyable accessor: copy this on the host and capture it into a
    // device ParallelFor, then call eval()/typeOf() below.
    AMREX_GPU_HOST_DEVICE
    const BcTypeArray& GetBcTypeArray() const { return m_bc_type; }

    // ----------------------------------------------------------------------
    // Static, device-callable boundary logic. These are shared by every BC
    // type (Constant/TensionTest/Expression) because they all dispatch purely
    // on the (face,direction) Type. Operating on a passed-in m_bc_type avoids
    // dereferencing a host object/vtable on the device.
    // ----------------------------------------------------------------------

    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    static Set::Vector set(const std::array<Type,AMREX_SPACEDIM> &bc_type,
                    const Set::Vector &u, const Set::Matrix &gradu, const Set::Matrix &sigma, Set::Vector n)
    {
        Set::Vector f = Set::Vector::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            if      (bc_type[i] == Type::Displacement)
                f(i) = u(i);
            else if (bc_type[i] == Type::Traction)
                f(i) = (sigma*n)(i);
            else if (bc_type[i] == Type::Neumann)
                f(i) = (gradu*n)(i);
            else
                continue;
        }
        return f;
    }

    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    static std::array<Type,AMREX_SPACEDIM> typeOf (
                const BcTypeArray &m_bc_type,
                const int &i, const int &j, [[maybe_unused]] const int &k,
                const amrex::Box &domain)
    {
        amrex::IntVect m(AMREX_D_DECL(i,j,k));
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

        // Corners
        #if AMREX_SPACEDIM == 2
        if (m[0] == lo.x && m[1] == lo.y)                  return m_bc_type[Face::XLO_YLO];
        if (m[0] == lo.x && m[1] == hi.y)                  return m_bc_type[Face::XLO_YHI];
        if (m[0] == hi.x && m[1] == lo.y)                  return m_bc_type[Face::XHI_YLO];
        if (m[0] == hi.x && m[1] == hi.y)                  return m_bc_type[Face::XHI_YHI];
        if (m[0] == lo.x)                                  return m_bc_type[Face::XLO];
        if (m[0] == hi.x)                                  return m_bc_type[Face::XHI];
        if (m[1] == lo.y)                                  return m_bc_type[Face::YLO];
        if (m[1] == hi.y)                                  return m_bc_type[Face::YHI];

        #elif AMREX_SPACEDIM == 3

        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return m_bc_type[Face::XLO_YLO_ZLO];
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return m_bc_type[Face::XLO_YLO_ZHI];
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return m_bc_type[Face::XLO_YHI_ZLO];
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return m_bc_type[Face::XLO_YHI_ZHI];
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return m_bc_type[Face::XHI_YLO_ZLO];
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return m_bc_type[Face::XHI_YLO_ZHI];
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return m_bc_type[Face::XHI_YHI_ZLO];
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return m_bc_type[Face::XHI_YHI_ZHI];

        if (m[1] == lo.y && m[2] == lo.z)                  return m_bc_type[Face::YLO_ZLO];
        if (m[1] == lo.y && m[2] == hi.z)                  return m_bc_type[Face::YLO_ZHI];
        if (m[1] == hi.y && m[2] == lo.z)                  return m_bc_type[Face::YHI_ZLO];
        if (m[1] == hi.y && m[2] == hi.z)                  return m_bc_type[Face::YHI_ZHI];
        if (m[2] == lo.z && m[0] == lo.x)                  return m_bc_type[Face::ZLO_XLO];
        if (m[2] == lo.z && m[0] == hi.x)                  return m_bc_type[Face::ZLO_XHI];
        if (m[2] == hi.z && m[0] == lo.x)                  return m_bc_type[Face::ZHI_XLO];
        if (m[2] == hi.z && m[0] == hi.x)                  return m_bc_type[Face::ZHI_XHI];
        if (m[0] == lo.x && m[1] == lo.y)                  return m_bc_type[Face::XLO_YLO];
        if (m[0] == lo.x && m[1] == hi.y)                  return m_bc_type[Face::XLO_YHI];
        if (m[0] == hi.x && m[1] == lo.y)                  return m_bc_type[Face::XHI_YLO];
        if (m[0] == hi.x && m[1] == hi.y)                  return m_bc_type[Face::XHI_YHI];

        if (m[0] == lo.x)                                  return m_bc_type[Face::XLO];
        if (m[0] == hi.x)                                  return m_bc_type[Face::XHI];
        if (m[1] == lo.y)                                  return m_bc_type[Face::YLO];
        if (m[1] == hi.y)                                  return m_bc_type[Face::YHI];
        if (m[2] == lo.z)                                  return m_bc_type[Face::ZLO];
        if (m[2] == hi.z)                                  return m_bc_type[Face::ZHI];

        #endif

        return {AMREX_D_DECL(Type::None,Type::None,Type::None)};
    }

    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    static Set::Vector eval (const BcTypeArray &m_bc_type,
                const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain)
    {
        (void)i; (void)j; (void)k; // Suppress "unused variable" warnings

        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

        amrex::IntVect m(AMREX_D_DECL(i,j,k));

        // Corners
        #if AMREX_SPACEDIM == 2

        if (m[0] == lo.x && m[1] == lo.y)                  return set(m_bc_type[Face::XLO_YLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(m_bc_type[Face::XLO_YHI],     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(m_bc_type[Face::XHI_YLO],     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(m_bc_type[Face::XHI_YHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV));

        if (m[0] == lo.x)                                  return set(m_bc_type[Face::XLO],         u, gradu, sigma, Set::Vector(-1, 0));
        if (m[0] == hi.x)                                  return set(m_bc_type[Face::XHI],         u, gradu, sigma, Set::Vector(+1, 0));
        if (m[1] == lo.y)                                  return set(m_bc_type[Face::YLO],         u, gradu, sigma, Set::Vector( 0,-1));
        if (m[1] == hi.y)                                  return set(m_bc_type[Face::YHI],         u, gradu, sigma, Set::Vector( 0,+1));

        #elif AMREX_SPACEDIM == 3

        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YLO_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YLO_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YHI_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YHI_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YLO_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YLO_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YHI_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YHI_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,+SQRT3INV));

        if (m[1] == lo.y && m[2] == lo.z)                  return set(m_bc_type[Face::YLO_ZLO],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,-SQRT2INV));
        if (m[1] == lo.y && m[2] == hi.z)                  return set(m_bc_type[Face::YLO_ZHI],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,+SQRT2INV));
        if (m[1] == hi.y && m[2] == lo.z)                  return set(m_bc_type[Face::YHI_ZLO],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,-SQRT2INV));
        if (m[1] == hi.y && m[2] == hi.z)                  return set(m_bc_type[Face::YHI_ZHI],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,+SQRT2INV));
        if (m[2] == lo.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZLO_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == lo.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZLO_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == hi.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZHI_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       +SQRT2INV));
        if (m[2] == hi.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZHI_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       +SQRT2INV));
        if (m[0] == lo.x && m[1] == lo.y)                  return set(m_bc_type[Face::XLO_YLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(m_bc_type[Face::XLO_YHI],     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(m_bc_type[Face::XHI_YLO],     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(m_bc_type[Face::XHI_YHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV,0       ));

        if (m[0] == lo.x)                                  return set(m_bc_type[Face::XLO],         u, gradu, sigma, Set::Vector(-1, 0, 0));
        if (m[0] == hi.x)                                  return set(m_bc_type[Face::XHI],         u, gradu, sigma, Set::Vector(+1, 0, 0));
        if (m[1] == lo.y)                                  return set(m_bc_type[Face::YLO],         u, gradu, sigma, Set::Vector( 0,-1, 0));
        if (m[1] == hi.y)                                  return set(m_bc_type[Face::YHI],         u, gradu, sigma, Set::Vector( 0,+1, 0));
        if (m[2] == lo.z)                                  return set(m_bc_type[Face::ZLO],         u, gradu, sigma, Set::Vector( 0, 0,-1));
        if (m[2] == hi.z)                                  return set(m_bc_type[Face::ZHI],         u, gradu, sigma, Set::Vector( 0, 0,+1));

        #endif

        Util::Abort(INFO,"Boundary condition error");
        return Set::Vector::Zero();
    }

#ifndef ALAMO_GPU
    virtual
    AMREX_GPU_HOST_DEVICE
    std::array<Type,AMREX_SPACEDIM> getType (
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) = 0;

    virtual
    AMREX_GPU_HOST_DEVICE
    Set::Vector operator () (const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) = 0;
#else
    // Non-static convenience wrappers for GPU call sites that hold a concrete
    // BC object without requiring virtual dispatch.
    AMREX_GPU_HOST_DEVICE
    std::array<Type,AMREX_SPACEDIM> getType (
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) const
    {
        return typeOf(m_bc_type, i, j, k, domain);
    }

    AMREX_GPU_HOST_DEVICE
    Set::Vector operator () (const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) const
    {
        return eval(m_bc_type, u, gradu, sigma, i, j, k, domain);
    }
#endif

protected:
    Set::Scalar m_time = 0.0;
    BcTypeArray m_bc_type;
};
}
}
}
#endif
