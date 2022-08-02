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

public:
    enum Type {Displacement, Traction, Periodic, Neumann}; 

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


    virtual
    Set::Vector operator () (const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) = 0;

protected:
    Set::Scalar m_time = 0.0;
};
}
}
}
#endif
