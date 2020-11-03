#ifndef MODEL_SOLID_AFFINE_H_
#define MODEL_SOLID_AFFINE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
    enum LatticeType {FCC, BCC, HCP};   // standard Bravais lattices. More to be added later
//class Affine
//{

struct PlasticState{
    Set::Matrix epsp;   // Plastic strain
    Set::Scalar alpha;  // Equivalent plastic strain
    Set::Matrix beta;

    static PlasticState Zero()
    {
        PlasticState ret;
        ret.epsp = Set::Matrix::Zero();
        ret.beta = Set::Matrix::Zero();
        ret.alpha = 0.0;
        return ret;
    }

//    AMREX_FORCE_INLINE
//    void operator = (const PlasticState &rhs)
//    {
//        epsp = rhs.epsp;
//        beta = rhs.beta;
//        alpha = rhs.alpha;
//    }

    AMREX_FORCE_INLINE
    void operator += (const PlasticState &rhs)
    {
        epsp += rhs.epsp;
        alpha += rhs.alpha;
        beta += rhs.beta;
    }
    AMREX_FORCE_INLINE
    PlasticState operator * (const Set::Scalar _alpha) const
    {
        PlasticState ret;
        ret.epsp = epsp*_alpha;
        ret.alpha = alpha*_alpha;
        ret.beta = beta*_alpha;
        return ret;
    }
    friend PlasticState operator * (const Set::Scalar alpha, const PlasticState b);
    friend PlasticState operator + (const PlasticState a, const PlasticState b);
    friend PlasticState operator - (const PlasticState a, const PlasticState b);
};

AMREX_FORCE_INLINE
PlasticState operator * (const Set::Scalar alpha, const PlasticState b)
{
    PlasticState ret;
    ret.epsp = b.epsp*alpha;
    ret.beta = b.beta*alpha;
    ret.alpha = b.alpha*alpha;
    return ret;
}

AMREX_FORCE_INLINE
PlasticState operator + (const PlasticState a, const PlasticState b)
{
    PlasticState ret;
    ret.epsp = a.epsp + b.epsp;
    ret.beta = a.beta + b.beta;
    ret.alpha = a.alpha + b.alpha;
    return ret;
}

AMREX_FORCE_INLINE
PlasticState operator - (const PlasticState a, const PlasticState b)
{
    PlasticState ret;
    ret.epsp = a.epsp - b.epsp;
    ret.beta = a.beta - b.beta;
    ret.alpha = a.alpha - b.alpha;
    return ret;
}


struct SlipSystem{
    bool on;        // slip system active or not
    Set::Vector n;  // slip plane normal
    Set::Vector s;  // slip direction
    Set::Scalar tcrss;      // critically resolved shear stress
    Set::Scalar tcrss0;     // pristine material
    Set::Scalar s_exp;        // power law exponent
    Set::Scalar gam0dot;    // gamma0 dot.

//    AMREX_FORCE_INLINE
//    void operator = (const SlipSystem &rhs)
//    {
//        n = rhs.n;
//        s = rhs.s;
//        tcrss = rhs.tcrss;
//        tcrss0 = rhs.tcrss0;
//        s_exp = rhs.s_exp;
//        gam0dot = rhs.gam0dot;
//    }

    AMREX_FORCE_INLINE
    void operator += (const SlipSystem &rhs)
    {
        n += rhs.n;
        s += rhs.s;
        tcrss += rhs.tcrss;
        tcrss0 += rhs.tcrss0;
        s_exp += rhs.s_exp;
        gam0dot += rhs.gam0dot;
    }
    AMREX_FORCE_INLINE
    SlipSystem operator * (const Set::Scalar alpha) const
    {
        SlipSystem ret;
        ret.n = n*alpha;
        ret.s = s*alpha;
        ret.tcrss = tcrss*alpha;
        ret.tcrss0 = tcrss0*alpha;
        ret.s_exp = s_exp*alpha;
        ret.gam0dot = gam0dot*alpha;
        
        return ret;
    }
    friend SlipSystem operator * (const Set::Scalar alpha, const SlipSystem b);
    friend SlipSystem operator + (const SlipSystem a, const SlipSystem b);
    friend SlipSystem operator - (const SlipSystem a, const SlipSystem b);
};
AMREX_FORCE_INLINE
SlipSystem operator * (const Set::Scalar alpha, const SlipSystem b)
{
    SlipSystem ret;
    ret.n = b.n*alpha;
    ret.s = b.s*alpha;
    ret.tcrss = b.tcrss*alpha;
    ret.tcrss0 = b.tcrss0*alpha;
    ret.s_exp = b.s_exp*alpha;
    ret.gam0dot = b.gam0dot*alpha;
    
    return ret;
}
AMREX_FORCE_INLINE
SlipSystem operator + (const SlipSystem a, const SlipSystem b)
{
    SlipSystem ret;
    ret.n = a.n + b.n;
    ret.s = a.s + b.s;
    ret.tcrss = a.tcrss + b.tcrss;
    ret.tcrss0 = a.tcrss0 + b.tcrss0;
    ret.s_exp = a.s_exp + b.s_exp;
    ret.gam0dot = a.gam0dot + b.gam0dot;
    return ret;
}
AMREX_FORCE_INLINE
SlipSystem operator - (const SlipSystem a, const SlipSystem b)
{
    SlipSystem ret;
    ret.n = a.n - b.n;
    ret.s = a.s - b.s;
    ret.tcrss = a.tcrss - b.tcrss;
    ret.tcrss0 = a.tcrss0 - b.tcrss0;
    ret.s_exp = a.s_exp - b.s_exp;
    ret.gam0dot = a.gam0dot - b.gam0dot;
    return ret;
}
    
//};
}
}
}
#endif