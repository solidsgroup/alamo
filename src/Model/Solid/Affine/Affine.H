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
//class Affine
//{
//public:
    struct PlasticState{
        Set::Matrix epsp;   // Plastic strain
        Set::Scalar alpha;  // Equivalent plastic strain
        Set::Matrix beta;

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

    /*struct{
        bool on;        // slip system active or not
        Set::Vector n;  // slip plane normal
        Set::Vector s;  // slip direction
        Set::Scalar tcrss;      // critically resolved shear stress
        Set::Scalar tcrss0;     // pristine material
        Set::Scalar exp;        // power law exponent
        Set::Scalar gam0dot;    // gamma0 dot.
    } slipsystem;*/
//};
}
}
}
#endif