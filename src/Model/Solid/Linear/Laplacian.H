#ifndef MODEL_SOLID_LINEAR_LAPLACIAN_H_
#define MODEL_SOLID_LINEAR_LAPLACIAN_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "IO/ParmParse.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
///
/// Isotropic implements the following elastic model
/// \f[ \mathbb{C}_{ijkl} = \alpha\delta_{ik}\delta_{jl}  \f]
/// corresponding to the elastic model
/// \f[ \sigma_{ij} = \alpha\delta_{ik}\delta_{jl}u_{k,l} = \alpha u_{i,j} \f]
/// \f[ f_{i} = \mathbb{C}_{ijkl}u_{k,jl} = \alpha\delta_{ik}\delta_{jl}u_{k,jl} = \alpha u_{i,jj} \f]
///
class Laplacian : public Model::Solid::Solid<Set::Sym::Diagonal>
{
public:
    Laplacian() { Define(); } ;
    Laplacian(Solid<Set::Sym::Diagonal> base) : Solid<Set::Sym::Diagonal>(base) {};
    ~Laplacian() {} ;
    void Define()
    {
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>(Set::Matrix::Ones());
    }
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return ddw;
    }
    using Model::Solid::Solid<Set::Sym::Diagonal>::operator();
public:
    void Randomize() {}
    static Laplacian Zero()
    {
        Laplacian ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>::Zero();
        return ret;
    }
    static void Parse(Laplacian,IO::ParmParse) {}
    
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    static Laplacian Random()
    {
        Laplacian ret;
        return ret;
    }
    AMREX_FORCE_INLINE
    void operator += (const Laplacian &rhs)
    {
        ddw += rhs.ddw;
    }    
    AMREX_FORCE_INLINE
        virtual Laplacian operator * (const Set::Scalar alpha) const
    {
        Laplacian ret;
        ret.ddw = ddw*alpha;
        return ret;
    }    
    friend Laplacian operator * (const Set::Scalar alpha, const Laplacian b);
    friend Laplacian operator + (const Laplacian a, const Laplacian b);
    friend Laplacian operator - (const Laplacian a, const Laplacian b);
};
AMREX_FORCE_INLINE
Laplacian operator * (const Set::Scalar alpha, const Laplacian b)
{
    Laplacian ret;
    ret.ddw = b.ddw*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
Laplacian operator + (const Laplacian a, const Laplacian b)
{
    Laplacian ret;
    ret.ddw = a.ddw + b.ddw;
    return ret;
}    
AMREX_FORCE_INLINE
Laplacian operator - (const Laplacian a, const Laplacian b)
{
    Laplacian ret;
    ret.ddw = a.ddw - b.ddw;
    return ret;
}    


}
}
}
#endif
