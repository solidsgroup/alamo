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
    Laplacian(Set::Scalar alpha) { Define(alpha); } ;
    Laplacian(Solid<Set::Sym::Diagonal> base) : Solid<Set::Sym::Diagonal>(base) {};
    ~Laplacian() {} ;
    void Define(Set::Scalar alpha = 1.0)
    {
        ddw = alpha*Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>(Set::Matrix::Ones());
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
public:
    static Laplacian Zero()
    {
        return Laplacian(0.0);
    }
    static void Parse(Laplacian value,IO::ParmParse pp)
    {
        Set::Scalar alpha = NAN;
        // Coefficient for the Laplacian
        pp_query_default("alpha",alpha,1.0);
        value.Define(alpha);
    }
    
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal> ddw;
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    static Laplacian Random()
    {
        Set::Scalar alpha = (Util::Random() + 0.1)/1.1;
        return Laplacian(alpha);
    }

#define OP_CLASS Laplacian
#define OP_VARS X(ddw)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}
#endif
