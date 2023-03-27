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
public:
    void Randomize() {}
    static Laplacian Zero()
    {
        Laplacian ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>::Zero();
        return ret;
    }
    static void Parse(Laplacian,IO::ParmParse) {}
    
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal> ddw;
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    static Laplacian Random()
    {
        Laplacian ret;
        return ret;
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
