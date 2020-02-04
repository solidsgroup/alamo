#ifndef MODEL_SOLID_LINEAR_LAPLACIAN_H_
#define MODEL_SOLID_LINEAR_LAPLACIAN_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"

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
    Set::Scalar W(Set::Matrix gradu)
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(Set::Matrix gradu)
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal> DDW(Set::Matrix /*gradu*/)
    {
        return ddw;
    }
	using Model::Solid::Solid<Set::Sym::Diagonal>::operator();
public:
	void Randomize() {}
    static const KinematicVariable kinvar = KinematicVariable::gradu;
	static Laplacian Random()
    {
        Laplacian ret;
        return ret;
    }

};

}
}
}
#endif


