//
// This model implements an affine hexagonal material:
//
// .. math::
//
//    W_{aff}(\nabla\mathbf{u}) = W_{lin}(\nabla\mathbf{u} - \mathbf{F}_0)
//
// where :math:`W_{aff}` is the present model, :math:`\mathbf{F}_0` is the eigenstrain,
// and :math:`W_{lin}` is the cubic model defined in :ref:`Model::Solid::Linear::Hexagonal`.
//
#ifndef MODEL_SOLID_AFFINE_HEXAGONAL_H_
#define MODEL_SOLID_AFFINE_HEXAGONAL_H_

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/Hexagonal.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Hexagonal : public Linear::Hexagonal
{
public:

    Hexagonal() {};
    Hexagonal(Linear::Hexagonal base) : Linear::Hexagonal(base) {};
    virtual ~Hexagonal() {};

    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2,Set::Matrix a_F0=Set::Matrix::Zero())
    {
        Linear::Hexagonal::Define(C11,C12,C13,C33,C44,phi1,Phi,phi2);
        F0 = a_F0;
    }
    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, Eigen::Matrix3d R, Set::Matrix a_F0=Set::Matrix::Zero())
    {
        Linear::Hexagonal::Define(C11,C12,C13,C33,C44,R);
        F0 = a_F0;
    }
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return Linear::Hexagonal::W(gradu - F0);
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return Linear::Hexagonal::DW(gradu - F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & gradu) const override
    {
        return Linear::Hexagonal::DDW(gradu - F0);
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << ddw;
    }
    AMREX_FORCE_INLINE
    void SetF0 (Set::Matrix &a_F0) {F0 = a_F0;}

public:
    Set::Matrix F0 = Set::Matrix::Zero();
    static const KinematicVariable kinvar = KinematicVariable::gradu;

    static Hexagonal Zero()
    {
        Hexagonal ret = Linear::Hexagonal::Zero();
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static Hexagonal Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random(),Util::Random(),Util::Random());
    }
    static Hexagonal Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44)
    {
        Hexagonal ret = Linear::Hexagonal::Random(C11,C12,C13,C33,C44);
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(Hexagonal & value, IO::ParmParse & pp)
    {
        Linear::Hexagonal::Parse(value,pp);
        pp.queryarr("F0",value.F0); // Eigenstrain matrix. Can be defined in 2D or 3D.
    }

    AMREX_FORCE_INLINE
    static Hexagonal Combine(const std::vector<Hexagonal> &models, const std::vector<Set::Scalar> &eta)
    {
        Hexagonal ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        ret.F0 = Set::Matrix::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret.ddw += models[n].ddw * (eta[n] / etasum);
            ret.F0  += models[n].F0  * (eta[n] / etasum);
        }
        return ret;
    }

    #define OP_CLASS Hexagonal
    #define OP_VARS  X(ddw) X(F0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

template<>
AMREX_ATTRIBUTE_WEAK
int Set::Field<Model::Solid::Affine::Hexagonal>::NComp() const 
{
    return 1;
}

template<>
AMREX_ATTRIBUTE_WEAK
std::string Set::Field<Model::Solid::Affine::Hexagonal>::Name(int i) const 
{
    if (i==0) return name + ".Cxxxx";
    return name;
}

template<>
AMREX_ATTRIBUTE_WEAK
void Set::Field<Model::Solid::Affine::Hexagonal>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::Hexagonal> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                dst(i,j,k,a_dstcomp + 0) = src(i,j,k).ddw(0,0,0,0);
            });
        }
    }    
}


#endif
