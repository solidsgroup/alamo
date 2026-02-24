#ifndef MODEL_SOLID_FINITE_NEOHOOKEANPREDEFORMED_H_
#define MODEL_SOLID_FINITE_NEOHOOKEANPREDEFORMED_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Finite/NeoHookean.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
class NeoHookeanPredeformed : public NeoHookean
{
public:
    NeoHookeanPredeformed() {};
    NeoHookeanPredeformed(NeoHookean base) : NeoHookean(base) {};
    virtual ~NeoHookeanPredeformed() {};

    Set::Scalar W(const Set::Matrix& F) const override
    {
        return NeoHookean::W(F * F0.inverse());
    }
    Set::Matrix DW(const Set::Matrix& F) const override
    {
        return NeoHookean::DW(F * F0.inverse());
    }
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& F) const override
    {
        return NeoHookean::DDW(F * F0.inverse());
    }
    virtual void Print(std::ostream& out) const override
    {
        out << "mu = " << mu << " kappa = " << kappa << " F0 = " << F0;
    }

public:
    Set::Matrix F0 = Set::Matrix::Identity();

public:
    static NeoHookeanPredeformed Zero()
    {
        NeoHookeanPredeformed ret = NeoHookean::Zero();
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static NeoHookean Random()
    {
        NeoHookeanPredeformed ret = NeoHookean::Random();
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(NeoHookeanPredeformed& value, IO::ParmParse& pp)
    {
        pp.queryclass<NeoHookean>(value);
        
        if (pp.contains("eps0") && pp.contains("F0"))
        {
            Util::Abort("Cannot specify both F0 and eps0");
        }
        else if (pp.contains("F0")) 
        {
            Set::Matrix F0;
            pp_queryarr("F0", F0); // Large-deformation eigendeformation (Identity = no deformation)
            value.F0 = F0;
        }
        else if (pp.contains("eps0"))
        {
            Set::Matrix eps0;
            pp_queryarr("eps0",eps0); // Small-deformation eigendeformation (Zero = no deformation)
            value.F0 = eps0 + Set::Matrix::Identity(); 
        }
        else
        {
            value.F0 = Set::Matrix::Identity();
        }
        Util::Assert(INFO,TEST(fabs(value.F0.determinant()) > 1E-8 ),"F0 must be non-singular");
    }

#define OP_CLASS NeoHookeanPredeformed
#define OP_VARS X(kappa) X(mu) X(F0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}



template<>
ALAMO_SINGLE_DEFINITION
int Set::Field<Model::Solid::Finite::NeoHookeanPredeformed>::NComp() const
{
    return 2 + AMREX_SPACEDIM * AMREX_SPACEDIM;
}

template<>
ALAMO_SINGLE_DEFINITION
std::string Set::Field<Model::Solid::Finite::NeoHookeanPredeformed>::Name(int i) const
{
    if (i == 0) return name + "_mu";
    if (i == 1) return name + "_kappa";
#if AMREX_SPACEDIM==2
    if (i == 2) return name + "_F0xx";
    if (i == 3) return name + "_F0xy";
    if (i == 4) return name + "_F0yx";
    if (i == 5) return name + "_F0yy";
#elif AMREX_SPACEDIM==3
    //Util::Abort(INFO, "Not implemented yet");
    if (i == 2) return name + "_F0xx";
    if (i == 3) return name + "_F0xy";
    if (i == 4) return name + "_F0yx";
    if (i == 5) return name + "_F0yy";
#endif
    return name;
}

template<>
ALAMO_SINGLE_DEFINITION
void Set::Field<Model::Solid::Finite::NeoHookeanPredeformed>::Copy(int a_lev, amrex::MultiFab& a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Finite::NeoHookeanPredeformed> const& src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                dst(i, j, k, a_dstcomp + 0) = src(i, j, k).mu;
                dst(i, j, k, a_dstcomp + 1) = src(i, j, k).kappa;
#if AMREX_SPACEDIM==2
                dst(i, j, k, a_dstcomp + 2) = src(i, j, k).F0(0, 0);
                dst(i, j, k, a_dstcomp + 3) = src(i, j, k).F0(0, 1);
                dst(i, j, k, a_dstcomp + 4) = src(i, j, k).F0(1, 0);
                dst(i, j, k, a_dstcomp + 5) = src(i, j, k).F0(1, 1);
#elif AMREX_SPACEDIM==3
                dst(i, j, k, a_dstcomp + 2) = src(i, j, k).F0(0, 0);
                dst(i, j, k, a_dstcomp + 3) = src(i, j, k).F0(0, 1);
                dst(i, j, k, a_dstcomp + 4) = src(i, j, k).F0(1, 0);
                dst(i, j, k, a_dstcomp + 5) = src(i, j, k).F0(1, 1);
                //Util::Abort(INFO, "Not implemented");
#endif

                //dst(i, j, k, a_dstcomp + 0) = src(i, j, k).ddw(0, 0, 0, 0);
            });
        }
    }
}


#endif
