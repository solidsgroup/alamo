#ifndef MODEL_SOLID_FINITE_PSEUDOLINEARPREDEFORMED_H_
#define MODEL_SOLID_FINITE_PSEUDOLINEARPREDEFORMED_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Finite/PseudoLinearCubic.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
class PseudoLinearCubicPredeformed : public PseudoLinearCubic
{
public:
    PseudoLinearCubicPredeformed() {};
    PseudoLinearCubicPredeformed(PseudoLinearCubic base) : PseudoLinearCubic(base) {};
    virtual ~PseudoLinearCubicPredeformed() {};

    Set::Scalar W(const Set::Matrix& F) const override
    {
        Set::Matrix F0inv = F0.inverse();
        return PseudoLinearCubic::W(F * F0inv);
    }
    Set::Matrix DW(const Set::Matrix& F) const override
    {
        Set::Matrix F0inv = F0.inverse();
        return PseudoLinearCubic::DW(F * F0inv) * F0inv.transpose();
    }
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& F) const override
    {
        Set::Matrix F0inv = F0.inverse();
        auto DDW = PseudoLinearCubic::DDW(F * F0inv);
        auto ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            for (int J = 0; J < AMREX_SPACEDIM; J++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    for (int L = 0; L < AMREX_SPACEDIM; L++)
                    {
                        ret(i,J,k,L) = 0.0;
                        for (int Q = 0; Q < AMREX_SPACEDIM; Q++)
                            for (int R = 0; R < AMREX_SPACEDIM; R++)
                                ret(i,J,k,L) += DDW(i,Q,k,R)*F0inv(J,Q)*F0inv(L,R);
                    }
        return ret;
    }

public:
    Set::Matrix F0 = Set::Matrix::Identity();

public:
    static PseudoLinearCubicPredeformed Zero()
    {
        PseudoLinearCubicPredeformed ret = PseudoLinearCubic::Zero();
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static PseudoLinearCubicPredeformed Random()
    {
        PseudoLinearCubicPredeformed ret = PseudoLinearCubic::Random();
        ret.F0 = Set::Matrix::Random();
        return ret;
    }

    AMREX_FORCE_INLINE
    static PseudoLinearCubicPredeformed Combine(const std::vector<PseudoLinearCubicPredeformed> &models, const std::vector<Set::Scalar> &eta, int order)
    {
        PseudoLinearCubicPredeformed ret = PseudoLinearCubicPredeformed::Zero();
        if (order == 1)
        {
            Set::Scalar etasum = 0.;
            for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
            for (unsigned int n = 0 ; n < models.size(); n++)
            {
                ret.C11 += models[n].C11 * (eta[n] / etasum);
                ret.C12 += models[n].C12 * (eta[n] / etasum);
                ret.C44 += models[n].C44 * (eta[n] / etasum);
                ret.q   += models[n].q   * (eta[n] / etasum);
                ret.F0  += models[n].F0  * (eta[n] / etasum);
            }
            return ret;
        }
        else if (order == 2)
        {
            Set::Scalar etasum = 0.;
            for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n]*eta[n];
            for (unsigned int n = 0 ; n < models.size(); n++)
            {
                ret.C11 += models[n].C11 * (eta[n]*eta[n] / etasum);
                ret.C12 += models[n].C12 * (eta[n]*eta[n] / etasum);
                ret.C44 += models[n].C44 * (eta[n]*eta[n] / etasum);
                ret.q   += models[n].q   * (eta[n]*eta[n] / etasum);
                ret.F0  += models[n].F0  * (eta[n]*eta[n] / etasum);
            }
            return ret;
        }

        Util::Exception(INFO,"invalid value used for order, ",order);
        return ret; // should never happen
    }


    static void Parse(PseudoLinearCubicPredeformed& value, IO::ParmParse& pp)
    {
        PseudoLinearCubic::Parse(value, pp);
        if (pp.contains("eps0") && pp.contains("F0"))
        {
            Util::Abort("Cannot specify both F0 and eps0");
        }
        else if (pp.contains("F0")) 
        {
            Set::Matrix F0;
            pp_queryarr("F0", F0); // Large-strain eigendeformation (Identity = no deformation)
            value.F0 = F0;
        }
        else if (pp.contains("eps0"))
        {
            Set::Matrix eps0;
            pp_queryarr("eps0",eps0); // Small-strain eigendeformation (Zero = no deformation)
            value.F0 = eps0 + Set::Matrix::Identity();
        }
        else
        {
            value.F0 = Set::Matrix::Identity();
        }
        Util::Assert(INFO,TEST(fabs(value.F0.determinant()) > 1E-8 ),"F0 must be non-singular");
    }

#define OP_CLASS PseudoLinearCubicPredeformed
#define OP_VARS X(C11) X(C12) X(C44) X(q) X(F0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}


template<>
ALAMO_SINGLE_DEFINITION
int Set::Field<Model::Solid::Finite::PseudoLinearCubicPredeformed>::NComp() const
{
    return 4;
}

template<>
ALAMO_SINGLE_DEFINITION
std::string Set::Field<Model::Solid::Finite::PseudoLinearCubicPredeformed>::Name(int i) const
{
    if (i == 0) return name + "_F0xx";
    if (i == 1) return name + "_F0xy";
    if (i == 2) return name + "_F0yx";
    if (i == 3) return name + "_F0yy";
    return name;
}

template<>
ALAMO_SINGLE_DEFINITION
void Set::Field<Model::Solid::Finite::PseudoLinearCubicPredeformed>::Copy(int a_lev, amrex::MultiFab& a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Finite::PseudoLinearCubicPredeformed> const& src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                dst(i, j, k, a_dstcomp + 0) = src(i, j, k).F0(0, 0);
                dst(i, j, k, a_dstcomp + 1) = src(i, j, k).F0(0, 1);
                dst(i, j, k, a_dstcomp + 2) = src(i, j, k).F0(1, 0);
                dst(i, j, k, a_dstcomp + 3) = src(i, j, k).F0(1, 1);
            });
        }
    }
}


#endif
