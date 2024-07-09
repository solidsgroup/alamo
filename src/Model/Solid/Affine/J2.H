//
// This models an isotropic elastic-perfectly-plastic, non-time-dependent solid model.
// 
// The energy and derivatives are:
//
// .. math::
//    :nowrap:
//
//    \begin{gather}
//    W = \frac{1}{2}(\varepsilon - \varepsilon_p):\mathbb{C}(\varepsilon-\varepsilon_p) \\
//    DW = \mathbb{C}(\varepsilon-\varepsilon_p) \\
//    DDW = \mathbb{C}
//    \end{gather}
//
// where :math:`\mathbb{C}` is an isotropic :code:`Set::Matrix4` and :math:`\varepsilon_p` is 
// is stored in the :code:`F0` eigenstrain.
//
// The plastic strain is evolved according to the following:
//
// #. Calculate the deviatoric stress :math:`\sigma_v=\sigma - \frac{1}{3}tr(\sigma)\mathbf{I}`
// #. Calculate :math:`J_2=\sqrt{\frac{3}{2}\sigma_v:\sigma_v}`
// #. If :math:`J_2<\sigma_0` then quit - no plasticity occurs
// #. Calculate :math:`\Delta\sigma = (1-\frac{\sigma_0}{J_2})`, which projects the stress 
//    back on the yield surface.
// #. Convert to change in plastic strain, :math:`\Delta\varepsilon=\mathbb{C}^{-1}\Delta\sigma`
// #. Update plastic strain: :math:`\varepsilon_p += \Delta\varepsilon`
//
// Notes:
//
// * This does not implement any kind of hardening model. Rate hardening, isotropic hardening,
//   and kinematic hardening have yet to be implemneted.
//
#ifndef MODEL_SOLID_AFFINE_J2_H_
#define MODEL_SOLID_AFFINE_J2_H_

#include "AMReX.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class J2 : public Isotropic
{
public:

    J2() {};

    void Advance(Set::Scalar dt, Set::Matrix eps, Set::Matrix /*sig*/) override
    {
        Set::Matrix sigma = DW(eps);
        Set::Matrix sig_dev = sigma - sigma.trace() * Set::Matrix::Identity() / (double)AMREX_SPACEDIM;
        if ((hardening<0.0) || (ratecoeff<0.0)) // rate independent
        {
            Set::Scalar J2 = sqrt(1.5 * (sig_dev.transpose() * sig_dev).trace());
            if (J2 < sigma0) return; // No plasticity happens.
            Set::Matrix dsigp = (1. - sigma0/J2) * sig_dev;
            F0 += ddw.Inverse() * dsigp; 
        }
        else // rate-dependent
        {
            Set::Matrix eps_dev = eps - eps.trace() * Set::Matrix::Identity() / (double)AMREX_SPACEDIM;
            Set::Scalar m_mu = ddw.Mu();
            Set::Matrix zeta_trial = sig_dev + 2.0*m_mu*eps_dev;
            Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(sigma0 + hardening*alpha);
            if (f_trial <= 0.0) return;
            Set::Matrix n_new = zeta_trial / zeta_trial.norm();
            Set::Scalar dGamma = f_trial / (2.0*m_mu*(1 + hardening/(3.0*m_mu)));
            dGamma *= ratecoeff * dt;
            alpha = alpha + SQ2O3*dGamma;
            F0 += dGamma * n_new;
        }
    }

    // Parameters
    Set::Scalar sigma0 = 1.0;
    Set::Scalar hardening = NAN;
    Set::Scalar ratecoeff = NAN;

    // Internal variables to track state
    Set::Scalar alpha = 0.0;
    static constexpr Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));

    static J2 Zero()
    {
        J2 ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>::Zero();
        ret.F0 = Set::Matrix::Zero();
        ret.sigma0 = 0.0;
        ret.hardening = 0.0;
        ret.ratecoeff = 0.0;
        return ret;
    }

    // See also inputs to :ref:`Model::Solid::Affine::Isotropic`
    static void Parse(J2 & value, IO::ParmParse & pp)
    {
        Isotropic::Parse(value,pp);
        // J2 Yield criterion
        pp_query_default("sigma0",value.sigma0,1.0);
        pp_query_default("hardening",value.hardening,-1.0);
        pp_query_default("ratecoeff",value.ratecoeff,-1.0);
    }

    #define OP_CLASS J2
    #define OP_VARS  X(ddw) X(F0) X(sigma0) X(hardening) X(ratecoeff)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}


template<>
AMREX_ATTRIBUTE_WEAK
int Set::Field<Model::Solid::Affine::J2>::NComp() const 
{
    return AMREX_SPACEDIM*AMREX_SPACEDIM;
}

template<>
AMREX_ATTRIBUTE_WEAK
std::string Set::Field<Model::Solid::Affine::J2>::Name(int i) const 
{
    if (i==0) return name + ".epsp_xx";
    if (i==1) return name + ".epsp_xy";
    if (i==2) return name + ".epsp_xz";
    if (i==3) return name + ".epsp_yx";
    if (i==4) return name + ".epsp_yy";
    if (i==5) return name + ".epsp_yz";
    if (i==6) return name + ".epsp_zx";
    if (i==7) return name + ".epsp_zy";
    if (i==8) return name + ".epsp_zz";
    return name;
}

template<>
AMREX_ATTRIBUTE_WEAK
void Set::Field<Model::Solid::Affine::J2>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::J2> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k).F0(n/AMREX_SPACEDIM,n%AMREX_SPACEDIM);
                });
            }
        }
    }    
}


#endif
