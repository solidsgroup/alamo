#ifndef MODEL_SOLID_PLASTIC_J2_H_
#define MODEL_SOLID_PLASTIC_J2_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class J2Plastic : public Isotropic
{
public:
    J2Plastic() {};
    
    void DefineSlipSystem(Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta)
    {
        theta = a_theta;
        yield_strength = a_yield;
        hardening_modulus = a_hardening;
    }

    Set::Scalar YieldSurface()
    {
        return yield_strength + alpha*hardening_modulus;
    }
    
    Set::Scalar PlasticEnergy()
    {
        return (yield_strength*alpha + 0.5*hardening_modulus*alpha*alpha);
    }

    void Advance(Set::Scalar /*dt*/, Set::Matrix eps, Set::Matrix /*sig*/) override
    {
        Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));
        
        Set::Matrix sig = DW(strain_prev);
        Set::Matrix sigdev = sig - (1.0/((double)AMREX_SPACEDIM))*sig.trace()*Set::Matrix::Identity();

        Set::Matrix del_eps = eps - strain_prev;
        Set::Matrix epsdev = del_eps - (1.0/((double)AMREX_SPACEDIM))*del_eps.trace()*Set::Matrix::Identity();

        Set::Matrix zeta_trial = sigdev  - beta + 2.0*ddw.Mu()*epsdev;
        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*alpha);
        if( f_trial <= 0.0) return;
    
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(ddw.Mu())*(1.0 + (hardening_modulus/(3.0*ddw.Mu()))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        beta += SQ2O3*dH*n_new;
        F0 += dGamma*n_new;
        alpha += SQ2O3*dGamma;

        strain_prev = eps;

        if(std::isnan(F0.norm()))
        {
            Util::Message(INFO, "f_trial = ", f_trial);
            Util::Abort(INFO, "Nans detected");
        }
    }

public:
    Set::Matrix beta = Set::Matrix::Zero();
    Set::Matrix strain_prev = Set::Matrix::Zero();
    Set::Scalar alpha = 0.0;
    Set::Scalar theta;                  // isotropic and kinematic hardening parameter
    Set::Scalar yield_strength;          // yield strength
    Set::Scalar hardening_modulus;      // hardening modulus

public:
    static J2Plastic Zero()
    {
        J2Plastic ret;
        ret.Define(0.,0.,Set::Matrix::Zero());
        ret.DefineSlipSystem(0.,0.,0.);
        return ret;
    }
    static J2Plastic Random()
    {
        J2Plastic ret;
        ret.Define(Util::Random(),Util::Random(),Set::Matrix::Random());
        ret.DefineSlipSystem(Util::Random(),Util::Random(), Util::Random());
        ret.beta = Set::Matrix::Random();
        ret.strain_prev = Set::Matrix::Random();
        ret.alpha = Util::Random();
        return ret;
    }
    static void Parse(J2Plastic & value, IO::ParmParse & pp)
    {
        Isotropic::Parse(value,pp);
        
        Set::Scalar yield = 1.0, hardening = 1.0, theta_tmp = 1.0;
        pp.query("yield", yield);
        pp.query("hardening", hardening);
        pp.query("theta", theta_tmp);

        value.DefineSlipSystem(yield, hardening, theta_tmp);
    }
    #define OP_CLASS J2Plastic
    #define OP_VARS X(ddw) X(F0) X(beta) X(theta) X(yield_strength) X(hardening_modulus) X(strain_prev) X(alpha)
    #include "Model/Solid/InClassOperators.H"

};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

template<>
inline int Set::Field<Model::Solid::Affine::J2Plastic>::NComp() const 
{
    return AMREX_SPACEDIM*AMREX_SPACEDIM;
}

template<>
inline std::string Set::Field<Model::Solid::Affine::J2Plastic>::Name(int i) const 
{

    if (i==0) return name + ".epsp_xx";
    if (i==1) return name + ".epsp_xy";
    if (i==2) return name + AMREX_D_PICK("",".epsp_yx",".epsp_xz");
    if (i==3) return name + AMREX_D_PICK("",".epsp_yy",".epsp_yx");
    if (i==4) return name + ".epsp_yy";
    if (i==5) return name + ".epsp_yz";
    if (i==6) return name + ".epsp_zx";
    if (i==7) return name + ".epsp_zy";
    if (i==8) return name + ".epsp_zz";
    return name;
}

template<>
inline void Set::Field<Model::Solid::Affine::J2Plastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::J2Plastic> const & src = ((*this)[a_lev])->array(mfi);
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
