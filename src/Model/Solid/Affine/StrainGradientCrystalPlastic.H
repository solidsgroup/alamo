#ifndef MODEL_SOLID_PLASTIC_STRAINGRADIENTCRYSTALPLASTIC_H_
#define MODEL_SOLID_PLASTIC_STRAINGRADIENTCRYSTALPLASTIC_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Cubic.H"

/*
    Note: everything is hardcoded for FCC right now. 
    We will deal with templating BCC and HCP later.
*/

namespace Model
{
namespace Solid
{
namespace Affine
{
class StrainGradientCrystalPlastic : public Cubic
{
public:
    StrainGradientCrystalPlastic(){};

    //StrainGradientCrystalPlastic(StrainGradientCrystalPlastic &alpha){}

    void DefineSlipGradientSystem()
    {
        // m_Ke = a_Ke;
        // m_b = a_b;
        // m_de = a_de;
        // m_ds = a_ds;
        // m_tau0 = a_tau0;
        // m_mu = a_mu;
        // m_c0 = a_c0;
        // m_w1 = a_w1;
        // m_w2 = a_w2;
        // m_h = a_h;
        // m_q1 = a_q1;
        // m_q2 = a_q2;
        // m_rho_s0 = a_rho_s0;
        // m_rho_g0 = a_rho_g0;
        // m_F0 = a_F0;
        // m_T = a_T;

        for (int i=0; i<12; i++)
        {
            rho_se[i] = m_rho_s0/24.0;
            rho_ss[i] = m_rho_s0/24.0;
            rho_gen[i] = m_rho_g0/36.0;
            rho_get[i] = m_rho_g0/36.0;
            rho_gsm[i] = m_rho_g0/36.0;
            gam[i] = 0.0; 
            gamdot[i] = 0.0;
        }
    }

    void ComputePlasticStrainRate(  Set::Matrix sig, 
                                    StrainGradientCrystalPlastic gradX,
                                    StrainGradientCrystalPlastic gradY,
                                    StrainGradientCrystalPlastic gradZ,
                                    Set::Scalar dt)
    {
        // gradX, gradY, and gradZ are the model gradient in the x, y, and z directions
        // This means that gradX.gamdotalpha[i] is the gradient of gamma_dot in the x direction for the i^{th} slip system.
        // Along the same lines gradY.gamalpha[i] is the gradient of gamma in the y direction for the i^{th} slip system.

        Eigen::Matrix3d sig_3D = Eigen::Matrix3d::Zero();
        sig_3D(Eigen::seq(0,AMREX_SPACEDIM-1),Eigen::seq(0,AMREX_SPACEDIM-1)) = sig;

        Eigen::Matrix3d R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        std::array<Eigen::Vector3d, 4> plane_rot;
        std::array<Eigen::Vector3d, 12> dir_rot;

        for (int i=0; i<4; i++)
        {
            plane_rot[i] = R_crystal*plane[i];
            for(int j=0; j<3; j++)
            {
                int alpha = 3*i + j;
                dir_rot[alpha] = R_crystal*dir[i][j];
            }
        }

        Set::Scalar total_rho = 0.0;
        for (int i = 0; i<12; i++) 
            total_rho += rho_se[i] + rho_ss[i] + std::abs(rho_gen[i]) + std::abs(rho_get[i]) + std::abs(rho_gsm[i]);
        
        // Let's first compute all dislocation densities
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                int alpha = 3*i + j;
                
                Set::Scalar anisotropic_rho = 0.0;
                for (int k=0; k<12; k++) anisotropic_rho += (m_w1 + (k == alpha ? 1.0 : 0.0)*(1.0-m_w2))*m_h*(rho_se[k] + rho_ss[k] + std::abs(rho_gen[k]) + std::abs(rho_get[k]) + std::abs(rho_gsm[k]));
                
                Eigen::Vector3d t_alpha = (dir_rot[alpha]).cross(plane_rot[i]);

                Set::Scalar rho_se_dot = ((m_Ke/m_b)*std::sqrt(anisotropic_rho) - m_de*rho_se[alpha]/m_b)*std::abs(gamdot[alpha]);

                Set::Scalar rho_ss_dot = (((2.0*m_Ke) - (Set::Constant::Pi*m_ds*m_ds*rho_ss[alpha]*0.25))*(1.0/m_b)*std::sqrt(anisotropic_rho) - m_ds*rho_ss[alpha]/m_b)*std::abs(gamdot[alpha]);

                Eigen::Vector3d curlTerm = Eigen::Vector3d::Zero();

                Eigen::Vector3d gamdotgradalpha(gradX.gamdot[alpha],gradY.gamdot[alpha],gradZ.gamdot[alpha]);
                for (int bi=0; bi<4; bi++)
                {
                    for (int bj=0; bj<3; bj++)
                    {
                        int beta = 3*bi + bj;
                        Eigen::Vector3d gamgradbeta(gradX.gam[beta],gradY.gam[beta],gradZ.gam[beta]);
                        
                        Eigen::Vector3d temp1 = (gamdotgradalpha*gam[beta] + gamdot[alpha]*gamgradbeta);
                        curlTerm += temp1.cross(dir_rot[beta])*((plane_rot[i]).dot(plane_rot[bi]));
                    }
                }
                
                Set::Scalar rho_gsm_dot = (1.0/m_b)*(dir_rot[alpha]).dot(curlTerm);
                Set::Scalar rho_get_dot = (1.0/m_b)*(t_alpha).dot(curlTerm);
                Set::Scalar rho_gen_dot = (1.0/m_b)*(plane_rot[i]).dot(curlTerm);

                rho_se[alpha] += rho_se_dot*dt;
                rho_ss[alpha] += rho_ss_dot*dt;
                rho_gsm[alpha] += rho_gsm_dot*dt;
                rho_get[alpha] += rho_get_dot*dt;
                rho_gen[alpha] += rho_gen_dot*dt;
            }
        }

        total_rho = 0.0;
        for (int k = 0; k<12; k++) total_rho += rho_se[k] + rho_ss[k] + std::abs(rho_gen[k]) + std::abs(rho_get[k]) + std::abs(rho_gsm[k]);
        
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                int alpha = 3*i+j;
                Set::Scalar tau = (dir_rot[alpha]).dot(sig_3D*(plane_rot[i]));
                Set::Scalar tau_prime = m_tau0 + m_mu*m_b*m_c0*std::sqrt(total_rho);

                Set::Scalar anisotropic_rho = 0.0;
                for (int k=0; k<12; k++) anisotropic_rho += (m_w1 + (k == alpha ? 1.0 : 0.0)*(1.0-m_w2))*m_h*(rho_se[k] + rho_ss[k] + std::abs(rho_gen[k]) + std::abs(rho_get[k]) + std::abs(rho_gsm[k]));
                Set::Scalar S_alpha = m_mu*m_b*std::sqrt(anisotropic_rho);

                Set::Scalar temp1 = (std::abs(tau) - S_alpha) / tau_prime;
                temp1 = temp1 > 0.0 ? std::pow(temp1,m_q1) : 0.0;
                Set::Scalar temp2 = 1. - temp1;
                temp2 = temp2 > 0.0 ? std::pow(temp2,m_q2) : 0.0;
                //if (outflag) Util::Message(INFO,"temp1 = ", temp1,". temp2 = ", temp2);
                gamdot[alpha] = m_gam0dot_sgcp[alpha]*std::exp(-m_F0*temp2/(m_T))*(tau>0.0 ? 1.0 : (tau < 0.0 ? -1.0 : 0.0));
            }
        }
        
    }

    void Advance(Set::Scalar dt, Set::Matrix /*eps*/, Set::Matrix /*sig*/) override
    {
        Eigen::Matrix3d R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        //std::array<Set::Scalar, 12> gamdot;
        
        for (int i=0; i< 4; i++)
        {
            for (int j=0; j<3; j++)
            {
                Eigen::Matrix3d dir_rot = ((R_crystal*dir[i][j]) * (R_crystal*plane[i]).transpose());
                F0 += dt * gamdot[3*i + j] * dir_rot(Eigen::seq(0,AMREX_SPACEDIM-1),Eigen::seq(0,AMREX_SPACEDIM-1));
                gam[3*i + j] += dt*gamdot[3*i + j];
            }
        }

    }

public:

    static const bool gradient_required = true;

    // For now FCC is hardcoded in here.
    inline static const std::array<Eigen::Vector3d, 4> plane = 
    {
        Eigen::Vector3d(1,1,1).normalized(),
        Eigen::Vector3d(-1,-1,1).normalized(),
        Eigen::Vector3d(-1,1,1).normalized(),
        Eigen::Vector3d(1,-1,1).normalized()
    };
    inline static const std::array<std::array<Eigen::Vector3d,3>,4> dir = 
    {
        Eigen::Vector3d(0,-1,1).normalized(), Eigen::Vector3d(1,0,-1).normalized(), Eigen::Vector3d(-1,1,0).normalized(),
        Eigen::Vector3d(0,1,1).normalized(), Eigen::Vector3d(-1,0,-1).normalized(), Eigen::Vector3d(1,-1,0).normalized(),
        Eigen::Vector3d(0,-1,1).normalized(), Eigen::Vector3d(-1,0,-1).normalized(), Eigen::Vector3d(1,1,0).normalized(),
        Eigen::Vector3d(0,1,1).normalized(), Eigen::Vector3d(1,0,-1).normalized(), Eigen::Vector3d(-1,-1,0).normalized()
    };

    std::array<Set::Scalar,12> gam;
    std::array<Set::Scalar,12> gamdot;
    std::array<Set::Scalar,12> rho_gsm;
    std::array<Set::Scalar,12> rho_get;
    std::array<Set::Scalar,12> rho_gen;
    std::array<Set::Scalar,12> rho_se;
    std::array<Set::Scalar,12> rho_ss;

    static std::array<Set::Scalar,12> m_gam0dot_sgcp;
    static Set::Scalar m_Ke;
    static Set::Scalar m_b;
    static Set::Scalar m_de;
    static Set::Scalar m_ds;
    static Set::Scalar m_tau0;
    static Set::Scalar m_mu;
    static Set::Scalar m_c0;
    static Set::Scalar m_rho_s0;
    static Set::Scalar m_rho_g0;
    static Set::Scalar m_w1, m_w2, m_h;
    static Set::Scalar m_q1, m_q2;
    static Set::Scalar m_F0, m_T;
    
public:
    static StrainGradientCrystalPlastic Zero()
    {
        StrainGradientCrystalPlastic ret;
        ret.q.coeffs() *= 0;
        ret.C11 = 0.0; ret.C12 = 0.0; ret.C44 = 0.0;
        ret.F0 = Set::Matrix::Zero();
        ret.DefineSlipGradientSystem();
        return ret;
    }

    static void Parse(StrainGradientCrystalPlastic & value, IO::ParmParse & pp)
    {
        Cubic::Parse(value,pp);

        value.m_Ke = 35.0e3;
        value.m_b = 0.248; //nm
        value.m_de = 5.e0;
        value.m_ds = 25.0e0;
        value.m_tau0 = 3.40e-1; //GPa
        value.m_mu = value.C44; //0.2*(3.0*value.C44 + value.C11 - value.C12);
        value.m_c0 = 8.7e-5;
        value.m_w1 = 1.0, value.m_w2 = 1.0, value.m_h = 0.084*0.084;
        value.m_q1 = 0.5, value.m_q2 = 1.25e2;
        value.m_rho_s0 = 6.e12;
        value.m_rho_g0 = 6.e12;
        value.m_F0 = 3.114e4; // divided by kB
        value.m_T = 298.;

        pp.query("Ke", value.m_Ke);
        pp.query("b", value.m_b);
        pp.query("de", value.m_de);
        pp.query("ds", value.m_ds);
        pp.query("tau0", value.m_tau0);
        pp.query("mu", value.m_mu);
        pp.query("c0", value.m_c0);
        pp.query("w1", value.m_w1);
        pp.query("w2", value.m_w2);
        pp.query("h", value.m_h);
        pp.query("q1", value.m_q1);
        pp.query("q2", value.m_q2);
        pp.query("rho_s0", value.m_rho_s0);
        pp.query("rho_g0", value.m_rho_g0);
        pp.query("m_F0", value.m_F0);
        pp.query("T", value.m_T);

        amrex::Vector<Set::Scalar> input_gam0dot;
        if(pp.contains("gam0dot"))
        {
            pp.queryarr("gam0dot",input_gam0dot);
            if(input_gam0dot.size() == 1)
            {
                for (int i = 0; i< 12; i++) 
                    value.m_gam0dot_sgcp[i] = input_gam0dot[0];
            }
            else if(input_gam0dot.size() <= 0 || input_gam0dot.size() != 12) Util::Abort(INFO, "Invalid number of values in gam0dot");
            else
            {
                for (int i = 0; i< 12; i++) 
                    value.m_gam0dot_sgcp[i] = input_gam0dot[i];
            }
        }
        else
            for (int i = 0; i< 12; i++) value.m_gam0dot_sgcp[i] = 0.001;

        value.DefineSlipGradientSystem();

        Util::Message(INFO, "modulus = ", value.C11, " ", value.C12, " ", value.C44);
        //Util::Message(INFO, "modulus = ", value.C110, " ", value.C120, " ", value.C440);
    }

    #define OP_CLASS StrainGradientCrystalPlastic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) 
    #define OP_LOOP_NUM 12
    #define OP_LOOP_VARS X(gam) X(gamdot) X(rho_se) X(rho_ss) X(rho_gen) X(rho_gsm) X(rho_get)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

Set::Scalar StrainGradientCrystalPlastic::m_Ke;
Set::Scalar StrainGradientCrystalPlastic::m_b;
Set::Scalar StrainGradientCrystalPlastic::m_de;
Set::Scalar StrainGradientCrystalPlastic::m_ds;
Set::Scalar StrainGradientCrystalPlastic::m_tau0;
Set::Scalar StrainGradientCrystalPlastic::m_mu;
Set::Scalar StrainGradientCrystalPlastic::m_c0;
Set::Scalar StrainGradientCrystalPlastic::m_rho_s0;
Set::Scalar StrainGradientCrystalPlastic::m_rho_g0;
Set::Scalar StrainGradientCrystalPlastic::m_w1, StrainGradientCrystalPlastic::m_w2, StrainGradientCrystalPlastic::m_h;
Set::Scalar StrainGradientCrystalPlastic::m_q1, StrainGradientCrystalPlastic::m_q2;
Set::Scalar StrainGradientCrystalPlastic::m_F0, StrainGradientCrystalPlastic::m_T;
std::array<Set::Scalar, 12> StrainGradientCrystalPlastic::m_gam0dot_sgcp;

}
}
}

template<>
inline int Set::Field<Model::Solid::Affine::StrainGradientCrystalPlastic>::NComp() const 
{
    return AMREX_SPACEDIM*AMREX_SPACEDIM;
}

template<>
inline std::string Set::Field<Model::Solid::Affine::StrainGradientCrystalPlastic>::Name(int i) const 
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
inline void Set::Field<Model::Solid::Affine::StrainGradientCrystalPlastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::StrainGradientCrystalPlastic> const & src = ((*this)[a_lev])->array(mfi);
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
