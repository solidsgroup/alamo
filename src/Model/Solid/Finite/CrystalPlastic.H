//  CrystalPlastic.H :finite-deformation FCC crystal plasticity (Newton implicit)
//  ------------------------------------------------------------------------
//  Inherits PseudoLinear::Cubic for anisotropic elasticity
//  12 {111}<110> systems, rate-dependent flow
//
#ifndef MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_

#include <array>
#include <vector>
#include <cmath>
#include </usr/include/eigen3/Eigen/Core>
#include </usr/include/eigen3/Eigen/LU>
#include </usr/include/eigen3/unsupported/Eigen/MatrixFunctions>
#ifndef AMREX_DEVICE_COMPILE
  #include <iostream>
#endif

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinear/Cubic.H"

namespace Model 
{
namespace Solid 
{
namespace Finite 
{

class CrystalPlastic : public PseudoLinear::Cubic
{
public:
    CrystalPlastic() {};

    virtual ~CrystalPlastic() {};
    // ------------------------------------------------------------------
    //  Elastic energy, stress, tangent
    // ------------------------------------------------------------------
    Set::Scalar W(const Set::Matrix &F) const override
    { 
        return PseudoLinear::Cubic::W(F * Set::reduce(Fp).inverse()); }

    Set::Matrix DW(const Set::Matrix &F) const override
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        return PseudoLinear::Cubic::DW(F * Fpinv) * Fpinv.transpose();
    }

    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix &F) const override
    {
        // compute inverse of plastic deformation
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();

        // pull back the 4th-order elasticity tensor
        auto DDW = PseudoLinear::Cubic::DDW(F * Fpinv);

        // initialize K to zero
        auto ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();

        // perform the six-fold contraction:
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            for (int J = 0; J < AMREX_SPACEDIM; J++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    for (int L = 0; L < AMREX_SPACEDIM; L++)
                    {
                        ret(i,J,k,L) = 0.0;
                        for (int Q = 0; Q < AMREX_SPACEDIM; Q++)
                            for (int R = 0; R < AMREX_SPACEDIM; R++)
                                ret(i,J,k,L) += DDW(i,Q,k,R)*Fpinv(J,Q)*Fpinv(L,R);
                    }
        return ret;
    }

virtual void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix P, Set::Scalar time) override
{
    if (time < tstart) return;

    const int MAX_IT = 10;
    const Set::Scalar RTOL = 1e-6;

    
    q.normalize();
    Set::Matrix3d Rmat = q.toRotationMatrix();
    std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ss = slipSystems(Rmat);

    
    #ifndef AMREX_DEVICE_COMPILE
    std::cout << "\n[DEBUG] Advance() start  time=" << time << "\n";
    std::cout << "  Input P dims: " << P.rows() << " x " << P.cols() << "\n";
    #endif

    Eigen::Matrix<Set::Scalar,12,1> dgamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Set::Matrix3d Fp_old = Fp, Fp_new = Fp;

    for(int it = 0; it < MAX_IT; ++it)
    {
        std::cout << "[DEBUG] === Newton iter " << it << " ===\n";

        //Set::Matrix3d Lp = Set::Matrix3d::Zero();
        Set::Matrix3d Lp = Set::Matrix3d::Zero();
        for(int a = 0; a < 12; ++a) {
            Lp += (dgamma[a] / dt) * ss[a].first * ss[a].second.transpose();
        }
        std::cout << "  Lp:\n" << Lp << "\n";

        Fp_new = (Lp * dt).exp() * Fp_old;
        Set::Matrix3d Fpinv = Fp_new.inverse();
        std::cout << "  Fp_old:\n" << Fp_old << "\n";
        std::cout << "  Fp_new:\n" << Fp_new << "\n";
        std::cout << "  Fpinv:\n" << Fpinv << "\n";
        std::cout << "  F:\n" << F << "\n";
        
        Set::Matrix3d Fe = Set::Matrix3d::Zero();
        Fe = Set::expand(F)* Fpinv;
        std::cout << "  Fe:\n" << Fe << "\n";

        Set::Matrix3d P3 = Set::expand(P);
        std::cout << "  P3:\n" << P3 << "\n";
        std::cout << "  P:\n" << P << "\n";

        Eigen::Matrix<Set::Scalar,12,1> Rvec;
        for(int a = 0; a < 12; ++a)
        {
            Set::Scalar tau   = ss[a].first.dot(P3 * ss[a].second);
            Set::Scalar abs_r = std::fabs(tau / tau_crss[a]);
            Set::Scalar sign  = (tau >= 0.0 ? 1.0 : -1.0);
            Set::Scalar gdot  = gammadot0 * std::pow(abs_r, m_rate_inv) * sign;
            Rvec[a] = dgamma[a] - dt * gdot;
        }

        std::cout << "  Rvec (size 12): [";
        for(int a = 0; a < 12; ++a) {
            std::cout << Rvec[a];
            if (a < 11) std::cout << ", ";
        }
        std::cout << "]   ||R|| = " << Rvec.norm() << "\n";

        if (Rvec.norm() < RTOL) {
            std::cout << "  Converged on iter " << it 
                      << " because ||R|| < " << RTOL << "\n";
            break;
        }

        Eigen::Matrix<Set::Scalar,12,12> J;
        J.setZero();

        for(int a = 0; a < 12; ++a)
        {
            Set::Scalar tau   = ss[a].first.dot(P3 * ss[a].second);
            Set::Scalar abs_r = std::fabs(tau / tau_crss[a]);

            const Set::Scalar EPS = 1e-8;
            if (abs_r < EPS) abs_r = EPS;

            Set::Scalar dgdottau = (gammadot0 * m_rate_inv / tau_crss[a])
                                 * std::pow(abs_r, m_rate_inv - 1.0);

            for(int b = 0; b < 12; ++b)
            {
                Set::Matrix3d dLp = (ss[b].first * ss[b].second.transpose()) / dt;

                Set::Matrix3d dFp = (dLp * dt) * Fp_new;

                Set::Matrix3d dFe = -Fe * dFp * Fpinv;

                Set::Matrix dP = PseudoLinear::Cubic::DW(Set::reduce(Fe + dFe)) - P;

                Set::Scalar dtau = ss[a].first.dot(Set::expand(dP) * ss[a].second);

                J(a,b) = (a == b ? Set::Scalar(1.0) : Set::Scalar(0.0))
                        - dt * dgdottau * dtau;
            }
        }

        std::cout << "  J (12x12):\n";
        for(int a = 0; a < 12; ++a)
        {
            std::cout << "    [ ";
            for(int b = 0; b < 12; ++b)
            {
                std::cout << J(a,b);
                if (b < 11) std::cout << ", ";
            }
            std::cout << " ]\n";
        }

        Eigen::PartialPivLU<Eigen::Matrix<Set::Scalar,12,12>> lu(J);
        Eigen::Matrix<Set::Scalar,12,1> delta = lu.solve(-Rvec);

        std::cout << "  delta (size 12): [";
        for(int a = 0; a < 12; ++a) {
            std::cout << delta[a];
            if (a < 11) std::cout << ", ";
        }
        std::cout << "]   ||delta|| = " << delta.norm() << "\n";

        dgamma += delta;
        if (delta.norm() < RTOL) {
            std::cout << "  Converged on iter " << it 
                      << " because ||delta|| < " << RTOL << "\n";
            break;
        }
    } 

    Fp    = Fp_new;
    gamma += dgamma;
}

public:
    //Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;
    //static constexpr KinematicVariable kinvar = KinematicVariable::F;
    //Set::Scalar C11=1.68, C12=1.21, C44=0.75;
    //mutable Set::Quaternion q = Set::Quaternion(1.0, 0.0, 0.0, 0.0);

    // Plastic state variables
    Set::Matrix3d Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,12,1> gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();

    // Parameters
    Set::Scalar m_rate_inv = NAN;
    //Set::Scalar tau_crss = NAN;
    Eigen::Matrix<Set::Scalar,12,1> tau_crss = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Set::Scalar gammadot0 = NAN;

    // Control variables
    Set::Scalar tstart = NAN;

private:
    static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> slipSystems(Set::Matrix3d R)
    {
        const std::array<Set::Vector3d,12> n0 =
            {
                Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
                Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
                Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
                Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
                Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
                Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
                Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
                Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
                Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
                Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0),
                Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0),
                Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0)
            };
        const std::array<Set::Vector3d,12> a0 =
            {
                Set::Vector3d::Constant( 1, 0,-1) / sqrt(2.0), //  1  1  1
                Set::Vector3d::Constant( 0,-1, 1) / sqrt(2.0), //  1  1  1
                Set::Vector3d::Constant( 1,-1, 0) / sqrt(2.0), //  1  1  1
                Set::Vector3d::Constant( 1, 0,-1) / sqrt(2.0), // -1  1 -1
                Set::Vector3d::Constant( 1, 1, 0) / sqrt(2.0), // -1  1 -1
                Set::Vector3d::Constant( 0, 1, 1) / sqrt(2.0), // -1  1 -1
                Set::Vector3d::Constant( 1, 1, 0) / sqrt(2.0), // 1  -1 -1
                Set::Vector3d::Constant( 0,-1, 1) / sqrt(2.0), // 1  -1 -1
                Set::Vector3d::Constant( 1, 0, 1) / sqrt(2.0), // 1  -1 -1
                Set::Vector3d::Constant( 0, 1, 1) / sqrt(2.0), // -1  -1 1
                Set::Vector3d::Constant( 1, 0, 1) / sqrt(2.0), // -1  -1 1
                Set::Vector3d::Constant( 1,-1, 0) / sqrt(2.0)  // -1  -1 1
            };


        static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ret;
        for (int n = 0 ; n < 12 ; n ++)
        {
            ret[n].first = R*a0[n];
            ret[n].second = R*n0[n];
        }
        return ret;
    }

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.C11 = 0.0;
        ret.C12 = 0.0;
        ret.C44 = 0.0;
        ret.q = Set::Quaternion(0.0,0.0,0.0,0.0);
        ret.Fp = Set::Matrix3d::Zero();
        ret.gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
        ret.m_rate_inv = 0.0;
        ret.tau_crss = Eigen::Matrix<Set::Scalar,12,1>::Zero();
        ret.gammadot0 = 0.0;
        return ret;
    }
    static CrystalPlastic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static CrystalPlastic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        CrystalPlastic ret = CrystalPlastic::Zero();
        ret.C11 = C11;
        ret.C12 = C12;
        ret.C44 = C44;
        ret.q = Set::Quaternion::UnitRandom();
        return ret;
    }
    static void Parse(CrystalPlastic & value, IO::ParmParse & pp)
    {
        PseudoLinear::Cubic::Parse(value,pp);
        // TODO Add inputs for the CP values

        std::vector<Set::Scalar> tau_crss;
        // Critical resolved shear stress :math:`\tau_{crss}`
        pp_queryarr("tau_crss",tau_crss);
        if (tau_crss.size() == 1)
            for (int n = 0; n < 12; n++) value.tau_crss[n] = tau_crss[0];
        else if (tau_crss.size() == 12)
            for (int n = 0; n < 12; n++) value.tau_crss[n] = tau_crss[n];
        else
            Util::Exception(INFO,"Expected either 1 or 12 values for tau_crss but got ",tau_crss.size());

        // Rate hardening coefficient :math:`\dot{\gamma}_0`
        pp_query_default("gammadot0",value.gammadot0,1.0);

        // Inverse of the hardening exponent :math:`\frac{1}{m}`
        pp_query_default("m_rate_inv",value.m_rate_inv,0.5);
        
        // Time to activate plastic slip
        pp_query_default("tstart",value.tstart,0.0);    
    }

#define OP_CLASS CrystalPlastic
#define OP_VARS X(C11) X(C12) X(C44) X(q) X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"





}
}
}


template<>
inline int Set::Field<Model::Solid::Finite::CrystalPlastic>::NComp() const 
{
    return 12;
}

template<>
inline std::string Set::Field<Model::Solid::Finite::CrystalPlastic>::Name(int i) const 
{
    if (i==0) return name + "gamma1";
    if (i==1) return name + "gamma2";
    if (i==2) return name + "gamma3";
    if (i==3) return name + "gamma4";
    if (i==4) return name + "gamma5";
    if (i==5) return name + "gamma6";
    if (i==6) return name + "gamma7";
    if (i==7) return name + "gamma8";
    if (i==8) return name + "gamma9";
    if (i==9) return name + "gamma10";
    if (i==10) return name + "gamma11";
    if (i==11) return name + "gamma12";
    return name;
}

template<>
inline void Set::Field<Model::Solid::Finite::CrystalPlastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Finite::CrystalPlastic> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < 12; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k).gamma(n);
                });
            }
        }
    }    
}

#endif 
