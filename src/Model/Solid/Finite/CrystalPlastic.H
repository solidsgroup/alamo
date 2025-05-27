//
//  Implicit finite‑deformation FCC crystal‑plasticity model (compile‑ready)
//  ------------------------------------------------------------------------
//  * Cubic elastic base (inherits PseudoLinear::Cubic)
//  * 12 FCC slip systems, rate‑dependent power‑law flow rule
//  * Backward‑Euler + exponential map for Fp
//  * No hardening yet (tau_crss constant)
//
#ifndef MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_

#include <array>
#include <vector>
#include <cmath>
// Eigen: adjust path if Eigen installed elsewhere
#include </usr/include/eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>  // for debug printouts

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinear/Cubic.H"

namespace Model {
namespace Solid {
namespace Finite {

class CrystalPlastic : public PseudoLinear::Cubic
{
public:
    // -----------------------------------------------------------------
    //  Elastic energy, stress, tangent (multiplicative split)
    // -----------------------------------------------------------------
    Set::Scalar W(const Set::Matrix &F) const override
    { return PseudoLinear::Cubic::W(F * Set::reduce(Fp).inverse()); }

    Set::Matrix DW(const Set::Matrix &F) const override
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        return PseudoLinear::Cubic::DW(F * Fpinv) * Fpinv.transpose();
    }

    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix &F) const override
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        auto DDW0 = PseudoLinear::Cubic::DDW(F * Fpinv);
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for (int i=0;i<AMREX_SPACEDIM;++i)
            for (int J=0;J<AMREX_SPACEDIM;++J)
                for (int k=0;k<AMREX_SPACEDIM;++k)
                    for (int L=0;L<AMREX_SPACEDIM;++L)
                        for (int Q=0;Q<AMREX_SPACEDIM;++Q)
                            for (int R=0;R<AMREX_SPACEDIM;++R)
                                ret(i,J,k,L)+=DDW0(i,Q,k,R)*Fpinv(J,Q)*Fpinv(L,R);
        return ret;
    }

    // -----------------------------------------------------------------
    //  Implicit Backward‑Euler update (per material point)
    // -----------------------------------------------------------------
    void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix /*P_in*/, Set::Scalar time) override
    {
        #ifndef AMREX_DEVICE_COMPILE
        amrex::Print() << ">>> entering Advance, level time = " << time << "\n";
        #endif

        if (time < tstart) return;                     // elastic until tstart

        constexpr int NSYS = 12;
        const int MAX_IT = 25;
        const Set::Scalar TOL = 1e-10;

        q.normalize();
        auto ss = slipSystems(q.toRotationMatrix());   // slip dirs/normals

        Eigen::Matrix<Set::Scalar,NSYS,1> dgamma = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        Set::Matrix3d Fp_old = Fp, Fp_new;

        int converged_it = -1;
        for (int it=0; it<MAX_IT; ++it)
        {
            Set::Matrix3d Lp = Set::Matrix3d::Zero();
            for (int a=0;a<NSYS;++a)
                Lp += (dgamma[a]/dt)*ss[a].first*ss[a].second.transpose();

            Fp_new = (Lp*dt).exp()*Fp_old;
            Set::Matrix3d Fpinv = Fp_new.inverse();

            Set::Matrix3d Fe = Set::reduce(F)*Fpinv;
            Set::Matrix   Pe = PseudoLinear::Cubic::DW(Fe);
            Set::Matrix3d P3 = Set::expand(Pe);

            Eigen::Matrix<Set::Scalar,NSYS,1> dgamma_next;
            for (int a=0;a<NSYS;++a)
            {
                Set::Scalar tau = ss[a].first.dot(P3*ss[a].second);
                Set::Scalar gdot = gammadot0*std::pow(std::fabs(tau/tau_crss[a]), m_rate_inv)*(tau>=0?1.0:-1.0);
                dgamma_next[a] = dt*gdot;
            }
            if ((dgamma_next-dgamma).norm()<TOL){ dgamma=dgamma_next; converged_it=it; break; }
            dgamma = dgamma_next;
        }

        Fp = Fp_new;
        gamma += dgamma;

        // Debug message
        std::cout << "CrystalPlastic implicit solver time=" << time
                  << " dt=" << dt
                  << (converged_it>=0 ? " converged in " : " did NOT converge in ")
                  << (converged_it>=0 ? converged_it+1 : MAX_IT)
                  << " iterations." << std::endl;
    }

    // ---------------- state & parameters ----------------
    Set::Matrix3d Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,12,1> gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Set::Scalar m_rate_inv = 0.05;        // 1/m
    Eigen::Matrix<Set::Scalar,12,1> tau_crss = Eigen::Matrix<Set::Scalar,12,1>::Constant(50.0);
    Set::Scalar gammadot0 = 1.0;
    Set::Scalar tstart = 0.0;

private:
    // ------------- slip system generator -------------
    static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12>
    slipSystems(const Set::Matrix3d &R)
    {
        using V3 = Set::Vector3d;
        static const std::array<V3,12> n0 = {V3(1,1,1),V3(1,1,1),V3(1,1,1), V3(-1,1,-1),V3(-1,1,-1),V3(-1,1,-1), V3(1,-1,-1),V3(1,-1,-1),V3(1,-1,-1), V3(-1,-1,1),V3(-1,-1,1),V3(-1,-1,1)};
        static const std::array<V3,12> a0 = {V3(1,0,-1),V3(0,-1,1),V3(1,-1,0), V3(1,0,-1),V3(1,1,0),V3(0,1,1), V3(1,1,0),V3(0,-1,1),V3(1,0,1), V3(0,1,1),V3(1,0,1),V3(1,-1,0)};
        std::array<std::pair<V3,V3>,12> out;
        for(int a=0;a<12;++a)
        {
            out[a].first  = R*(a0[a]/std::sqrt(2.0));
            out[a].second = R*(n0[a]/std::sqrt(3.0));
        }
        return out;
    }

public:
    // --------------------- utility factories --------------------
    static CrystalPlastic Zero(){return CrystalPlastic();}

// -- member arithmetic operators --
#define OP_CLASS CrystalPlastic
#define OP_VARS  X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)
#include "Model/Solid/InClassOperators.H"

};

// ---------------- arithmetic/operator macros (free functions) ----------------
#define OP_CLASS CrystalPlastic
#define OP_VARS  X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)

#include "Model/Solid/ExtClassOperators.H"

} // namespace Finite
} // namespace Solid
} // namespace Model

#endif // MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_

