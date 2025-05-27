//  CrystalPlastic.H — finite‑deformation FCC crystal plasticity (Newton implicit???)
//  ------------------------------------------------------------------------
//  Inherits PseudoLinear::Cubic for anisotropic elasticity
//  12 {111}<110> systems, rate‑dependent flow
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

namespace Model {
namespace Solid {
namespace Finite {

class CrystalPlastic : public PseudoLinear::Cubic
{
public:
    // ------------------------------------------------------------------
    //  Elastic energy, stress, tangent
    // ------------------------------------------------------------------
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
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> K = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for(int i=0;i<AMREX_SPACEDIM;i++)
          for(int J=0;J<AMREX_SPACEDIM;J++)
            for(int k=0;k<AMREX_SPACEDIM;k++)
              for(int L=0;L<AMREX_SPACEDIM;L++)
                for(int Q=0;Q<AMREX_SPACEDIM;Q++)
                  for(int R=0;R<AMREX_SPACEDIM;R++)
                    K(i,J,k,L) += DDW0(i,Q,k,R)*Fpinv(J,Q)*Fpinv(L,R);
        return K;
    }

    // ------------------------------------------------------------------
    //  Implicit Backward‑Euler update (local Newton)
    // ------------------------------------------------------------------
    void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix /*P_in*/, Set::Scalar time) override;

    // ---------------- state & parameters ----------------
    Set::Matrix3d                           Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,12,1>         gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Set::Scalar                             m_rate_inv = 0.05;  // 1/m
    Eigen::Matrix<Set::Scalar,12,1>         tau_crss  = Eigen::Matrix<Set::Scalar,12,1>::Constant(50.0);
    Set::Scalar                             gammadot0 = 1.0;    // s⁻¹
    Set::Scalar                             tstart    = 0.0;    // activate plasticity

    // ---------------- utility ----------------
    static CrystalPlastic Zero() { return CrystalPlastic(); }

private:
    static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12>
    slipSystems(const Set::Matrix3d &R);

    // ------------ member operator macros (needed by AMReX) ------------
public:
#define OP_CLASS CrystalPlastic
#define OP_VARS  X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)
#include "Model/Solid/InClassOperators.H"
#undef OP_VARS
#undef OP_CLASS
}; // class CrystalPlastic

// ------------------------------------------------------------------
//  Implementation of Advance 
// ------------------------------------------------------------------
void CrystalPlastic::Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix, Set::Scalar time)
{
    if (time < tstart) return;

    constexpr int NSYS = 12;
    const int MAX_IT = 15;
    const Set::Scalar RTOL = 1e-10;

    q.normalize();
    auto ss = slipSystems(q.toRotationMatrix());

    Eigen::Matrix<Set::Scalar,NSYS,1> dgamma = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Set::Matrix3d Fp_old = Fp, Fp_new = Fp;

    int converged = -1;
    for(int it=0; it<MAX_IT; ++it)
    {
        // Lp from current guess
        Set::Matrix3d Lp = Set::Matrix3d::Zero();
        for(int a=0;a<NSYS;++a)
            Lp += (dgamma[a]/dt) * ss[a].first * ss[a].second.transpose();

        Fp_new = (Lp*dt).exp() * Fp_old;
        Set::Matrix3d Fpinv = Fp_new.inverse();

        Set::Matrix3d Fe = Set::reduce(F) * Fpinv;
        Set::Matrix   Pe = PseudoLinear::Cubic::DW(Fe);
        Set::Matrix3d P3 = Set::expand(Pe);

        Eigen::Matrix<Set::Scalar,NSYS,1> R;
        Eigen::Matrix<Set::Scalar,NSYS,NSYS> J; J.setZero();

        for(int a=0;a<NSYS;++a)
        {
            Set::Scalar tau = ss[a].first.dot(P3 * ss[a].second);
            Set::Scalar abs_r = std::fabs(tau / tau_crss[a]);
            Set::Scalar sign  = tau >= 0.0 ? 1.0 : -1.0;
            Set::Scalar gdot  = gammadot0 * std::pow(abs_r, m_rate_inv) * sign;
            R[a] = dgamma[a] - dt * gdot;

            Set::Scalar dgdottau = (gammadot0 * m_rate_inv / tau_crss[a]) * std::pow(abs_r, m_rate_inv - 1.0) * sign;

            for(int b=0;b<NSYS;++b)
            {
                Set::Matrix3d dLp = (ss[b].first * ss[b].second.transpose()) / dt;
                Set::Matrix3d dFp = (dLp*dt) * Fp_new;
                Set::Matrix3d dFe = -Fe * dFp * Fpinv;
                Set::Matrix   dP  = PseudoLinear::Cubic::DW(Fe + dFe) - Pe;
                Set::Scalar   dtau = ss[a].first.dot(Set::expand(dP) * ss[a].second);
                J(a,b) = (a==b ? 1.0 : 0.0) - dt * dgdottau * dtau;
            }
        }

        if (R.norm() < RTOL) { converged = it; break; }

        Eigen::PartialPivLU<Eigen::Matrix<Set::Scalar,NSYS,NSYS>> lu(J);
        Eigen::Matrix<Set::Scalar,NSYS,1> delta = lu.solve(-R);
        dgamma += delta;
        if (delta.norm() < RTOL) { converged = it; break; }
    }

    Fp    = Fp_new;
    gamma += dgamma;

}

// ------------------------------------------------------------------
//  Slip‑system generator
// ------------------------------------------------------------------
std::array<std::pair<Set::Vector3d,Set::Vector3d>,12>
CrystalPlastic::slipSystems(const Set::Matrix3d &R)
{
    using V3 = Set::Vector3d;
    static const std::array<V3,12> n0 = {
        V3( 1, 1, 1),V3( 1, 1, 1),V3( 1, 1, 1),
        V3(-1, 1,-1),V3(-1, 1,-1),V3(-1, 1,-1),
        V3( 1,-1,-1),V3( 1,-1,-1),V3( 1,-1,-1),
        V3(-1,-1, 1),V3(-1,-1, 1),V3(-1,-1, 1)};
    static const std::array<V3,12> a0 = {
        V3( 1, 0,-1),V3( 0,-1, 1),V3( 1,-1, 0),
        V3( 1, 0,-1),V3( 1, 1, 0),V3( 0, 1, 1),
        V3( 1, 1, 0),V3( 0,-1, 1),V3( 1, 0, 1),
        V3( 0, 1, 1),V3( 1, 0, 1),V3( 1,-1, 0)};

    std::array<std::pair<V3,V3>,12> ss;
    for(int a=0;a<12;++a)
    {
        ss[a].first  = R * (a0[a]/std::sqrt(2.0));
        ss[a].second = R * (n0[a]/std::sqrt(3.0));
    }
    return ss;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------
#define OP_CLASS CrystalPlastic
#define OP_VARS  X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)
#include "Model/Solid/ExtClassOperators.H"
#undef OP_VARS
#undef OP_CLASS

} // namespace Finite bracket (will remove the comment later, using to make sure no brackets are missed)
} // namespace Solid
} // namespace Model

#endif // MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_

