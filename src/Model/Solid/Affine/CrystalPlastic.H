//
// Small-strain crystal plasticity with orthotropic elasticity and a thermal
// eigenstrain, for quasi-static solves through Integrator::Base::Mechanics.
//
// .. math::
//
//    \sigma = w\,\mathbb{C}:(\varepsilon - \varepsilon_{th} - \varepsilon_p),\qquad
//    \varepsilon_{th} = \alpha (T - T_0)\mathbf{I},\qquad
//    \dot\varepsilon_p = \sum_n \dot\gamma_n\,\mathrm{sym}(\mathbf{a}_n\otimes\mathbf{N}_n)
//
// The stiffness tensor, slip systems, and rate/hardening/thermal-softening
// parameters are taken from :ref:`Model::Solid::Finite::CrystalPlastic`, so the
// same inputs (C11..C66, slip_systems, tau_crss, ...) apply. The weight
// :math:`w` is set by the integrator (e.g. from the remaining solid fraction).
// In 2D the solve is plane strain: :math:`\varepsilon_{zz}=0`, and the full
// 3D stress is used to resolve shear on the 3D slip systems.
//

#ifndef MODEL_SOLID_AFFINE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_AFFINE_CRYSTALPLASTIC_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Finite/CrystalPlastic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class CrystalPlastic : public Solid<Set::Sym::MajorMinor>
{
public:
    static constexpr const char* name = "crystalplastic";
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    using Voigt = Eigen::Matrix<Set::Scalar,6,6>;
    using Slip = Eigen::Matrix<Set::Scalar,12,1>;
    using SlipMatrix = Eigen::Matrix<Set::Scalar,3,12>;

    CrystalPlastic() = default;
    ~CrystalPlastic() = default;

    // Build from a parsed finite CP model: rotated stiffness and lab-frame
    // slip systems are computed once here.
    void Define(const Finite::CrystalPlastic& cp)
    {
        Set::Quaternion q = cp.q;
        q.normalize();
        const Eigen::Matrix3d R = q.toRotationMatrix();
        const Set::Scalar Cm[6][6] = {
            {cp.C11, cp.C12, cp.C13, 0, 0, 0},
            {cp.C12, cp.C22, cp.C23, 0, 0, 0},
            {cp.C13, cp.C23, cp.C33, 0, 0, 0},
            {0, 0, 0, cp.C44, 0, 0},
            {0, 0, 0, 0, cp.C55, 0},
            {0, 0, 0, 0, 0, cp.C66}};
        auto material = [&](int a, int b, int c, int d) -> Set::Scalar
        { return Cm[V(a,b)][V(c,d)]; };
        for (int I = 0; I < 6; ++I)
            for (int J = 0; J < 6; ++J)
            {
                const int i = Vi(I), j = Vj(I), k = Vi(J), l = Vj(J);
                Set::Scalar c = 0.0;
                for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                for (int e = 0; e < 3; ++e) for (int f = 0; f < 3; ++f)
                    c += R(i,a)*R(j,b)*R(k,e)*R(l,f)*material(a,b,e,f);
                C(I,J) = c;
            }
        for (int n = 0; n < 12; ++n)
        {
            slip_directions.col(n) = R * cp.slip_directions.col(n);
            slip_normals.col(n) = R * cp.slip_normals.col(n);
        }
        tau_crss = cp.tau_crss;
        gammadot0 = cp.gammadot0;
        m_rate_inv = cp.m_rate_inv;
        hardening_self = cp.hardening_self;
        hardening_latent = cp.hardening_latent;
        temperature_reference = cp.temperature_reference;
        thermal_softening_coefficient = cp.thermal_softening_coefficient;
        minimum_crss_fraction = cp.minimum_crss_fraction;
        tstart = cp.tstart;
    }

    // Unweighted 3D stress for a given in-plane displacement gradient.
    Set::Matrix3d Stress3(const Set::Matrix& gradu) const
    {
        Set::Matrix3d eps = Set::Matrix3d::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                eps(i,j) = 0.5*(gradu(i,j) + gradu(j,i));
        eps -= eps_p;
        eps -= thermal_strain * Set::Matrix3d::Identity();
        Eigen::Matrix<Set::Scalar,6,1> e, s;
        for (int I = 0; I < 6; ++I)
            e(I) = (I < 3 ? 1.0 : 2.0) * eps(Vi(I),Vj(I));
        s = C * e;
        Set::Matrix3d sig;
        for (int I = 0; I < 6; ++I)
        {
            sig(Vi(I),Vj(I)) = s(I);
            sig(Vj(I),Vi(I)) = s(I);
        }
        return sig;
    }

    Set::Scalar W(const Set::Matrix& gradu) const
    {
        const Set::Matrix3d sig = Stress3(gradu);
        Set::Matrix3d eps = Set::Matrix3d::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                eps(i,j) = 0.5*(gradu(i,j) + gradu(j,i));
        eps -= eps_p + thermal_strain * Set::Matrix3d::Identity();
        return 0.5 * weight * (sig.cwiseProduct(eps)).sum();
    }
    Set::Matrix DW(const Set::Matrix& gradu) const
    {
        const Set::Matrix3d sig = Stress3(gradu);
        Set::Matrix ret;
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                ret(i,j) = weight * sig(i,j);
        return ret;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix& /*gradu*/) const
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> ret;
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
        for (int k = 0; k < AMREX_SPACEDIM; ++k)
        for (int l = 0; l < AMREX_SPACEDIM; ++l)
            ret(i,j,k,l) = weight * C(V(i,j),V(k,l));
        return ret;
    }

    Set::Scalar Resistance(int n) const
    {
        const Set::Scalar dT = Util::Max(temperature - temperature_reference, 0.0);
        const Set::Scalar thermal_factor = Util::Max(minimum_crss_fraction,
            std::exp(-thermal_softening_coefficient * dT));
        Set::Scalar r = tau_crss(n) * thermal_factor +
            hardening_self * gamma(n);
        for (int b = 0; b < 12; ++b)
            if (b != n) r += hardening_latent * gamma(b);
        return r;
    }

    // Called by Base::Mechanics after each solve with the solved
    // displacement gradient. Explicit update of slip; dt is subcycled so
    // that no system slips more than max_slip_increment per substep.
    void Advance(Set::Scalar dt, Set::Matrix gradu, Set::Matrix /*sigma*/, Set::Scalar time)
    {
        if (time < tstart || !(weight > active_weight) || !(dt > 0.0)) return;
        Set::Scalar remaining = dt;
        for (int iter = 0; iter < 1000 && remaining > 0.0; ++iter)
        {
            const Set::Matrix3d sig = Stress3(gradu);
            Slip rate;
            Set::Scalar max_rate = 0.0;
            for (int n = 0; n < 12; ++n)
            {
                const Set::Scalar tau = slip_directions.col(n).dot(sig * slip_normals.col(n));
                rate(n) = gammadot0 * std::pow(std::fabs(tau / Resistance(n)), m_rate_inv) *
                    (tau > 0 ? 1.0 : -1.0);
                max_rate = Util::Max(max_rate, std::fabs(rate(n)));
            }
            Set::Scalar h = remaining;
            if (max_rate * h > max_slip_increment) h = max_slip_increment / max_rate;
            for (int n = 0; n < 12; ++n)
            {
                const Set::Matrix3d an = slip_directions.col(n) * slip_normals.col(n).transpose();
                eps_p += rate(n) * h * 0.5 * (an + an.transpose());
                gamma(n) += std::fabs(rate(n)) * h;
            }
            remaining -= h;
        }
    }

    Set::Scalar TotalSlip() const { return gamma.sum(); }

    void Print(std::ostream& out) const { out << C; }

public:
    // Elastic and slip-system data (constant per node)
    Voigt C = Voigt::Zero();
    SlipMatrix slip_directions = SlipMatrix::Zero();
    SlipMatrix slip_normals = SlipMatrix::Zero();
    Slip tau_crss = Slip::Zero();
    Set::Scalar gammadot0 = 0.0, m_rate_inv = 1.0;
    Set::Scalar hardening_self = 0.0, hardening_latent = 0.0;
    Set::Scalar temperature_reference = 300.0;
    Set::Scalar thermal_softening_coefficient = 0.0, minimum_crss_fraction = 1.0;
    Set::Scalar tstart = 0.0;
    Set::Scalar max_slip_increment = 1.0e-4;
    Set::Scalar active_weight = 0.5;

    // Set by the integrator every step
    Set::Scalar weight = 1.0;
    Set::Scalar temperature = 300.0;
    Set::Scalar thermal_strain = 0.0;

    // History
    Set::Matrix3d eps_p = Set::Matrix3d::Zero();
    Slip gamma = Slip::Zero();

    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.weight = 0.0; ret.temperature = 0.0; ret.temperature_reference = 0.0;
        ret.m_rate_inv = 0.0; ret.minimum_crss_fraction = 0.0;
        ret.max_slip_increment = 0.0; ret.active_weight = 0.0;
        return ret;
    }

private:
    static int V(int i, int j)
    {
        if (i == j) return i;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
        return 5;
    }
    static int Vi(int I) { const int t[6] = {0,1,2,1,0,0}; return t[I]; }
    static int Vj(int I) { const int t[6] = {0,1,2,2,2,1}; return t[I]; }

public:
#define OP_CLASS CrystalPlastic
#define OP_VARS X(C) X(slip_directions) X(slip_normals) X(tau_crss) X(gammadot0) X(m_rate_inv) X(hardening_self) X(hardening_latent) X(temperature_reference) X(thermal_softening_coefficient) X(minimum_crss_fraction) X(tstart) X(max_slip_increment) X(active_weight) X(weight) X(temperature) X(thermal_strain) X(eps_p) X(gamma)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

template<>
inline int Set::Field<Model::Solid::Affine::CrystalPlastic>::NComp() const
{
    return 13;
}

template<>
inline std::string Set::Field<Model::Solid::Affine::CrystalPlastic>::Name(int i) const
{
    if (i < 12) return name + "_gamma" + std::to_string(i + 1);
    return name + "_total_slip";
}

template<>
inline void Set::Field<Model::Solid::Affine::CrystalPlastic>::Copy(int a_lev, amrex::MultiFab& a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (!bx.ok()) continue;
        amrex::Array4<const Model::Solid::Affine::CrystalPlastic> const& src = ((*this)[a_lev])->array(mfi);
        amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int n = 0; n < 12; ++n)
                dst(i,j,k,a_dstcomp + n) = src(i,j,k).gamma(n);
            dst(i,j,k,a_dstcomp + 12) = src(i,j,k).TotalSlip();
        });
    }
}

#endif
