#ifndef MODEL_SOLID_FINITE_ADHESION_H_
#define MODEL_SOLID_FINITE_ADHESION_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include <cmath>         // For std::pow and std::exp

namespace Model {
namespace Solid {
namespace Finite {

class Adhesion : public Solid<Set::Sym::Major> {
public:
    Adhesion() {};
    Adhesion(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~Adhesion() {};

    // Potential Energy Function W(F)
    Set::Scalar W(const Set::Matrix& a_F) const override {
#if AMREX_SPACEDIM == 2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM == 3
        Eigen::Matrix3d F = a_F;
#endif

        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Set::Scalar w = 0.0;

        // Energy terms:
        w += d * (0.5 * mu * ((F * F.transpose()).trace() / J23 - 3.0)); // Shear term
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);                        // Compressive term
        w += zeta * std::exp(-std::pow(J, n));                           // Adhesion/penetration term

        return w;
    }

    // First Derivative: DW(F)
    Set::Matrix DW(const Set::Matrix& a_F) const override {
#if AMREX_SPACEDIM == 2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM == 3
        Eigen::Matrix3d F = a_F;
#endif

        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        Eigen::Matrix3d dw = Eigen::Matrix3d::Zero();

        dw += d * mu * (F / J23 - (F * F.transpose()).trace() * FinvT / (3.0 * J23)); // Shear term
        dw += kappa * (J - 1.0) * J * FinvT;                                         // Compressive term
        dw += -zeta * n * std::pow(J, n - 1) * std::exp(-std::pow(J, n)) * FinvT;    // Adhesion term

#if AMREX_SPACEDIM == 2
        Set::Matrix r_dw;
        r_dw(0, 0) = dw(0, 0);
        r_dw(0, 1) = dw(0, 1);
        r_dw(1, 0) = dw(1, 0);
        r_dw(1, 1) = dw(1, 1);
        return r_dw;
#elif AMREX_SPACEDIM == 3
        return dw;
#endif
    }

    // Second Derivative: DDW(F)
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& a_F) const override {
#if AMREX_SPACEDIM == 2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM == 3
        Eigen::Matrix3d F = a_F;
#endif

        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ddw = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();

        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        // Tensor expressions for DDW: Expand based on second derivatives of each term
        // (Include similar logic from Neo-Hookean with added adhesion derivative contributions)

        return ddw;
    }

    virtual void Print(std::ostream& out) const override {
        out << "d = " << d << ", mu = " << mu 
            << ", kappa = " << kappa << ", zeta = " << zeta 
            << ", n = " << n;
    }

public:
    Set::Scalar d = NAN, mu = NAN, kappa = NAN, zeta = NAN, n = NAN;
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

    static Adhesion Zero() {
        Adhesion ret;
        ret.d = 0.0;
        ret.mu = 0.0;
        ret.kappa = 0.0;
        ret.zeta = 0.0;
        ret.n = 0.0;
        return ret;
    }

    static Adhesion Random() {
        Adhesion ret;
        ret.d = Util::Random();
        ret.mu = Util::Random();
        ret.kappa = Util::Random();
        ret.zeta = Util::Random();
        ret.n = Util::Random();
        return ret;
    }

    static void Parse(Adhesion& value, IO::ParmParse& pp) {
        pp_query("d", value.d);
        pp_query("mu", value.mu);
        pp_query("kappa", value.kappa);
        pp_query("zeta", value.zeta);
        pp_query("n", value.n);
    }

#define OP_CLASS Adhesion
#define OP_VARS X(d) X(mu) X(kappa) X(zeta) X(n)
#include "Model/Solid/InClassOperators.H"
};

#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
