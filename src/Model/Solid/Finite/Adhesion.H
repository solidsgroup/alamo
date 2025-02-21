#ifndef MODEL_SOLID_FINITE_ADHESION_H_
#define MODEL_SOLID_FINITE_ADHESION_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "Set/Set.H"
#include </usr/include/eigen3/Eigen/Dense> // For Eigen::Matrix3d
#include <array> // For std::array to represent fourth-order tensors
#include <cmath> // For std::pow and std::fabs

namespace Model {
namespace Solid {
namespace Finite {

class Adhesion : public Solid<Set::Sym::Major> {
  public:
    Adhesion() {}
    Adhesion(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {}
    virtual ~Adhesion() {}

    // User-specified material parameters:
    Set::Scalar d = NAN;    ///< Adhesion strength
    Set::Scalar E = NAN;    ///< Young's modulus
    Set::Scalar nu = NAN;   ///< Poisson's ratio
    Set::Scalar zeta = NAN; ///< Adhesion decay parameter
    Set::Scalar n = NAN;    ///< Decay exponent

    // Computed parameters (from E and nu):
    Set::Scalar mu = NAN;    ///< Shear modulus = E/(2(1+nu))
    Set::Scalar kappa = NAN; ///< Bulk modulus  = E/(3(1-2*nu))

    static Adhesion Zero() {
        Adhesion adhesion;
        adhesion.d = 0;
        adhesion.E = 0;
        adhesion.nu = 0;
        adhesion.zeta = 0;
        adhesion.n = 0;
        adhesion.mu = 0;
        adhesion.kappa = 0;
        return adhesion;
    }

    static Adhesion Random() {
        Adhesion adhesion;
        adhesion.d = Util::Random();
        adhesion.E = Util::Random();
        adhesion.nu = Util::Random();
        adhesion.zeta = Util::Random();
        adhesion.n = Util::Random();
        adhesion.mu = adhesion.E / (2.0 * (1.0 + adhesion.nu));
        adhesion.kappa = adhesion.E / (3.0 * (1.0 - 2.0 * adhesion.nu));
        return adhesion;
    }

    // Parse function: Only E, nu, d, zeta, and n are expected.
    // The shear and bulk moduli are computed from E and nu.
    static void Parse(Adhesion &value, IO::ParmParse &pp) {
        pp_query("d", value.d);
        pp_query("E", value.E);
        pp_query("nu", value.nu);
        pp_query("zeta", value.zeta);
        pp_query("n", value.n);
        // Compute moduli from E and nu:
        value.mu = value.E / (2.0 * (1.0 + value.nu));
        value.kappa = value.E / (3.0 * (1.0 - 2.0 * value.nu));
    }

    Set::Scalar W(const Set::Matrix &a_F) const override {
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
        Set::Scalar J23 = std::pow(std::fabs(J), 2.0 / 3.0);
        Set::Scalar w = 0.0;
        // Compute shear and bulk contributions (scaled by adhesion strength d)
        w += d * (0.5 * mu * (((F.transpose() * F).trace()) - 3) / J23 +
                  0.5 * kappa * (J - 1.0) * (J - 1.0));
        // Adhesion decay term:
        w += zeta / std::pow(J, n);
        return w;
    }

    Set::Matrix DW(const Set::Matrix &a_F) const override {
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
        Set::Scalar J23 = std::pow(std::fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();
        Eigen::Matrix3d dw =
            d * mu *
                (F / J23 - (1.0 / 3.0) * (((F.transpose() * F).trace()) - 3) *
                               FinvT / J23) +
            d * kappa * (J - 1.0) * J * FinvT -
            zeta * n * std::pow(J, -n) * FinvT;
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

    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>
    DDW(const Set::Matrix &a_F) const override {
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
        Set::Scalar J23 = std::pow(std::fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();
        Set::Scalar trace_FTF = (F.transpose() * F).trace();

        Set::Matrix4<3, Set::Sym::Major> ddw;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    for (int l = 0; l < 3; ++l) {
                        ddw(i, j, k, l) = 0.0;
                        Set::Scalar t1 = 0.0, t2 = 0.0, t3 = 0.0;
                        // Shear-related term (mu)
                        if (i == k && j == l)
                            t1 += 1.0;
                        t1 -= (2.0 / 3.0) * F(i, j) * FinvT(k, l);
                        t1 -= (2.0 / 3.0) * FinvT(i, j) * F(k, l);
                        t1 += (2.0 / 9.0) * (trace_FTF - 3) * FinvT(i, j) *
                              FinvT(k, l);
                        t1 += (1.0 / 3.0) * (trace_FTF - 3) * FinvT(i, l) *
                              FinvT(k, j);
                        // Bulk-related term (kappa)
                        t2 += (2.0 * J * J - J) * FinvT(i, j) * FinvT(k, l) +
                              (J - J * J) * FinvT(i, l) * FinvT(k, j);
                        // Adhesion-related term (zeta)
                        t3 += (n * n) * std::pow(J, -n) * FinvT(i, j) *
                                  FinvT(k, l) +
                              n * std::pow(J, -n) * FinvT(i, l) * FinvT(k, j);
                        ddw(i, j, k, l) =
                            d * mu * t1 / J23 + d * kappa * t2 + zeta * t3;
                    }
                }
            }
        }
#if AMREX_SPACEDIM == 2
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> r_ddw;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                    for (int l = 0; l < 2; ++l)
                        r_ddw(i, j, k, l) = ddw(i, j, k, l);
        return r_ddw;
#elif AMREX_SPACEDIM == 3
        return ddw;
#endif
    }

    virtual void Print(std::ostream &out) const override {
        out << "d = " << d << ", E = " << E << ", nu = " << nu
            << ", mu (computed) = " << mu << ", kappa (computed) = " << kappa
            << ", zeta = " << zeta << ", n = " << n;
    }

#define OP_CLASS Adhesion
#define OP_VARS X(d) X(E) X(nu) X(zeta) X(n)
#include "Model/Solid/InClassOperators.H"
};

#include "Model/Solid/ExtClassOperators.H"

} // namespace Finite
} // namespace Solid
} // namespace Model

#endif