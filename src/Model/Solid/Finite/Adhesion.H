#ifndef MODEL_SOLID_FINITE_ADHESION_H_
#define MODEL_SOLID_FINITE_ADHESION_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include <cmath> // For std::pow and std::exp
#include </usr/include/eigen3/Eigen/Dense> // For Eigen::Matrix3d
#include <array> // For std::array to represent fourth-order tensors

namespace Set {

// Outer product for two 3x3 matrices producing a 3x3x3x3 tensor
inline std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> 
Outer(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B) {
    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> result = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    result[i][j][k][l] = A(i, j) * B(k, l);
    return result;
}

// Derivative of inverse transpose: d(F^-T) / dF
inline std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> 
Derivative(const Eigen::Matrix3d& FinvT) {
    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> result = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    result[i][j][k][l] = -FinvT(i, j) * FinvT(k, l);
    return result;
}

} // namespace Set

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
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0); // Compressive term
        w += zeta / std::pow(J, n); // Adhesion term

        return w;
    }

    // First Derivative: DW(F)
    Set::Matrix DW(const Set::Matrix& a_F) const override {
        Eigen::Matrix3d F = a_F;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        Eigen::Matrix3d dw = d * mu * (F / J23 - (F * F.transpose()).trace() * FinvT / (3.0 * J23))
                            + kappa * (J - 1.0) * J * FinvT
                            - zeta * n / std::pow(J, n + 1) * FinvT;

        return dw;
    }

    // Second Derivative: DDW(F)
    auto DDW(const Set::Matrix& a_F) const override {
        Eigen::Matrix3d F = a_F;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        auto ddw = Set::Outer(F / J23, F / J23);
        auto FinvT_outer = Set::Outer(FinvT, FinvT);
        auto dFinvT = Set::Derivative(FinvT);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        ddw[i][j][k][l] += -((F * F.transpose()).trace() * FinvT_outer[i][j][k][l] / (3.0 * J23)
                                           + kappa * J * J * FinvT_outer[i][j][k][l]
                                           + zeta * n * (n * std::pow(J, -n - 2) * FinvT_outer[i][j][k][l]
                                           + std::pow(J, -n - 1) * dFinvT[i][j][k][l]));

        return ddw;
    }

    virtual void Print(std::ostream& out) const override {
        out << "d = " << d << ", mu = " << mu 
            << ", kappa = " << kappa << ", zeta = " << zeta 
            << ", n = " << n;
    }

public:
    Set::Scalar d = NAN, mu = NAN, kappa = NAN, zeta = NAN, n = NAN;
};

} // namespace Finite
} // namespace Solid
} // namespace Model

#endif
