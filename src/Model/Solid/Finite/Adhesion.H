#ifndef MODEL_SOLID_FINITE_ADHESION_H_
#define MODEL_SOLID_FINITE_ADHESION_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include <cmath> // For std::pow and std::exp
#include </usr/include/eigen3/Eigen/Dense> // For Eigen::Matrix3d
#include <array> // For std::array to represent fourth-order tensors
#include "Set/Set.H"

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
                    result[i][j][k][l] = -FinvT(i, j) * FinvT(j,i);
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

public:
    Set::Scalar d = NAN, mu = NAN, kappa = NAN, zeta = NAN, n = NAN;

    // Static methods for Zero and Random initialization
    static Adhesion Zero() {
        Adhesion adhesion;
        adhesion.d = 0;
        adhesion.mu = 0;
        adhesion.kappa = 0;
        adhesion.zeta = 0;
        adhesion.n = 0;
        return adhesion;
    }

    static Adhesion Random() {
        Adhesion adhesion;
        adhesion.d = 1;//Util::Random();
        adhesion.mu = 0;//Util::Random();
        adhesion.kappa = 1;//Util::Random();
        adhesion.zeta =1; //Util::Random();
        adhesion.n =2; //Util::Random();
        return adhesion;
    }

    // Parsing method
    static void Parse(Adhesion& value, IO::ParmParse& pp) {
        pp.query("d", value.d);
        pp.query("mu", value.mu);
        pp.query("kappa", value.kappa);
        pp.query("zeta", value.zeta);
        pp.query("n", value.n);
    }

    // Potential Energy Function W(F)
    Set::Scalar W(const Set::Matrix& a_F) const override {
        Eigen::Matrix3d F = a_F;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Set::Scalar w = 0.0;

        w += d * (0.5 * mu * ((F * F.transpose()).trace() - 3.0) / J23 
                  + 0.5 * kappa * (J - 1.0) * (J - 1.0));
        w += zeta / std::pow(J, n);
        return w;
    }

    // First Derivative: DW(F)
    Set::Matrix DW(const Set::Matrix& a_F) const override {
        Eigen::Matrix3d F = a_F;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Set::Scalar J53 = std::pow(fabs(J), 5.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        Eigen::Matrix3d dw = d * mu * ( F / J23 
                          - (1.0 / 3.0) * ((F * F.transpose()).trace() - 3.0) * J * FinvT / J53)
                          + d * kappa * (J - 1.0) * J * FinvT
                          - zeta * n * std::pow(J, -n) * FinvT;

        return dw;
    }

    // Second Derivative: DDW(F)
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& a_F) const override {
        Eigen::Matrix3d F = a_F;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
        Set::Scalar J49 = std::pow(fabs(J), 4.0 / 9.0);
        Set::Scalar J53 = std::pow(fabs(J), 5.0 / 3.0);
        Set::Scalar J83 = std::pow(fabs(J), 8.0 / 3.0);
        Eigen::Matrix3d FinvT = F.inverse().transpose();

        auto Deriv = ::Set::Derivative(FinvT);

        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ddw;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    for (int l = 0; l < 3; ++l) {
                        // Shear-related term (mu)
                        Set::Scalar t0 = 0.0;
                        Set::Scalar ms = 0.0;
                        Set::Scalar t1 = 0.0;
                        Set::Scalar t11 = 0.0;

                        if (i == k && j == l) {
                            t0 = 1.0 / J23;  // Identity term contribution
                        }

                        t0 -= (2.0 / 3.0) * F(i, j) * FinvT(k, l) / J23;

                        // Compute the derivative of the second term in the formula
                        ms -= (1.0 / 3.0) * FinvT(k, j) * FinvT(i, l) * ((F * F.transpose()).trace() - 3) 
                            - (2.0 / 3.0) * F(i, j) * FinvT(i, j);

                        // Compute the final term that accounts for J^-2/3 scaling
                        t1 += (ms - (2.0 / 3.0) * FinvT(k, l) * FinvT(i, j) * ((F * F.transpose()).trace() - 3)) / (J23 * J23);


                        // **Fixed Parentheses Here**
                        t11 = t0 + t1;

                        // Bulk-related term (kappa)
                        Set::Scalar t2 = 0.0;
                        t2 += (2 * J * J * FinvT(k, l) - J * FinvT(k, l)) * FinvT(i, j) + (J - (J * J)) * (FinvT(k, j) * FinvT(i, l));

                        
                        

                        // Adhesion-related term (zeta)
                        Set::Scalar t3 = 0.0;
                        t3 += std::pow(n, 2) * std::pow(J, -n ) * FinvT(k,l) * FinvT(i,j) + n * std::pow(J, -n ) * FinvT(k, j ) * FinvT(i,l);

                        // Combine terms
                        ddw(i, j, k, l) = d * mu * t11 + d * kappa * t2 +zeta* t3;
                    }
                }
            }
        }

        return ddw;
    }

    virtual void Print(std::ostream& out) const override {
        out << "d = " << d << ", mu = " << mu 
            << ", kappa = " << kappa << ", zeta = " << zeta 
            << ", n = " << n;
    }
    #define OP_CLASS Adhesion
    #define OP_VARS X(kappa) X(mu) X(n) X(zeta) X(d)
    #include "Model/Solid/InClassOperators.H"

};
#include "Model/Solid/ExtClassOperators.H"

} // namespace Finite
} // namespace Solid
} // namespace Model

#endif
