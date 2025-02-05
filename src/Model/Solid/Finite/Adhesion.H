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
adhesion.d =Util::Random();
adhesion.mu = Util::Random();
adhesion.kappa = Util::Random();
adhesion.zeta =Util::Random();
adhesion.n =Util::Random();
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
#if AMREX_SPACEDIM==2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;
#endif
Set::Scalar J = F.determinant();
Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
Set::Scalar w = 0.0;
//All terms of w
w += d * (0.5 * mu * ((F.transpose()*F).trace() - 3) / J23 
+ 0.5 * kappa * (J - 1.0) * (J - 1.0));
w += zeta / std::pow(J, n);


return w;
}

// First Derivative: DW(F)
Set::Matrix DW(const Set::Matrix& a_F) const override {
#if AMREX_SPACEDIM==2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;
#endif
Set::Scalar J = F.determinant();
Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
Set::Scalar J53 = std::pow(fabs(J), 5.0 / 3.0);
Eigen::Matrix3d FinvT = F.inverse().transpose();

Eigen::Matrix3d dw = d * mu * ( F / J23  - (1.0 / 3.0) * ((F.transpose()*F).trace() - 3) * FinvT / J23) + d * kappa * (J - 1.0) * J * FinvT- zeta * n * std::pow(J, -n) * FinvT;


#if AMREX_SPACEDIM==2
        Set::Matrix r_dw;
        r_dw(0, 0) = dw(0, 0);
        r_dw(0, 1) = dw(0, 1);
        r_dw(1, 0) = dw(1, 0);
        r_dw(1, 1) = dw(1, 1);
        return r_dw;
#elif AMREX_SPACEDIM==3
        return dw;
#endif
}

// Second Derivative: DDW(F)
Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& a_F) const override {
#if AMREX_SPACEDIM==2
        Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
        F(0, 0) = a_F(0, 0);
        F(0, 1) = a_F(0, 1);
        F(1, 0) = a_F(1, 0);
        F(1, 1) = a_F(1, 1);
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;
#endif
Set::Scalar J = F.determinant();
Set::Scalar J23 = std::pow(fabs(J), 2.0 / 3.0);
Eigen::Matrix3d FinvT = F.inverse().transpose();
Set::Scalar trace_FT_F = (F.transpose() * F).trace(); // Compute tr(F^TF )

auto Deriv = ::Set::Derivative(FinvT);

Set::Matrix4<3, Set::Sym::Major> ddw;

for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++j) {
for (int k = 0; k < 3; ++k) {
for (int l = 0; l < 3; ++l) {
ddw(i, j, k, l) = 0.0;

Set::Scalar t1 = 0.0, t2 = 0.0, t3=0;
// Shear-related term (mu)
if (i == k && j == l) t1 += 1.0;
t1 -= (2. / 3.) * F(i, j) * FinvT(k, l);
t1 -= (2. / 3.) * FinvT(i, j) * F(k, l);
t1 += (2. / 9.) * ((F.transpose()*F).trace()-3) * FinvT(i, j) * FinvT(k, l);
t1 += (1. / 3.) * ((F.transpose()*F).trace()-3) * FinvT(i, l) * FinvT(k, j);

// Bulk-related term (kappa)
t2 += (2 * J * J * FinvT(k, l) - J * FinvT(k, l)) * FinvT(i, j) + (J - (J * J)) * (FinvT(k, j) * FinvT(i, l));


// Adhesion-related term (zeta)
t3 += std::pow(n, 2) * std::pow(J, -n ) * FinvT(k,l) * FinvT(i,j) + n * std::pow(J, -n ) * FinvT(k, j ) * FinvT(i,l);

// Combine terms
ddw(i, j, k, l) = d * mu * t1/J23 + d * kappa * t2 +zeta* t3;
}
}
}
}
#if AMREX_SPACEDIM==2
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> r_ddw;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++)
                        r_ddw(i, j, k, l) = ddw(i, j, k, l);
        return r_ddw;
#elif AMREX_SPACEDIM==3
        return ddw;
#endif
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

