// This models a finite-strain, solid with adhesion behavior governed by
// a combination of elastic energy and an adhesion decay term.
//
// The energy and derivatives are:
//
// .. math::
//    :nowrap:
//
//    \begin{gather}
//    W = d \left[ \frac{\mu}{2J^{2/3}}(F^T F : I - 3) + \frac{\kappa}{2}(J - 1)^2 \right] + \frac{\zeta}{J^n} \\
//    DW = \text{(first Piola-Kirchhoff stress)} \\
//    DDW = \text{(tangent stiffness tensor)}
//    \end{gather}
//
// where:
// - :math:`d` The coefficient of degradable part,
// - :math:`\mu` is the shear modulus,
// - :math:`\kappa` is the bulk modulus,
// - :math:`\zeta` is the adhesion decay coefficient or the coefficient of non-degradable part,
// - :math:`n` is the decay exponent,
// - :math:`F` is the deformation gradient,
// - :math:`J = \det(F)` is the volume ratio,
// - :math:`F^T F : I = \text{trace}(F^T F)` is used to capture isochoric deformation.
//
// Notes:
//
// Adhesion energy model W(F):
// Consists of two parts:
//   1. Degradable (elastic shear + volumetric terms)
//   2. Non-degradable (adhesion term decaying with volume expansion)
// This model is rate-independent, reversible, and purely elastic.

#ifndef MODEL_SOLID_FINITE_ADHESION_H_
#define MODEL_SOLID_FINITE_ADHESION_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Model/Solid/Solid.H"
#include "Set/Set.H"
#include <array>              // For std::array to represent fourth-order tensors
#include <cmath>              // For std::pow and std::abs
#include <eigen3/Eigen/Dense> // For Eigen::Matrix3d

namespace Model
{
namespace Solid
{
namespace Finite
{

class Adhesion : public NeoHookean
{
public:
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

    Adhesion() {}
    Adhesion(Solid<Set::Sym::Major> base) : NeoHookean(base) {}
    virtual ~Adhesion() {}

    // Material parameters (declared in order)
    Set::Scalar d = NAN;    ///< Adhesion strength (degradable part coefficient)
    Set::Scalar zeta = NAN; ///< Adhesion decay parameter (non-degradable part)
    Set::Scalar n = NAN;    ///< Decay exponent

    // Return a zero-parameter adhesion model.
    static Adhesion
    Zero()
    {
        Adhesion adhesion;
        adhesion.d = 0.0;    ///< Adhesion strength (degradable part coefficient)
        adhesion.zeta = 0.0; ///< Adhesion decay parameter (non-degradable part)
        adhesion.n = 0.0;    ///< Decay exponent
        adhesion.mu = 0.0;
        adhesion.kappa = 0.0;
        return adhesion;
    }

    // Return a random-parameter adhesion model.
    static Adhesion
    Random()
    {
        Adhesion adhesion;
        adhesion.d = Util::Random() ; // Util::Random();    ///< Adhesion strength (degradable part coefficient)
        adhesion.zeta = Util::Random(); ///< Adhesion decay parameter (non-degradable part)
        adhesion.n = Util::Random();    ///< Decay exponent
        adhesion.mu = Util::Random();
        adhesion.kappa = Util::Random();// Util::Random();
        return adhesion;
    }

    // Parse function: Reads material parameters from the input file.
    // Here, we assume the input file provides: mu, kappa, d, zeta, and n.
    static void
    Parse(Adhesion &value, IO::ParmParse &pp)
    {
        pp.queryclass<NeoHookean>(value);
        // Initial value for damage
        pp_query_default("d", value.d, 1.0);
        // Adhesion decay parameter (non-degradable part)
        pp_query_default("zeta", value.zeta, 0.0);
        // Decay exponent
        pp_query_default("n", value.n, 1.0);
    }

    // Strain energy function:
    // W = d * [ 0.5*mu*(trace(F^T F)-3)/J^(2/3) + 0.5*kappa*(J-1)^2 ] + zeta/J^n
    Set::Scalar
    W(const Set::Matrix &F) const override
    {
        Set::Scalar w_nh = (1.0 - d) * NeoHookean::W(F);
        Set::Scalar J = F.determinant();
        Set::Scalar w_compress = zeta * std::pow(J, -n);
        return w_nh + w_compress;
    }

    // First Piola-Kirchhoff stress (derivative of W)
    Set::Matrix
    DW(const Set::Matrix &F) const override
    {
        Set::Matrix dw_nh = (1.0 - d) * NeoHookean::DW(F);
        Set::Scalar J = F.determinant();
        Set::Matrix FinvT = F.inverse().transpose();
        Set::Matrix dw_compress = zeta * (-n * std::pow(J, -n) * FinvT);
        return dw_nh + dw_compress;
    }

    // Tangent stiffness tensor (second derivative of W)
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>
    DDW(const Set::Matrix &F) const override
    {
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ddw_nh =
            (1.0 - d) * NeoHookean::DDW(F);

        Set::Scalar J = F.determinant();
        Set::Matrix FinvT = F.inverse().transpose();

        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ddw_compress;

        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
            {
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                {
                    for (int l = 0; l < AMREX_SPACEDIM; l++)
                    {
                        // Adhesion-related term (zeta)

                        // dFinvT(i,j) / dF(k,l)
                        //   = dFinv(j,i) / dF(k,l)
                        //   = -dFinvT(k,j) * dFinvT(iL)

                        Set::Scalar t3 = n * n * std::pow(J, -n) * FinvT(i, j) * FinvT(k, l)
                                         + n * std::pow(J, -n) * FinvT(k, j) * FinvT(i, l);
                        ddw_compress(i, j, k, l) = zeta * t3;
                    }
                }
            }
        }
        return ddw_nh + ddw_compress;
    }

    // Print material parameters
    virtual void
    Print(std::ostream &out) const override
    {
        out << "d = " << d << ", mu = " << mu << ", kappa = " << kappa
            << ", zeta = " << zeta << ", n = " << n;
    }

#define OP_CLASS Adhesion
#define OP_VARS X(d) X(mu) X(kappa) X(zeta) X(n)
#include "Model/Solid/InClassOperators.H"
};

#include "Model/Solid/ExtClassOperators.H"

} // namespace Finite
} // namespace Solid
} // namespace Model

#endif // MODEL_SOLID_FINITE_ADHESION_H_