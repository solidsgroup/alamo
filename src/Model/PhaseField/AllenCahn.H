//
// This header defines an Allen-Cahn evolution operator for an order parameter
// eta with free energy
//
// .. math::
//
//    F[\eta] = \int_\Omega \left[\lambda w(\eta)
//              + \frac{\kappa}{2}|\nabla\eta|^2\right]d\mathbf{x}.
//
// The quartic chemical potential w is specified by its values w0, w12, and w1
// at eta=0, eta=1/2, and eta=1, with zero slope at eta=0 and eta=1.
// The resulting nonconservative source is
//
// .. math::
//
//    \dot{\eta} = -L\left[\lambda w'(\eta)-\kappa\Delta\eta\right].
//

#ifndef MODEL_PHASEFIELD_ALLENCAHN_H
#define MODEL_PHASEFIELD_ALLENCAHN_H

#include "AMReX_Array4.H"
#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"
#include "Unit/Unit.H"
#include "Util/Util.H"

namespace Model
{
namespace PhaseField
{
class AllenCahn
{
public:
    static constexpr const char* name = "allencahn";

    static void Parse(AllenCahn& value, IO::ParmParse& pp)
    {
        const bool configured =
            pp.contains("mobility") ||
            pp.contains("lambda") || pp.contains("kappa") ||
            pp.contains("sigma") || pp.contains("epsilon") ||
            pp.contains("w0") || pp.contains("w12") || pp.contains("w1");
        if (!configured) return;

        pp.query_required("mobility", value.mobility,
                        1.0 / Unit::Pressure() / Unit::Time());

        std::pair<std::string, Set::Scalar> coefficients[2];
        pp.query_exactly<2>({"lambda", "kappa", "sigma", "epsilon"}, coefficients,
                            {Unit::Pressure(), Unit::Energy() / Unit::Length(),
                            Unit::Energy() / Unit::Area(), Unit::Length()});

        if (coefficients[0].first == "lambda" && coefficients[1].first == "kappa")
        {
            value.lambda = coefficients[0].second;
            value.kappa = coefficients[1].second;
        }
        else if (coefficients[0].first == "sigma" && coefficients[1].first == "epsilon")
        {
            const Set::Scalar sigma = coefficients[0].second;
            const Set::Scalar epsilon = coefficients[1].second;
            value.lambda = 4.0 * sigma / epsilon;
            value.kappa = 0.5 * sigma * epsilon;
        }
        else
        {
            Util::Exception(INFO,
                "Allen-Cahn requires either (lambda,kappa) or (sigma,epsilon)");
        }

        if (value.mobility < 0.0)
            Util::Exception(INFO, "mobility must be nonnegative");
        if (value.lambda < 0.0 || value.kappa < 0.0)
            Util::Exception(INFO, "Allen-Cahn energy coefficients must be nonnegative");

        Set::Scalar w0, w12, w1;
        pp.query_default("w0", w0, 0.0);
        pp.query_default("w12", w12, 0.0625);
        pp.query_default("w1", w1, 0.0);

        value.a2 = -11.0 * w0 + 16.0 * w12 - 5.0 * w1;
        value.a3 =  18.0 * w0 - 32.0 * w12 + 14.0 * w1;
        value.a4 =  -8.0 * w0 + 16.0 * w12 -  8.0 * w1;

        value.potential_curvature = Util::Max(
            Util::Abs(2.0 * value.a2),
            Util::Abs(2.0 * value.a2 + 6.0 * value.a3 + 12.0 * value.a4));
        if (value.a4 != 0.0)
        {
            const Set::Scalar eta = -value.a3 / (4.0 * value.a4);
            if (eta > 0.0 && eta < 1.0)
                value.potential_curvature = Util::Max(
                    value.potential_curvature,
                    Util::Abs(2.0 * value.a2 + 6.0 * value.a3 * eta +
                            12.0 * value.a4 * eta * eta));
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Local(const amrex::Array4<const Set::Scalar>& field,
                      const int component, const Set::Scalar scale,
                      const int i, const int j, const int k) const
    {
        const Set::Scalar value = scale * field(i,j,k,component);
        const Set::Scalar value2 = value * value;
        return -mobility * lambda *
            (2.0 * a2 * value + 3.0 * a3 * value2 +
             4.0 * a4 * value2 * value);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const amrex::Array4<const Set::Scalar>& field,
                           const int component, const Set::Scalar scale,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        return Local(field, component, scale, i, j, k) +
            GradientCoefficient() * scale *
                Numeric::Laplacian(field, i, j, k, component, DX);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient() const
    {
        return mobility * kappa;
    }

    Set::Scalar LocalStabilityRate() const
    {
        return mobility * lambda * potential_curvature;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx) const
    {
        return LocalStabilityRate() +
            4.0 * AMREX_SPACEDIM * GradientCoefficient() / (dx * dx);
    }

private:
    Set::Scalar mobility = 0.0;
    Set::Scalar lambda = 0.0;
    Set::Scalar kappa = 0.0;
    Set::Scalar a2 = 1.0;
    Set::Scalar a3 = -2.0;
    Set::Scalar a4 = 1.0;
    Set::Scalar potential_curvature = 2.0;
};
}
}

#endif
