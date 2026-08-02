//
// This header defines an Allen-Cahn evolution operator for an order parameter
// eta with free energy
//
// .. math::
//
//    F[\eta] = \int_\Omega \left[\lambda w(\eta)
//              + \frac{\kappa}{2}|\nabla\eta|^2\right]d\mathbf{x}.
//
// The chemical potential w is specified by its values w0, w12, and w1
// at eta=0, eta=1/2, and eta=1, with zero slope at eta=0 and eta=1 (so
// eta=0 and eta=1 are always local minima, regardless of the other
// parameters below, including the barrier skew). The optional barrier skew
// (w12_skew, default 0) tilts the interior maximum away from eta=1/2
// without moving w(1/2)=w12 itself -- it only redistributes where that
// same bump amplitude sits, sliding its peak toward eta=0 or eta=1. At
// skew=0 this reduces exactly to the original quartic; a nonzero skew
// makes w a quintic. See Parse() below. The resulting nonconservative
// source is
//
// .. math::
//
//    \dot{\eta} = -L\left[\lambda w'(\eta)-\kappa\Delta\eta\right].
//

#ifndef MODEL_PHASEFIELD_ALLENCAHN_H
#define MODEL_PHASEFIELD_ALLENCAHN_H

#include <cmath>

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
            pp.contains("w0") || pp.contains("w12") || pp.contains("w1") ||
            pp.contains("w12_skew");
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

        Set::Scalar w0, w12, w1, w12_skew;
        pp.query_default("w0", w0, 0.0);
        pp.query_default("w12", w12, 0.0625);
        pp.query_default("w1", w1, 0.0);
        pp.query_default("w12_skew", w12_skew, 0.0);

        // w(eta) = w0 + (w1-w0)*h(eta) + A*f(eta)*(1 + skew*(1-2*eta))
        // with h(eta) = 3eta^2 - 2eta^3 (smoothstep: h(0)=0, h(1)=1,
        // h'(0)=h'(1)=0) and f(eta) = eta^2*(1-eta)^2 (bump: f(0)=f(1)=0,
        // f'(0)=f'(1)=0). Both h and f -- and hence w -- have zero value
        // and zero slope at eta=0 and eta=1 for *any* skew, so the two
        // minima stay pinned there regardless of skew. A is fixed so that
        // w(1/2)=w12 for any skew (the (1-2*eta) tilt vanishes at
        // eta=1/2, so w(1/2) is skew-independent); skew only slides the
        // location of the interior maximum away from eta=1/2. At skew=0
        // this reduces exactly to the original quartic.
        const Set::Scalar A = 16.0 * (w12 - 0.5 * (w0 + w1));

        value.a2 =  3.0 * (w1 - w0) + A * (1.0 + w12_skew);
        value.a3 = -2.0 * (w1 - w0) - A * (2.0 + 4.0 * w12_skew);
        value.a4 =  A * (1.0 + 5.0 * w12_skew);
        value.a5 = -2.0 * A * w12_skew;

        // Max |w''(eta)| over eta in [0,1]: check the endpoints plus any
        // interior roots of w'''(eta) (a quadratic in eta, since w'' is a
        // cubic once a5 != 0).
        auto wpp = [&](Set::Scalar eta) {
            return 2.0 * value.a2 + 6.0 * value.a3 * eta +
                12.0 * value.a4 * eta * eta + 20.0 * value.a5 * eta * eta * eta;
        };
        value.potential_curvature = Util::Max(
            Util::Abs(wpp(0.0)), Util::Abs(wpp(1.0)));

        const Set::Scalar qa = 60.0 * value.a5;
        const Set::Scalar qb = 24.0 * value.a4;
        const Set::Scalar qc =  6.0 * value.a3;
        if (Util::Abs(qa) < 1.0e-14)
        {
            if (Util::Abs(qb) > 1.0e-14)
            {
                const Set::Scalar eta = -qc / qb;
                if (eta > 0.0 && eta < 1.0)
                    value.potential_curvature =
                        Util::Max(value.potential_curvature, Util::Abs(wpp(eta)));
            }
        }
        else
        {
            const Set::Scalar disc = qb * qb - 4.0 * qa * qc;
            if (disc >= 0.0)
            {
                const Set::Scalar sq = std::sqrt(disc);
                const Set::Scalar eta1 = (-qb + sq) / (2.0 * qa);
                const Set::Scalar eta2 = (-qb - sq) / (2.0 * qa);
                if (eta1 > 0.0 && eta1 < 1.0)
                    value.potential_curvature =
                        Util::Max(value.potential_curvature, Util::Abs(wpp(eta1)));
                if (eta2 > 0.0 && eta2 < 1.0)
                    value.potential_curvature =
                        Util::Max(value.potential_curvature, Util::Abs(wpp(eta2)));
            }
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Local(const amrex::Array4<const Set::Scalar>& field,
                      const int component, const Set::Scalar scale,
                      const int i, const int j, const int k) const
    {
        const Set::Scalar value = scale * field(i,j,k,component);
        const Set::Scalar value2 = value * value;
        const Set::Scalar value3 = value2 * value;
        return -mobility * lambda *
            (2.0 * a2 * value + 3.0 * a3 * value2 +
             4.0 * a4 * value3 + 5.0 * a5 * value3 * value);
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate() const
    {
        return mobility * lambda * potential_curvature;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
    Set::Scalar a5 = 0.0;
    Set::Scalar potential_curvature = 2.0;
};
}
}

#endif
