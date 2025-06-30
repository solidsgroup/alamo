#ifndef UTIL_DISTRIBUTION_H
#define UTIL_DISTRIBUTION_H

#include "Set/Set.H"

#include "AMReX_GpuQualifiers.H"

namespace Util
{
namespace Distribution
{
/// Generate a gamma-distributed random number using Marsaglia-Tsang method
AMREX_GPU_DEVICE
inline Set::Scalar Gamma(Set::Scalar shape, Set::Scalar scale = 1.0)
{
    if (shape <= 0.0)
        return 0.0;

    if (shape < 1.0)
    {
        // For shape < 1, use Johnk's generator
        Set::Scalar u = Util::Random();
        return scale * Gamma(shape + 1.0, 1.0) * std::pow(u, 1.0 / shape);
    }

    // Marsaglia-Tsang method for shape >= 1
    Set::Scalar d = shape - 1.0 / 3.0;
    Set::Scalar c = 1.0 / std::sqrt(9.0 * d);

    while (true)
    {
        Set::Scalar x, v;
        do
        {
            x = Util::Gaussian(0.0, 1.0);
            v = 1.0 + c * x;
        } while (v <= 0.0);

        v = v * v * v;
        Set::Scalar u = Util::Random();

        if (u < 1.0 - 0.0331 * x * x * x * x)
            return scale * d * v;

        if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v)))
            return scale * d * v;
    }
}

/// Generate a beta-distributed random number
AMREX_GPU_DEVICE
inline Set::Scalar Beta(Set::Scalar alpha, Set::Scalar beta)
{
    if (alpha <= 0.0 || beta <= 0.0)
        return 0.5;

    if (alpha == 1.0 && beta == 1.0)
        return Util::Random();

    Set::Scalar x = Gamma(alpha);
    Set::Scalar y = Gamma(beta);

    Set::Scalar sum = x + y;
    if (sum < 1e-30)
        return 0.5;

    return x / sum;
}
}
}
#endif
