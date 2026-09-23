#ifndef MODEL_MECHANISM_SURFACETEMPERATURE_H
#define MODEL_MECHANISM_SURFACETEMPERATURE_H

#include "Set/Set.H"

namespace Model::Mechanism
{
// Opt-in kinetic closures; these do not replace the caloric temperature field.
struct SurfaceTemperatureOptions
{
    enum Mode { Contour, SolidMass, Conversion, ArrheniusEquivalent };
    Mode mode = Contour;
    Set::Scalar contour_fraction = 0.5; // Total condensed volume fraction, not gas fraction.
    Set::Scalar band_cutoff = 0.01;
    Set::Scalar activation_temperature = 0.0;
    int solid_component = -1;
};

// Constant-cp phase transfer at a common enthalpy reference temperature.
// mass_change is the signed change in the reactant condensed density.
// Composition changes C as well as T; Q alone must not also create sensible H.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar ReferenceEnthalpyDeltaT(Set::Scalar temperature, Set::Scalar heat,
    Set::Scalar capacity, Set::Scalar mass_change,
    Set::Scalar reactant_minus_product_cp, Set::Scalar reference_temperature)
{
    const Set::Scalar delta_capacity = mass_change * reactant_minus_product_cp;
    const Set::Scalar new_capacity = capacity + delta_capacity;
    if (!(new_capacity > 0.0) || !std::isfinite(new_capacity)) return NAN;
    return (heat - delta_capacity * (temperature-reference_temperature)) / new_capacity;
}

// Frozen normal-line fields: total condensed fraction, caloric temperature,
// predicted converted mass density, and reacting-solid partial mass density.
struct SurfaceTemperatureRay
{
    amrex::Array4<const Set::Scalar> field;
    amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM> dx;
    amrex::IntVect lower, upper;
    amrex::GpuArray<int, AMREX_SPACEDIM> periodic;
    Set::Vector direction;
    int i, j, k;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(Set::Scalar distance, int component) const
    {
        int base[3] = {0,0,0};
        Set::Scalar alpha[AMREX_SPACEDIM];
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            Set::Scalar offset = distance * direction(d) / dx[d];
            if (periodic[d])
            {
                const int period = upper[d] - lower[d] + 1;
                offset -= period * std::floor(offset / period + 0.5);
            }
            Set::Scalar coordinate = (d == 0 ? i : (d == 1 ? j : k)) + offset;
            if (!periodic[d])
                coordinate = amrex::min(amrex::max(coordinate,
                    Set::Scalar(lower[d])), Set::Scalar(upper[d]));
            base[d] = static_cast<int>(std::floor(coordinate));
            alpha[d] = coordinate - base[d];
        }
        Set::Scalar value = 0.0;
        for (int corner = 0; corner < (1 << AMREX_SPACEDIM); ++corner)
        {
            int index[3] = {base[0],base[1],base[2]};
            Set::Scalar weight = 1.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                const int bit = (corner >> d) & 1;
                index[d] += bit;
                weight *= bit ? alpha[d] : 1.0-alpha[d];
            }
            value += weight * field(index[0],index[1],index[2],component);
        }
        return value;
    }
};

// First exit from the declared diffuse band on each side of the contour.
// Missing support is an error, not an implicit hot-gas or cold-solid fallback.
template<class Ray>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar WeightedSurfaceTemperature(const Ray& ray, Set::Scalar contour,
    Set::Scalar contour_temperature, Set::Scalar step, Set::Scalar reach,
    const SurfaceTemperatureOptions options)
{
    Set::Scalar ends[2] = {NAN,NAN};
    for (int side = 0; side < 2; ++side)
    {
        const int sign = side == 0 ? -1 : 1;
        Set::Scalar previous = 0.0;
        for (Set::Scalar distance = step; distance <= reach; distance += step)
        {
            const Set::Scalar eta = ray(contour + sign*distance, 0);
            if (eta <= options.band_cutoff || eta >= 1.0-options.band_cutoff)
            {
                Set::Scalar a = previous, b = distance;
                for (int iteration = 0; iteration < 24; ++iteration)
                {
                    const Set::Scalar middle = 0.5*(a+b);
                    const Set::Scalar f = ray(contour + sign*middle, 0);
                    if (f > options.band_cutoff && f < 1.0-options.band_cutoff)
                        a = middle;
                    else b = middle;
                }
                ends[side] = contour + sign*0.5*(a+b);
                break;
            }
            previous = distance;
        }
    }
    if (!std::isfinite(ends[0]) || !std::isfinite(ends[1])) return NAN;
    const int count = amrex::max(1, static_cast<int>(std::ceil((ends[1]-ends[0])/step)));
    Set::Scalar sum = 0.0, sum_temperature = 0.0;
    Set::Scalar log_rate_sum = -INFINITY;
    for (int n = 0; n <= count; ++n)
    {
        const Set::Scalar position = ends[0] + (ends[1]-ends[0])*n/count;
        const Set::Scalar temperature = ray(position, 1);
        const int component = options.mode == SurfaceTemperatureOptions::SolidMass ? 3 : 2;
        const Set::Scalar weight = amrex::max(Set::Scalar(0), ray(position, component)) *
            ((n == 0 || n == count) ? 0.5 : 1.0);
        if (!(weight > 0.0)) continue;
        if (!(temperature > 0.0) || !std::isfinite(temperature)) return NAN;
        sum += weight;
        sum_temperature += weight * temperature;
        if (options.mode == SurfaceTemperatureOptions::ArrheniusEquivalent &&
            options.activation_temperature > 0.0)
        {
            const Set::Scalar term = std::log(weight) - options.activation_temperature/temperature;
            if (!std::isfinite(log_rate_sum)) log_rate_sum = term;
            else
            {
                const Set::Scalar high = amrex::max(term, log_rate_sum);
                log_rate_sum = high + std::log(std::exp(term-high) + std::exp(log_rate_sum-high));
            }
        }
    }
    // An identically inactive predictor (e.g. A=0) has no conversion average.
    if (!(sum > 0.0)) return contour_temperature;
    if (options.mode == SurfaceTemperatureOptions::ArrheniusEquivalent &&
        options.activation_temperature > 0.0)
        return -options.activation_temperature/(log_rate_sum-std::log(sum));
    return sum_temperature/sum;
}
}
#endif
