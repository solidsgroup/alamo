#ifndef MODEL_PHASEFIELD_PHASEFIELD_H
#define MODEL_PHASEFIELD_PHASEFIELD_H

#include <tuple>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace PhaseField
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar H(Set::Scalar eta)
{
    // For eta=(1-tanh(2s/ell))/2, these bounds confine the
    // constitutive transition to |s|<ell/2.
    constexpr Set::Scalar eta_lo = 0.11920292202211755;
    constexpr Set::Scalar eta_hi = 1.0 - eta_lo;
    return Util::SmootherStep((eta - eta_lo) / (eta_hi - eta_lo));
}

template<typename... OBJ>
class PhaseField
{
    friend IO::ParmParse;

public:
    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Local(const amrex::Array4<const Set::Scalar>& field,
                      const int component, const Set::Scalar scale,
                      const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).Local(
                    field, component, scale, i, j, k);
            return this->template Local<I + 1>(
                field, component, scale, i, j, k);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const amrex::Array4<const Set::Scalar>& field,
                            const int component, const Set::Scalar scale,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj)(field, component, scale, i, j, k, DX);
            return this->template operator()<I + 1>(
                field, component, scale, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).GradientCoefficient();
            return GradientCoefficient<I + 1>();
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    Set::Scalar LocalStabilityRate() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).LocalStabilityRate();
            return LocalStabilityRate<I + 1>();
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    Set::Scalar StabilityRate(const Set::Scalar dx) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).StabilityRate(dx);
            return StabilityRate<I + 1>(dx);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

private:
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = {OBJ::name...};
    int selected = -1;
};
}
}

#endif
