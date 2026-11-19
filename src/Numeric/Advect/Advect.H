#ifndef NUMERIC_ADVECT_ADVECT_H
#define NUMERIC_ADVECT_ADVECT_H

#include <string>
#include <tuple>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{

enum class Form
{
    Advective,
    Conservative
};

struct Options
{
    Form form = Form::Advective;
};


template<typename... OBJ>
class Advect
{
private:
    friend IO::ParmParse;

    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:
    template<int I = 0, class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()( Phi const& phi,
                            Velocity const& velocity,
                            int i, int j, int k,
                            int comp,
                            const Set::Scalar* DX,
                            Options options = {},
                            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType()) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj)(phi, velocity, i, j, k, comp, DX, options, stencil);
            else
                return this->template operator()<I + 1>(phi, velocity, i, j, k, comp, DX, options, stencil);
        }
        else
        {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0, class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar Scalar( Phi const& phi,
                        Velocity const& velocity,
                        int i, int j, int k,
                        int comp,
                        const Set::Scalar* DX,
                        Options options = {},
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType()) const
    {
        return this->template operator()<I>(phi, velocity, i, j, k, comp, DX, options, stencil);
    }

    template<class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Vector Vector( Phi const& phi,
                        Velocity const& velocity,
                        int i, int j, int k,
                        int comp,
                        const Set::Scalar* DX,
                        Options options = {},
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType()) const
    {
        Set::Vector ret = Set::Vector::Zero();
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            ret(d) = (*this)(phi, velocity, i, j, k, comp + d, DX, options, stencil);
        return ret;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    int NGhost() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).NGhost();
            else
                return NGhost<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return 0;
        }
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool IsDefined() const
    {
        return selected >= 0;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::HC PhiLocation() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).phi_location;
            else
                return PhiLocation<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return Set::HC::Cell;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::HC VelocityLocation() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).velocity_location;
            else
                return VelocityLocation<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return Set::HC::Cell;
        }
    }
};

}
}

#endif
