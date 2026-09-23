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

// LowMach's pressure projection currently constructs face velocities with an
// arithmetic average of the two neighboring cell-centered values.  Every
// advection scheme declares its face-velocity interpolation so LowMach can
// warn when that projection/advection compatibility assumption is not met.
enum class FaceVelocityInterpolation
{
    ArithmeticCellAverage,
    Other
};

struct Options
{
    Form form = Form::Advective;
    // When supplied, these are the actual projected normal face velocities.
    // A face (i,j,k,d) separates cells (i-e_d) and i.
    amrex::GpuArray<amrex::Array4<const Set::Scalar>, AMREX_SPACEDIM> faces{};
    bool use_faces = false;
};


template<class Velocity>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Set::Scalar FaceVelocity(Velocity const& velocity, int i, int j, int k,
                         int d, Options const& options)
{
    if (options.use_faces) return options.faces[d](i,j,k);
    return 0.5 * (velocity(i-(d==0),j-(d==1),k-(d==2),d) +
                  velocity(i,j,k,d));
}


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

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    FaceVelocityInterpolation FaceVelocityInterpolationType() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).face_velocity_interpolation;
            return FaceVelocityInterpolationType<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return FaceVelocityInterpolation::Other;
        }
    }
};

}
}

#endif
