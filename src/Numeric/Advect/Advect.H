#ifndef NUMERIC_ADVECT_ADVECT_H
#define NUMERIC_ADVECT_ADVECT_H

#include <array>
#include <cmath>
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

struct FaceStates
{
    Set::Scalar left = 0.0;
    Set::Scalar right = 0.0;
};

template<int N>
struct CompressibleFlux
{
    std::array<Set::Scalar, N> mass{};
    Set::Vector momentum = Set::Vector::Zero();
    Set::Scalar energy = 0.0;
};


template<typename... OBJ>
class Advect
{
private:
    friend IO::ParmParse;

    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

    template<int N, class Scheme, class Rho, class Momentum, class Energy,
             class Density, class Velocity, class Pressure, class Temperature,
             class MoleFraction, class Gas>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    CompressibleFlux<N> ComputeCompressibleFaceFlux(
        Scheme const& scheme,
        Rho const& rho,
        Momentum const& momentum,
        Energy const& energy,
        Density const& density,
        Velocity const& velocity,
        Pressure const& pressure,
        Temperature const& temperature,
        MoleFraction mole_fraction,
        Gas const& gas,
        int i, int j, int k,
        int direction,
        Set::Scalar small) const
    {
        const int di = direction == 0 ? 1 : 0;
        const int dj = direction == 1 ? 1 : 0;
        const int dk = direction == 2 ? 1 : 0;
        const int ir = i + di;
        const int jr = j + dj;
        const int kr = k + dk;

        FaceStates rho_face = scheme.ReconstructFace(
            density, i, j, k, direction, 0);
        FaceStates pressure_face = scheme.ReconstructFace(
            pressure, i, j, k, direction, 0);
        FaceStates velocity_face = scheme.ReconstructFace(
            velocity, i, j, k, direction, direction);

        const Set::Scalar rho_left = Util::Max(small, rho_face.left);
        const Set::Scalar rho_right = Util::Max(small, rho_face.right);
        const Set::Scalar pressure_left = Util::Max(small, pressure_face.left);
        const Set::Scalar pressure_right = Util::Max(small, pressure_face.right);
        const Set::Scalar gamma_left = Util::Max(
            small,
            gas.gamma(temperature(i,j,k), mole_fraction, i, j, k));
        const Set::Scalar gamma_right = Util::Max(
            small,
            gas.gamma(temperature(ir,jr,kr), mole_fraction, ir, jr, kr));

        const Set::Scalar impedance_left = std::sqrt(
            gamma_left * pressure_left * rho_left);
        const Set::Scalar impedance_right = std::sqrt(
            gamma_right * pressure_right * rho_right);
        const Set::Scalar impedance_sum = Util::Max(
            small, impedance_left + impedance_right);

        const Set::Scalar face_velocity =
            (impedance_left * velocity_face.left +
             impedance_right * velocity_face.right +
             pressure_left - pressure_right) / impedance_sum;
        const Set::Scalar face_pressure = Util::Max(
            small,
            (impedance_right * pressure_left +
             impedance_left * pressure_right +
             impedance_left * impedance_right *
                 (velocity_face.left - velocity_face.right)) / impedance_sum);

        CompressibleFlux<N> flux;
        // The acoustic solve above needs both face states.  Once its velocity
        // is known, upwind schemes only need the selected transported state.
        for (int n = 0; n < N; ++n)
        {
            flux.mass[n] = face_velocity *
                scheme.ReconstructTransportedFaceValue(
                    rho, i, j, k, direction, n, face_velocity);
        }

        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            flux.momentum(d) = face_velocity *
                scheme.ReconstructTransportedFaceValue(
                    momentum, i, j, k, direction, d, face_velocity);
        }
        flux.momentum(direction) += face_pressure;

        flux.energy = face_velocity *
            scheme.ReconstructTransportedFaceValue(
                energy, i, j, k, direction, 0, face_velocity) +
            face_pressure * face_velocity;

        return flux;
    }

    template<int N, class Scheme, class Rho, class Momentum, class Energy,
             class Density, class Velocity, class Pressure, class Temperature,
             class MoleFraction, class Gas>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    CompressibleFlux<N> ComputeCompressibleFluxDivergenceForScheme(
        Scheme const& scheme,
        Rho const& rho,
        Momentum const& momentum,
        Energy const& energy,
        Density const& density,
        Velocity const& velocity,
        Pressure const& pressure,
        Temperature const& temperature,
        MoleFraction mole_fraction,
        Gas const& gas,
        int i, int j, int k,
        Set::Scalar eta,
        const Set::Scalar* DX,
        Set::Scalar small) const
    {
        CompressibleFlux<N> result;

        const CompressibleFlux<N> flux_xlo = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i-1, j, k, 0, small);
        const CompressibleFlux<N> flux_xhi = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i, j, k, 0, small);

        for (int n = 0; n < N; ++n)
        {
            result.mass[n] = (flux_xlo.mass[n] - flux_xhi.mass[n]) / DX[0];
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            result.momentum(d) =
                (flux_xlo.momentum(d) - flux_xhi.momentum(d)) / DX[0];
        }
        result.energy = (flux_xlo.energy - flux_xhi.energy) / DX[0];

#if AMREX_SPACEDIM >= 2
        const CompressibleFlux<N> flux_ylo = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i, j-1, k, 1, small);
        const CompressibleFlux<N> flux_yhi = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i, j, k, 1, small);

        for (int n = 0; n < N; ++n)
        {
            result.mass[n] += (flux_ylo.mass[n] - flux_yhi.mass[n]) / DX[1];
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            result.momentum(d) +=
                (flux_ylo.momentum(d) - flux_yhi.momentum(d)) / DX[1];
        }
        result.energy += (flux_ylo.energy - flux_yhi.energy) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
        const CompressibleFlux<N> flux_zlo = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i, j, k-1, 2, small);
        const CompressibleFlux<N> flux_zhi = ComputeCompressibleFaceFlux<N>(
            scheme, rho, momentum, energy, density, velocity, pressure,
            temperature, mole_fraction, gas, i, j, k, 2, small);

        for (int n = 0; n < N; ++n)
        {
            result.mass[n] += (flux_zlo.mass[n] - flux_zhi.mass[n]) / DX[2];
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            result.momentum(d) +=
                (flux_zlo.momentum(d) - flux_zhi.momentum(d)) / DX[2];
        }
        result.energy += (flux_zlo.energy - flux_zhi.energy) / DX[2];
#endif

        for (int n = 0; n < N; ++n) result.mass[n] *= eta;
        result.momentum *= eta;
        result.energy *= eta;
        return result;
    }

public:
    template<int N, int I = 0, class Rho, class Momentum, class Energy,
             class Density, class Velocity, class Pressure, class Temperature,
             class MoleFraction, class Gas>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    CompressibleFlux<N> ComputeCompressibleFluxDivergence(
        Rho const& rho,
        Momentum const& momentum,
        Energy const& energy,
        Density const& density,
        Velocity const& velocity,
        Pressure const& pressure,
        Temperature const& temperature,
        MoleFraction mole_fraction,
        Gas const& gas,
        int i, int j, int k,
        Set::Scalar eta,
        const Set::Scalar* DX,
        Set::Scalar small) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                return ComputeCompressibleFluxDivergenceForScheme<N>(
                    std::get<I>(obj), rho, momentum, energy, density, velocity,
                    pressure, temperature, mole_fraction, gas, i, j, k, eta,
                    DX, small);
            }
            return ComputeCompressibleFluxDivergence<N, I + 1>(
                rho, momentum, energy, density, velocity, pressure,
                temperature, mole_fraction, gas, i, j, k, eta, DX, small);
        }
        else
        {
            Util::Abort(INFO, selected);
            return {};
        }
    }

    template<int I = 0, class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()( Phi const& phi,
                            Velocity const& velocity,
                            int i, int j, int k,
                            int comp,
                            const Set::Scalar* DX,
                            Options options = {},
                            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType) const
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
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType) const
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
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType) const
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
