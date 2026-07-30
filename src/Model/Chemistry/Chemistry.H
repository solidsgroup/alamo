//
// This header defines the :code:Model::Chemistry::Chemistry class, which introduces a static
// dispatch mechanism to replicate traditional dynamic polymorphism without incurring the performance
// cost of virtual table lookups.
//
// The :code:Chemistry class acts as the static equivalent of a pure abstract interface for
// chemical reaction models. The function :code:ComputeChemistrySources statically dispatches to the
// appropriate underlying implementation, encapsulating key properties for reaction rates
// and heat release. The :code:ComputeChemistrySources function is resolved at compile-time for GPU
// compatibility and runtime performance, and returns the species production rates and heat release.
//
// This system enables polymorphic behavior without relying on runtime virtual functions, making it
// suitable for high-performance computing environments.
//

#ifndef MODEL_CHEMISTRY_H
#define MODEL_CHEMISTRY_H

#include <cmath>
#include <limits>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Model/Gas/Gas.H"
#include "Solver/Local/ODE/BackwardEuler.H"
#include "Solver/Local/ODE/ForwardEuler.H"
#include "Solver/Local/ODE/ODE.H"

namespace Model
{
namespace Chemistry
{
static constexpr int MAX_SPECIES = 32;
using SpeciesArray = std::array<Set::Scalar, MAX_SPECIES>;
using Source = std::pair<SpeciesArray, Set::Scalar>;

struct Advance
{
    struct Result
    {
        bool converged = true;
        Set::Scalar volume_fraction_change = 0.0;
        Set::Scalar failed_dt = 0.0;
        Set::Scalar residual_norm = 0.0;
    };
};

template<typename... OBJ>
class Chemistry
{
public:
    using ODESolver = Solver::Local::ODE::ODE<
        Solver::Local::ODE::ForwardEuler,
        Solver::Local::ODE::BackwardEuler>;

private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;
    ODESolver solver;
    bool split_chemistry = false;
    int nspecies = 0;

public:
    const char* model_name() const
    {
        return names[selected];
    }

    void Define(int a_nspecies)
    {
        nspecies = a_nspecies;
    }

    template<typename T>
    const T& Get() const
    {
        return std::get<T>(obj);
    }

    const ODESolver& Solver() const
    {
        return solver;
    }

    static void Parse(Chemistry& value, IO::ParmParse& pp)
    {
        Util::Assert(INFO, TEST(value.nspecies > 0 && value.nspecies <= MAX_SPECIES),
            "Chemistry requires between 1 and ", MAX_SPECIES, " gas species.");
        value.split_chemistry = false;
        pp.select<OBJ...>("model", value, pp.forward_args(value.nspecies));
        pp.select<
            Solver::Local::ODE::ForwardEuler,
            Solver::Local::ODE::BackwardEuler>("solver", value.solver);
        value.split_chemistry = value.Reactive();

        if (value.split_chemistry && value.solver.RequiresJacobian() &&
            !value.HasJacobian())
            Util::Exception(INFO, "The selected chemistry model does not provide the "
                "Jacobian required by the selected ODE solver.");
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool Reactive() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::tuple_element_t<I, std::tuple<OBJ...>>::reactive;
            return Reactive<I+1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return false;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool HasJacobian() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::tuple_element_t<I, std::tuple<OBJ...>>::has_jacobian;
            return HasJacobian<I + 1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return false;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool ComputeODEJacobian(
        const Set::Scalar* state,
        Set::Scalar* derivative,
        Set::Scalar pressure,
        Set::Scalar mixture_density,
        Set::Scalar gas_density,
        int nspecies,
        int dependent_species,
        Set::Scalar temperature_scale,
        const Model::Gas::Gas::DeviceData* gas) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            using Model = std::tuple_element_t<I, std::tuple<OBJ...>>;
            if (selected == I)
            {
                if constexpr (Model::has_jacobian)
                    return std::get<I>(obj).ComputeODEJacobian(
                        state, derivative, pressure, mixture_density,
                        gas_density, nspecies, dependent_species,
                        temperature_scale, gas);
                else
                    return false;
            }
            return ComputeODEJacobian<I + 1>(
                state, derivative, pressure, mixture_density, gas_density,
                nspecies, dependent_species, temperature_scale, gas);
        }
        else
        {
            Util::Abort(INFO,selected);
            return false;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Source
    ComputeChemistrySources (const Set::Scalar P, const Set::Scalar T, const SpeciesArray& rhoY,
                const Set::Scalar dt, const Model::Gas::Gas* gas = nullptr) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).ComputeChemistrySources(P,T,rhoY,dt,gas);
            else return ComputeChemistrySources<I+1>(P,T,rhoY,dt,gas);
        }
        else
        {
            Util::Abort(INFO,selected);
            return {};
        }
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool Split() const
    {
        return split_chemistry;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool Implicit() const
    {
        return split_chemistry && solver.Implicit();
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Model::Chemistry::Advance::Result
    Advance(
        const Set::Scalar dt,
        const Set::Scalar pressure,
        const Set::Scalar mixture_density,
        const int nspecies,
        SpeciesArray& rhoY,
        Set::Scalar& temperature,
        const Model::Gas::Gas* gas) const
    {
        Model::Chemistry::Advance::Result advance;
#ifdef ALAMO_GPU
        Util::IgnoreUnused(
            dt, pressure, mixture_density, nspecies, rhoY, temperature, gas);
        advance.converged = false;
        return advance;
#else
        if (!(dt > 0.0)) return advance;
        if (nspecies <= 0 || nspecies > MAX_SPECIES || gas == nullptr ||
            !(pressure > 0.0) || !(mixture_density > 0.0) ||
            !(temperature > 0.0))
        {
            advance.converged = false;
            advance.residual_norm = std::numeric_limits<Set::Scalar>::infinity();
            return advance;
        }

        SpeciesArray mass_fraction{};
        Set::Scalar gas_density = 0.0;
        int dependent_species = 0;
        for (int n = 0; n < nspecies; ++n)
        {
            rhoY[n] = rhoY[n] > 0.0 ? rhoY[n] : 0.0;
            gas_density += rhoY[n];
            if (rhoY[n] > rhoY[dependent_species]) dependent_species = n;
        }
        if (!(gas_density > 0.0))
        {
            advance.converged = false;
            advance.residual_norm = std::numeric_limits<Set::Scalar>::infinity();
            return advance;
        }

        Set::Scalar independent_mass_fraction = 0.0;
        for (int n = 0; n < nspecies; ++n)
        {
            mass_fraction[n] = rhoY[n] / gas_density;
            if (n != dependent_species)
                independent_mass_fraction += mass_fraction[n];
        }
        mass_fraction[dependent_species] = 1.0 - independent_mass_fraction;

        Set::Scalar inverse_mw = 0.0;
        for (int n = 0; n < nspecies; ++n)
            inverse_mw += mass_fraction[n] / gas->MW[n];
        Set::Scalar rho_eos = pressure / (gas->Rg * inverse_mw * temperature);
        amrex::Gpu::DeviceVector<Set::Scalar> gas_device_storage;
        const Model::Gas::Gas::DeviceData gas_data =
            gas->GetDeviceData(gas_device_storage);
        Set::Scalar remaining = dt;
        Set::Scalar chemistry_dt = dt;
        int rejected_steps = 0;

        while (remaining > 0.0)
        {
            const SpeciesArray mass_fraction_old = mass_fraction;
            const Set::Scalar temperature_old = temperature;
            const Set::Scalar rho_eos_old = rho_eos;
            const Set::Scalar temperature_scale =
                temperature_old > 1.0 ? temperature_old : 1.0;
            Set::Scalar state[MAX_SPECIES]{};
            for (int n = 0; n < nspecies - 1; ++n)
            {
                const int species = n < dependent_species ? n : n + 1;
                state[n] = mass_fraction_old[species];
            }
            state[nspecies - 1] = temperature_old / temperature_scale;

            auto rhs = [=,this] AMREX_GPU_DEVICE(
                const Set::Scalar* z, Set::Scalar* rate)
            {
                SpeciesArray Y{};
                Set::Scalar sum = 0.0;
                for (int n = 0; n < nspecies - 1; ++n)
                {
                    const int species = n < dependent_species ? n : n + 1;
                    Y[species] = z[n];
                    if (!(Y[species] >= 0.0) || !std::isfinite(Y[species]))
                        return false;
                    sum += Y[species];
                }
                Y[dependent_species] = 1.0 - sum;
                if (!(Y[dependent_species] >= 0.0) ||
                    !std::isfinite(Y[dependent_species])) return false;

                const Set::Scalar trial_temperature =
                    z[nspecies - 1] * temperature_scale;
                if (!(trial_temperature > 0.0) ||
                    !std::isfinite(trial_temperature)) return false;

                Set::Scalar trial_inverse_mw = 0.0;
                Set::Scalar cp = 0.0;
                for (int n = 0; n < nspecies; ++n)
                {
                    trial_inverse_mw += Y[n] / gas->MW[n];
                    cp += Y[n] * gas->cp_mol_species(trial_temperature, n) /
                        gas->MW[n];
                }
                if (!(trial_inverse_mw > 0.0) || !(cp > 0.0)) return false;

                const Set::Scalar rho_state = pressure /
                    (gas->Rg * trial_inverse_mw * trial_temperature);
                SpeciesArray intrinsic_rhoY{};
                for (int n = 0; n < nspecies; ++n)
                    intrinsic_rhoY[n] = rho_state * Y[n];
                const Source source = ComputeChemistrySources(
                    pressure, trial_temperature, intrinsic_rhoY, 0.0, gas);

                for (int n = 0; n < nspecies - 1; ++n)
                {
                    const int species = n < dependent_species ? n : n + 1;
                    rate[n] = source.first[species] / rho_state;
                }
                rate[nspecies - 1] = gas_density * source.second /
                    (rho_state * mixture_density * cp * temperature_scale);

                for (int n = 0; n < nspecies; ++n)
                    if (!std::isfinite(rate[n])) return false;
                return true;
            };

            auto jacobian = [=,this] AMREX_GPU_DEVICE(
                const Set::Scalar* z, Set::Scalar* derivative)
            {
                return ComputeODEJacobian(
                    z, derivative, pressure, mixture_density, gas_density,
                    nspecies, dependent_species, temperature_scale,
                    &gas_data);
            };

            const Solver::Local::ODE::Advance::Result result =
                solver.template Advance<MAX_SPECIES>(
                    state, nspecies, chemistry_dt, rhs, jacobian,
                    HasJacobian() && solver.SupportsJacobian());
            if (!result.converged)
            {
                chemistry_dt *= 0.5;
                ++rejected_steps;
                if (rejected_steps > 40 ||
                    !(chemistry_dt > dt *
                        std::numeric_limits<Set::Scalar>::epsilon()))
                {
                    advance.converged = false;
                    advance.failed_dt = chemistry_dt;
                    advance.residual_norm = result.residual_norm;
                    return advance;
                }
                continue;
            }

            Set::Scalar sum = 0.0;
            for (int n = 0; n < nspecies - 1; ++n)
            {
                const int species = n < dependent_species ? n : n + 1;
                mass_fraction[species] = state[n];
                sum += state[n];
            }
            mass_fraction[dependent_species] = 1.0 - sum;
            temperature = state[nspecies - 1] * temperature_scale;
            inverse_mw = 0.0;
            for (int n = 0; n < nspecies; ++n)
                inverse_mw += mass_fraction[n] / gas->MW[n];
            rho_eos = pressure / (gas->Rg * inverse_mw * temperature);
            advance.volume_fraction_change +=
                gas_density * (1.0 / rho_eos - 1.0 / rho_eos_old);

            remaining -= chemistry_dt;
            if (remaining <= dt * 1.0e-14) break;
            chemistry_dt = 2.0 * chemistry_dt < remaining ?
                2.0 * chemistry_dt : remaining;
            rejected_steps = 0;
        }

        for (int n = 0; n < nspecies; ++n)
            rhoY[n] = gas_density * mass_fraction[n];
        return advance;
#endif
    }
};
}
}


#endif
