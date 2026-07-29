#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include "AMReX_GpuContainers.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Thermo/Rocfire.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/Transport/Rocfire.H"
#include "Model/Gas/EOS/CPG.H"
#include "Model/Gas/EOS/TPG.H"
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {

/// Cell-local mole fractions reconstructed lazily from gas partial densities.
class MoleFraction
{
public:
    MoleFraction(const amrex::Array4<const Set::Scalar>& a_partial_density,
                const double* a_molecular_weight,
                const int a_nspecies,
                const int i, const int j, const int k)
        : partial_density(a_partial_density),
        molecular_weight(a_molecular_weight),
        moles(0.0)
    {
        for (int n = 0; n < a_nspecies; ++n)
            moles += partial_density(i,j,k,n) / molecular_weight[n];
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()(const int i, const int j, const int k, const int n) const
    {
        return moles > 0.0 ?
            (partial_density(i,j,k,n) / molecular_weight[n]) / moles :
            (n == 0 ? 1.0 : 0.0);
    }

private:
    amrex::Array4<const Set::Scalar> partial_density;
    const double* molecular_weight;
    Set::Scalar moles;
};

class Gas {
protected:
public:
    static constexpr int max_device_species = 32;
    using DeviceArray =
        amrex::GpuArray<Set::Scalar, max_device_species>;
    using DeviceData = amrex::GpuTuple<
        int, Set::Scalar, const Set::Scalar*,
        const Set::Scalar*, const Set::Scalar*, const Set::Scalar*,
        int, int, const Set::Scalar*, const Set::Scalar*,
        Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar,
        Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar,
        const Gas*>;

    int nspecies;
    double Rg;
    std::vector<double> MW;      // Species molecular weights, kg/kmol
    Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire> thermo;
    Transport::Transport<Transport::Mixture_Averaged,Transport::Rocfire> transport;
    EOS::EOS<EOS::CPG,EOS::TPG> eos;

public:
    Gas() 
    {
        this->Rg = Set::Constant::Rg;
    }
    Gas(IO::ParmParse& pp, std::string name) : Gas()
    {
        pp.queryclass(name, *this);
    }

    static void Parse(Gas & value, IO::ParmParse & pp)
    {
        // number of species
        pp.queryarr_required("mw",value.MW,Unit::Mass()/Unit::Amount());
        value.nspecies = value.MW.size();

        // thermo model
        pp.select<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire>("thermo", value.thermo, pp.forward_args(value.nspecies));

        // transport model
        pp.select<Transport::Mixture_Averaged,Transport::Rocfire>("transport", value.transport, pp.forward_args(value.nspecies, value.MW, value.thermo));

        // eos model
        pp.select<EOS::CPG,EOS::TPG>("eos", value.eos);
    }
    
public:
    // Thermodynamic quantities
    double cp_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double entropy_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double cp_mol_species(double T, int n) const;
    double dcp_mol_species_dT(double T, int n) const;
    double enthalpy_mol_species(double T, int n) const;
    double entropy_mol_species(double T, int n) const;

    // Transport quantities
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    void diffusion_coeffs(Set::Patch<Set::Scalar> DKM, double T, double P, Set::Patch<const Set::Scalar> X, int i, int j, int k);

    // EOS
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar> X, int i, int j, int k, double rtol=1e-12) const;
    double ComputeT(
            double pressure, double density,
            Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar> X, int i, int j, int k) const;

    MoleFraction MoleFractions(const amrex::Array4<const Set::Scalar>& partial_density,
                                const int i, const int j, const int k) const
    {
        return MoleFraction(partial_density, MW.data(), nspecies, i, j, k);
    }

    DeviceData GetDeviceData(
        amrex::Gpu::DeviceVector<Set::Scalar>& storage) const
    {
#ifndef ALAMO_GPU
        Util::IgnoreUnused(storage);
        return {nspecies, Rg, nullptr, nullptr, nullptr, nullptr,
                0, 0, nullptr, nullptr,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, this};
#else
        DeviceArray molecular_weight{}, cp_mass{}, h0_mass{}, Tref{};
        DeviceArray transport_value1{}, transport_value2{};
        for (int n = 0; n < nspecies; ++n)
            molecular_weight[n] = MW[n];

        if (std::strcmp(thermo.model_name(), Thermo::CpConstant::name) == 0)
        {
            const auto& [cp_moles, h0, reference_temperature] =
                thermo.Get<Thermo::CpConstant>().DeviceData();
            for (int n = 0; n < nspecies; ++n)
            {
                cp_mass[n] = cp_moles[n] / MW[n];
                h0_mass[n] = h0[n] / MW[n];
                Tref[n] = reference_temperature[n];
            }
        }
        else if (std::strcmp(thermo.model_name(), Thermo::Rocfire::name) == 0)
        {
            const auto& [model_mw, model_cp_mass, reference_temperature] =
                thermo.Get<Thermo::Rocfire>().DeviceData();
            for (int n = 0; n < nspecies; ++n)
            {
                cp_mass[n] = model_cp_mass * model_mw / MW[n];
                h0_mass[n] = 0.0;
                Tref[n] = reference_temperature;
            }
        }
        else
            Util::Exception(INFO,
                "CUDA LowMach supports cpconstant and rocfire thermodynamics");

        int transport_type = -1;
        Set::Scalar diffusion_factor = 0.0;
        Set::Scalar atmospheric_pressure = 0.0;
        Set::Scalar angstrom = 0.0;
        Set::Scalar square_centimeter = 0.0;
        Set::Scalar lewis = 0.0;
        Set::Scalar prandtl = 0.0;
        Set::Scalar lambda_a = 0.0;
        Set::Scalar lambda_b = 0.0;
        if (std::strcmp(
                transport.model_name(), Transport::Mixture_Averaged::name) == 0)
        {
            const auto& [type, value1, value2, d_fact, P_1atm, L_1ang,
                         A_1cmsq] =
                transport.Get<Transport::Mixture_Averaged>().DeviceData();
            transport_type = type == "constant" ? 0 : 1;
            for (int n = 0; n < nspecies; ++n)
            {
                transport_value1[n] = value1[n];
                transport_value2[n] = value2[n];
            }
            diffusion_factor = d_fact;
            atmospheric_pressure = P_1atm;
            angstrom = L_1ang;
            square_centimeter = A_1cmsq;
        }
        else if (std::strcmp(
                     transport.model_name(), Transport::Rocfire::name) == 0)
        {
            transport_type = 2;
            std::tie(lewis, prandtl, lambda_a, lambda_b) =
                transport.Get<Transport::Rocfire>().DeviceData();
        }
        else
            Util::Exception(INFO,
                "CUDA LowMach supports mixture_averaged and rocfire transport");

        const int eos_type =
            std::strcmp(eos.model_name(), EOS::TPG::name) == 0 ? 1 : 0;
        std::array<Set::Scalar, 6 * max_device_species> host{};
        for (int n = 0; n < max_device_species; ++n)
        {
            host[n] = molecular_weight[n];
            host[max_device_species + n] = cp_mass[n];
            host[2 * max_device_species + n] = h0_mass[n];
            host[3 * max_device_species + n] = Tref[n];
            host[4 * max_device_species + n] = transport_value1[n];
            host[5 * max_device_species + n] = transport_value2[n];
        }
        storage.resize(host.size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, host.begin(), host.end(),
            storage.begin());
        const Set::Scalar* device = storage.data();
        return {nspecies, Rg, device,
                device + max_device_species,
                device + 2 * max_device_species,
                device + 3 * max_device_species,
                eos_type, transport_type,
                device + 4 * max_device_species,
                device + 5 * max_device_species,
                diffusion_factor, atmospheric_pressure, angstrom,
                square_centimeter, lewis, prandtl, lambda_a, lambda_b,
                nullptr};
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar GasDensity(
        const DeviceData& data, const Composition& partial_density,
        int i, int j, int k)
    {
        const int count = amrex::get<0>(data);
        Set::Scalar density = 0.0;
        for (int n = 0; n < count; ++n)
            density += partial_density(i,j,k,n);
        return density;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar UniversalGasConstant(const DeviceData& data)
    {
#ifndef ALAMO_GPU
        return amrex::get<18>(data)->Rg;
#else
        return amrex::get<1>(data);
#endif
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar MolecularWeight(
        const DeviceData& data, int species)
    {
#ifndef ALAMO_GPU
        return amrex::get<18>(data)->MW[species];
#else
        return amrex::get<2>(data)[species];
#endif
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar CpMassSpecies(
        const DeviceData& data, Set::Scalar temperature, int species)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        return gas->cp_mol_species(temperature, species) /
            gas->MW[species];
#else
        Util::IgnoreUnused(temperature);
        return amrex::get<3>(data)[species];
#endif
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar DCpMassSpeciesDT(
        const DeviceData& data, Set::Scalar temperature, int species)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        return gas->dcp_mol_species_dT(temperature, species) /
            gas->MW[species];
#else
        Util::IgnoreUnused(data, temperature, species);
        return 0.0;
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar GasConstant(
        const DeviceData& data, const Composition& partial_density,
        int i, int j, int k)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->R(X, i, j, k);
#else
        const auto& [count, universal_gas_constant, molecular_weight,
                     cp_mass, h0_mass, Tref, eos_type, transport_type,
                     value1, value2, diffusion_factor, atmospheric_pressure,
                     angstrom, square_centimeter, lewis, prandtl, lambda_a,
                     lambda_b, host_gas] = data;
        Set::Scalar density = 0.0;
        Set::Scalar moles = 0.0;
        for (int n = 0; n < count; ++n)
        {
            density += partial_density(i,j,k,n);
            moles += partial_density(i,j,k,n) / molecular_weight[n];
        }
        return density > 0.0 ?
            universal_gas_constant * moles / density :
            universal_gas_constant / molecular_weight[0];
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar CpMass(
        const DeviceData& data, Set::Scalar temperature,
        const Composition& partial_density,
        int i, int j, int k)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->cp_mass(temperature, X, i, j, k);
#else
        (void)temperature;
        const int count = amrex::get<0>(data);
        const auto& cp_mass = amrex::get<3>(data);
        Set::Scalar density = 0.0;
        Set::Scalar cp = 0.0;
        for (int n = 0; n < count; ++n)
        {
            const Set::Scalar rhoY = partial_density(i,j,k,n);
            density += rhoY;
            cp += rhoY * cp_mass[n];
        }
        return density > 0.0 ? cp / density : cp_mass[0];
#endif
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar EnthalpyMassSpecies(
        const DeviceData& data, Set::Scalar temperature, int species)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        return gas->enthalpy_mol_species(temperature, species) /
            gas->MW[species];
#else
        return amrex::get<3>(data)[species] *
                   (temperature - amrex::get<5>(data)[species]) +
               amrex::get<4>(data)[species];
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar ThermalConductivity(
        const DeviceData& data, Set::Scalar temperature,
        const Composition& partial_density, int i, int j, int k)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->thermal_conductivity(temperature, X, i, j, k);
#else
        const auto& [count, universal_gas_constant, molecular_weight,
                     cp_mass, h0_mass, Tref, eos_type, transport_type,
                     value1, value2, diffusion_factor, atmospheric_pressure,
                     angstrom, square_centimeter, lewis, prandtl, lambda_a,
                     lambda_b, host_gas] = data;
        if (transport_type == 2)
            return lambda_a * temperature + lambda_b;

        Set::Scalar moles = 0.0;
        for (int n = 0; n < count; ++n)
            moles += partial_density(i,j,k,n) / molecular_weight[n];
        auto mole_fraction = [=] AMREX_GPU_HOST_DEVICE(int n)
        {
            return moles > 0.0 ?
                partial_density(i,j,k,n) /
                    molecular_weight[n] / moles :
                (n == 0 ? 1.0 : 0.0);
        };
        auto viscosity = [=] AMREX_GPU_HOST_DEVICE(int n)
        {
            if (transport_type == 0) return value1[n];
            const Set::Scalar sigma = value1[n] / angstrom;
            const Set::Scalar Tstar = temperature / value2[n];
            const Set::Scalar collision =
                1.16145 * std::pow(Tstar, -0.14874) +
                0.52487 / std::exp(0.77320 * Tstar) +
                2.16178 / std::exp(2.43787 * Tstar);
            return 2.669570e-6 *
                std::sqrt(molecular_weight[n] * temperature) /
                (sigma * sigma * collision);
        };

        Set::Scalar conductivity = 0.0;
        for (int a = 0; a < count; ++a)
        {
            const Set::Scalar mu_a = viscosity(a);
            Set::Scalar phi = 0.0;
            for (int b = 0; b < count; ++b)
            {
                const Set::Scalar ratio =
                    std::pow(molecular_weight[b] / molecular_weight[a], 0.25);
                const Set::Scalar mu_ratio =
                    std::sqrt(mu_a / viscosity(b));
                phi += mole_fraction(b) *
                    std::pow(1.0 + mu_ratio * ratio, 2.0) /
                    (std::sqrt(8.0) *
                     std::sqrt(1.0 + molecular_weight[a] /
                                         molecular_weight[b]));
            }
            const Set::Scalar species_conductivity =
                transport_type == 0 ? value2[a] :
                mu_a * (cp_mass[a] +
                        1.25 * universal_gas_constant /
                            molecular_weight[a]);
            conductivity += mole_fraction(a) * species_conductivity / phi;
        }
        return conductivity;
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar DynamicViscosity(
        const DeviceData& data, Set::Scalar temperature,
        const Composition& partial_density, int i, int j, int k)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->dynamic_viscosity(temperature, X, i, j, k);
#else
        if (amrex::get<7>(data) == 2)
        {
            const Set::Scalar cp =
                CpMass(data, temperature, partial_density, i, j, k);
            return cp > 0.0 ?
                amrex::get<15>(data) *
                    ThermalConductivity(
                        data, temperature, partial_density, i, j, k) / cp :
                0.0;
        }

        const auto& molecular_weight = amrex::get<2>(data);
        const auto& value1 = amrex::get<8>(data);
        const auto& value2 = amrex::get<9>(data);
        const int count = amrex::get<0>(data);
        const int transport_type = amrex::get<7>(data);
        const Set::Scalar angstrom = amrex::get<12>(data);
        Set::Scalar moles = 0.0;
        for (int n = 0; n < count; ++n)
            moles += partial_density(i,j,k,n) / molecular_weight[n];
        auto mole_fraction = [=] AMREX_GPU_HOST_DEVICE(int n)
        {
            return moles > 0.0 ?
                partial_density(i,j,k,n) /
                    molecular_weight[n] / moles :
                (n == 0 ? 1.0 : 0.0);
        };
        auto viscosity = [=] AMREX_GPU_HOST_DEVICE(int n)
        {
            if (transport_type == 0) return value1[n];
            const Set::Scalar sigma = value1[n] / angstrom;
            const Set::Scalar Tstar = temperature / value2[n];
            const Set::Scalar collision =
                1.16145 * std::pow(Tstar, -0.14874) +
                0.52487 / std::exp(0.77320 * Tstar) +
                2.16178 / std::exp(2.43787 * Tstar);
            return 2.669570e-6 *
                std::sqrt(molecular_weight[n] * temperature) /
                (sigma * sigma * collision);
        };
        Set::Scalar mixture_viscosity = 0.0;
        for (int a = 0; a < count; ++a)
        {
            const Set::Scalar mu_a = viscosity(a);
            Set::Scalar phi = 0.0;
            for (int b = 0; b < count; ++b)
            {
                const Set::Scalar ratio =
                    std::pow(molecular_weight[b] / molecular_weight[a], 0.25);
                phi += mole_fraction(b) *
                    std::pow(1.0 + std::sqrt(mu_a / viscosity(b)) * ratio,
                             2.0) /
                    (std::sqrt(8.0) *
                     std::sqrt(1.0 + molecular_weight[a] /
                                         molecular_weight[b]));
            }
            mixture_viscosity += mole_fraction(a) * mu_a / phi;
        }
        return std::isfinite(mixture_viscosity) ? mixture_viscosity : 0.0;
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar DiffusionCoefficient(
        const DeviceData& data, Set::Scalar temperature,
        Set::Scalar pressure, const Composition& partial_density,
        int i, int j, int k, int species)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->diffusion_coefficient(
            temperature, pressure, X, i, j, k, species);
#else
        const int count = amrex::get<0>(data);
        if (count == 1) return 0.0;
        const int transport_type = amrex::get<7>(data);
        const Set::Scalar gas_density =
            GasDensity(data, partial_density, i, j, k);
        if (transport_type == 2)
        {
            const auto& molecular_weight = amrex::get<2>(data);
            Set::Scalar moles = 0.0;
            for (int n = 0; n < count; ++n)
                moles += partial_density(i,j,k,n) /
                    molecular_weight[n];
            const Set::Scalar mixture_molecular_weight =
                moles > 0.0 ? gas_density / moles :
                molecular_weight[0];
            const Set::Scalar intrinsic_density =
                pressure > 0.0 && temperature > 0.0 ?
                pressure * mixture_molecular_weight /
                    (amrex::get<1>(data) * temperature) : 0.0;
            const Set::Scalar cp =
                CpMass(data, temperature, partial_density, i, j, k);
            return intrinsic_density > 0.0 && cp > 0.0 ?
                ThermalConductivity(
                    data, temperature, partial_density, i, j, k) /
                    (intrinsic_density * cp * amrex::get<14>(data)) :
                0.0;
        }

        const auto& molecular_weight = amrex::get<2>(data);
        const auto& value1 = amrex::get<8>(data);
        const auto& value2 = amrex::get<9>(data);
        Set::Scalar moles = 0.0;
        for (int n = 0; n < count; ++n)
            moles += partial_density(i,j,k,n) / molecular_weight[n];
        auto mole_fraction = [=] AMREX_GPU_HOST_DEVICE(int n)
        {
            return moles > 0.0 ?
                partial_density(i,j,k,n) /
                    molecular_weight[n] / moles :
                (n == 0 ? 1.0 : 0.0);
        };
        if (transport_type == 0)
        {
            const Set::Scalar mixture_molecular_weight =
                moles > 0.0 ? gas_density / moles :
                molecular_weight[0];
            const Set::Scalar intrinsic_density =
                pressure > 0.0 && temperature > 0.0 ?
                pressure * mixture_molecular_weight /
                    (amrex::get<1>(data) * temperature) : 0.0;
            if (!(intrinsic_density > 0.0)) return 0.0;
            const Set::Scalar self_a =
                amrex::get<10>(data) * value1[species] /
                    intrinsic_density;
            Set::Scalar inverse_diffusion = 0.0;
            for (int b = 0; b < count; ++b)
            {
                if (species == b) continue;
                const Set::Scalar self_b =
                    amrex::get<10>(data) * value1[b] /
                        intrinsic_density;
                const Set::Scalar binary =
                    2.0 * self_a * self_b / (self_a + self_b);
                inverse_diffusion += mole_fraction(b) / binary;
            }
            return inverse_diffusion > 0.0 ?
                (1.0 - mole_fraction(species)) / inverse_diffusion : 0.0;
        }

        Set::Scalar inverse_diffusion = 0.0;
        for (int b = 0; b < count; ++b)
        {
            if (species == b) continue;
            const Set::Scalar sigma =
                0.5 * (value1[species] + value1[b]) / amrex::get<12>(data);
            const Set::Scalar epsilon =
                std::sqrt(value2[species] * value2[b]);
            const Set::Scalar Tstar = temperature / epsilon;
            const Set::Scalar collision =
                1.06036 / std::pow(Tstar, 0.15610) +
                0.19300 / std::exp(0.47635 * Tstar) +
                1.03587 / std::exp(1.52996 * Tstar) +
                1.76474 / std::exp(3.89411 * Tstar);
            Set::Scalar binary =
                0.0018583 *
                std::sqrt(temperature * temperature * temperature *
                          (1.0 / molecular_weight[species] +
                           1.0 / molecular_weight[b])) /
                ((pressure / amrex::get<11>(data)) *
                 sigma * sigma * collision);
            binary *= amrex::get<13>(data);
            inverse_diffusion += mole_fraction(b) / binary;
        }
        return inverse_diffusion > 0.0 ?
            (1.0 - mole_fraction(species)) / inverse_diffusion : 0.0;
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar ComputeEnergy(
        const DeviceData& data, Set::Scalar density,
        Set::Scalar momentumx, Set::Scalar momentumy,
        Set::Scalar temperature, const Composition& partial_density,
        int i, int j, int k)
    {
#ifndef ALAMO_GPU
        const Gas* gas = amrex::get<18>(data);
        const MoleFraction X =
            gas->MoleFractions(partial_density, i, j, k);
        return gas->ComputeE(
            density, momentumx, momentumy, temperature, X, i, j, k);
#else
        const Set::Scalar gas_constant =
            GasConstant(data, partial_density, i, j, k);
        const Set::Scalar cp =
            CpMass(data, temperature, partial_density, i, j, k);
        const Set::Scalar reference_temperature =
            amrex::get<6>(data) == 1 ? 298.15 : 0.0;
        const Set::Scalar internal_energy =
            (cp - gas_constant) *
            (temperature - reference_temperature);
        return density * internal_energy +
            0.5 * (momentumx*momentumx + momentumy*momentumy) / density;
#endif
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double cp_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.cp_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double enthalpy_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.enthalpy_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double entropy_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.entropy_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double dynamic_viscosity(double T, const Composition& X, int i, int j, int k) const
    {
        return transport.dynamic_viscosity(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double thermal_conductivity(double T, const Composition& X, int i, int j, int k) const
    {
        return transport.thermal_conductivity(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double diffusion_coefficient(double T, double P, const Composition& X,
                                int i, int j, int k, int species) const
    {
        return transport.diffusion_coefficient(T, P, X, i, j, k, species);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    const Composition& X, int i, int j, int k) const
    {
        return eos.ComputeE(*this, density, momentumx, momentumy, T, X, i, j, k);
    }

// Universal methods inherent to Gas
public:
    void ComputeLocalFractions(
            Set::Patch<const Set::Scalar> density_mf,
            Set::Patch<Set::Scalar> mass_fraction_mf,
            Set::Patch<Set::Scalar> mole_fraction_mf,
            const int i, const int j, const int k)
    {
        // In place update for mass and mole fractions
        // This needs to be called before any other functions if cell composition 
        // has changed, since they are dependent on mole fraction
        double density = 0.0;
        double moles = 0.0;
        for (int n=0; n<nspecies; ++n) {
            density += density_mf(i,j,k,n);
            moles += density_mf(i,j,k,n) / MW[n];
        }
        for (int n=0; n<nspecies; ++n) {
            mass_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / density;
            mole_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / MW[n] / moles;
        }
    }
    template<class Composition>
    double GetMW(const Composition& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<nspecies; ++n) mw += X(i,j,k,n) * MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar> rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<nspecies; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    template<class Composition>
    double R(const Composition& X, int i, int j, int k) const {
        // Specific gas constant, J/(kg-K)
        return Rg/GetMW(X, i, j, k);
    }
    template<class Composition>
    double cp_mass(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double cv_mol(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Rg;
    }
    template<class Composition>
    double cv_mass(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double gamma(double T, const Composition& X, int i, int j, int k) const  {
        // Specific heat ratio
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    double enthalpy_mass(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double enthalpy_mass(double T, const Composition& X, int i, int j, int k) const {
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    double entropy_mass(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif
