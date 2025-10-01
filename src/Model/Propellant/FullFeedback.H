//
// This SCP model implements the method described in the following paper
//
// .. bibliography::
//    :list: none
//    :filter: False
//
//    meier2024diffuse
//
// This method imitates inheritance :code:`Propellant::Propellant` using static polymorphism
//

#ifndef MODEL_PROPELLANT_FULLFEEDBACK_H
#define MODEL_PROPELLANT_FULLFEEDBACK_H

#include "AMReX_Extension.H"
#include "AMReX_GpuQualifiers.H"
#include "IO/ParmParse.H"
#include "Model/Propellant/Propellant.H"
#include "Set/Base.H"
#include "Unit/Unit.H"

#include <cmath>
#include <cstddef>
#include <string>

namespace Model
{
namespace Propellant
{
class FullFeedback
{
public:
    static constexpr const char *name = "fullfeedback";
    static constexpr std::size_t n_species = 2;

    FullFeedback() = default;
    FullFeedback(IO::ParmParse &pp, std::string name)
    {
        pp_queryclass(name, *this);
    }

    void set_pressure(const Set::Scalar a_P)
    {
        this->Pbar = a_P / Pref;

        zeta_2 = zeta_fit_a + Pbar * zeta_fit_b;

        if (arrhenius_dependency == 1)
            zeta_1 = zeta_2;
        else
            zeta_1 = zeta_0;

        // Pressure-dependent post calculation variables
        k1 = a1 * Pbar + b1 - zeta_1 / zeta;
        k2 = a2 * Pbar + b2 - zeta_1 / zeta;
        k3 = 4.0 * std::log((c1 * Pbar * Pbar + a3 * Pbar + b3) - k1 / 2.0 - k2 / 2.0);

        return;
    }

    AMREX_FORCE_INLINE
    Set::Scalar get_qdot(const Set::Scalar mdot, const Model::Propellant::PhiArray<n_species> &phi)
    {
        Set::Scalar qflux = k1 * phi[0] + k2 * (1.0 - phi[0]) + (zeta_1 / zeta) * std::exp(k3 * phi[0] * (1.0 - phi[0]));
        Set::Scalar mlocal = (mlocal_ap)*phi[0] + (mlocal_htpb) * (1.0 - phi[0]) + mlocal_comb * phi[0] * (1.0 - phi[0]);
        Set::Scalar mdota = std::fabs(mdot);
        Set::Scalar mbase = std::tanh(4.0 * mdota / (mlocal));

        return mbase * qflux;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_K(const Model::Propellant::PhiArray<n_species> &phi)
    {
        return k_ap * phi[0] + k_htpb * (1.0 - phi[0]);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_rho(const Model::Propellant::PhiArray<n_species> &phi)
    {
        return rho_ap * phi[0] + rho_htpb * (1.0 - phi[0]);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_cp(const Model::Propellant::PhiArray<n_species> &phi)
    {
        return cp_ap * phi[0] + cp_htpb * (1.0 - phi[0]);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_L(const Model::Propellant::PhiArray<n_species> &phi, Set::Scalar T)
    {
        Set::Scalar L = NAN;
        if (mob_ap == 1)
            L = m_ap * Pbar * std::exp(-E_ap / T) * phi[0];
        else
            L = m_ap * std::exp(-E_ap / T) * phi[0];
        L += m_htpb * std::exp(-E_htpb / T) * (1.0 - phi[0]);
        if (T <= bound)
            L = 0.0;
        return L;
    }

    static void Parse(FullFeedback &value, IO::ParmParse &pp)
    {
        // AP/HTPB interface length
        pp.query_default("phi.zeta", value.zeta, "10.0_um", Unit::Length());
        // Reference interface length for heat integration
        pp.query_default("phi.zeta_0", value.zeta_0, "10.0_um", Unit::Length());

        // Constant offset for $\zeta$ adjustment with pressure
        pp.query_default("phi.zeta_fit_a", value.zeta_fit_a, "45.0_um", Unit::Length());
        // Scaling factors for $\zeta$ adjustment with pressure
        pp.query_default("phi.zeta_fit_b", value.zeta_fit_b, "-6.42_um", Unit::Length());

        // Surrogate heat flux model paramater - AP
        pp_query_required("a1", value.a1, Unit::Less());
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("a2", value.a2, Unit::Less());
        // Surrogate heat flux model paramater - Total
        pp_query_required("a3", value.a3, Unit::Less());
        // Surrogate heat flux model paramater - AP
        pp_query_required("b1", value.b1, Unit::Less());
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("b2", value.b2, Unit::Less());
        // Surrogate heat flux model paramater - Total
        pp_query_required("b3", value.b3, Unit::Less());
        // Surrogate heat flux model paramater - Total
        pp_query_required("c1", value.c1, Unit::Less());

        // Whether to use pressure to determined the reference Zeta
        pp_query_default("pressure.dependency", value.arrhenius_dependency, true);

        // Reference pressure for pressure nondimensionalization
        pp_query_default("Pref", value.Pref, "1_MPa", Unit::Pressure());

        // AP volumetric mass flux reference value
        pp.query_default("mlocal_ap", value.mlocal_ap, "0.0", Unit::Density() / Unit::Time());
        // combined volumetric mass flux reference value
        pp_query_default("mlocal_comb", value.mlocal_comb, "0.0", Unit::Density() / Unit::Time());
        // HTPB volumetric mass flux reference value
        pp.query_default("mlocal_htpb", value.mlocal_htpb, "5000_kg/m^3/s", Unit::Density() / Unit::Time());

        // Set::Scalar massfraction = 0.8;
        // value.mlocal_htpb = 685000.0 - 850e3 * massfraction; // Replicating but needs to be revised!
        // Util::Message(INFO,value.mlocal_htpb);

        // AP Density
        pp_query_required("rho_ap", value.rho_ap, Unit::Density());
        // HTPB Density
        pp_query_required("rho_htpb", value.rho_htpb, Unit::Density());
        // AP Thermal Conductivity
        pp_query_required("k_ap", value.k_ap, Unit::ThermalConductivity());
        // HTPB Thermal Conductivity
        pp_query_required("k_htpb", value.k_htpb, Unit::ThermalConductivity());
        // AP Specific Heat
        pp_query_required("cp_ap", value.cp_ap, Unit::SpecificHeatCapacity());
        // HTPB Specific Heat
        pp_query_required("cp_htpb", value.cp_htpb, Unit::SpecificHeatCapacity());

        // AP Pre-exponential factor for Arrhenius Law
        pp_query_required("m_ap", value.m_ap, Unit::Velocity());
        // HTPB Pre-exponential factor for Arrhenius Law
        pp_query_required("m_htpb", value.m_htpb, Unit::Velocity());
        // COMBINED AP Activation Energy, divided by Rg, for Arrhenius Law (has units of temperature)
        pp_query_required("E_ap", value.E_ap, Unit::Temperature());
        // COMBINED HTPB Activation Energy, divided by Rg, for Arrhenius Law(has units of temperature)
        pp_query_required("E_htpb", value.E_htpb, Unit::Temperature());
        // Whether to include pressure to the arrhenius law [??]n
        pp_query_default("mob_ap", value.mob_ap, false);
        // HTPB Activation Energy for Arrhenius Law
        pp_query_default("bound", value.bound, "0.0", Unit::Temperature());
    }

    // IO Variables
    Set::Scalar zeta = NAN, zeta_0 = NAN;
    Set::Scalar a1 = NAN, a2 = NAN, a3 = NAN;
    Set::Scalar b1 = NAN, b2 = NAN, b3 = NAN;
    Set::Scalar c1 = NAN;
    // Set::Scalar h1 = NAN, h2 = NAN;
    bool arrhenius_dependency = false;
    Set::Scalar mlocal_ap = NAN, mlocal_comb = NAN, mlocal_htpb = NAN;

    // Thermal properties
    Set::Scalar k_ap = NAN, k_htpb = NAN;
    Set::Scalar rho_ap = NAN, rho_htpb = NAN;
    Set::Scalar cp_ap = NAN, cp_htpb = NAN;

    // Regression variables
    int mob_ap = -1;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;
    Set::Scalar Pref = NAN;

    // Calculated variables
    Set::Scalar Pbar = NAN;
    Set::Scalar zeta_1 = NAN, zeta_2 = NAN;
    Set::Scalar zeta_fit_a = NAN, zeta_fit_b = NAN;
    Set::Scalar k1 = NAN, k2 = NAN, k3 = NAN, k4 = NAN;
};

}
}

#endif
