//
//
//
//
//

#ifndef MODEL_PROPELLANT_HOMOGENIZE_H
#define MODEL_PROPELLANT_HOMOGENIZE_H

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
class Homogenize
{
public:
    static constexpr const char *name = "homogenize";
    static constexpr std::size_t n_species = 2;

    Homogenize() = default;
    Homogenize(IO::ParmParse &pp, std::string name)
    {
        pp_queryclass(name, *this);
    }

    Set::Scalar P = NAN;
    void set_pressure(const Set::Scalar P)
    {
        k4 = h1 * P + h2;
    }

    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Model::Propellant::PhiArray<n_species> &phi)
    {
        Set::Scalar qflux = k4 * phi[0];
        Set::Scalar mlocal = (mlocal_ap)*massfraction + (mlocal_htpb) * (1.0 - massfraction);
        Set::Scalar mdota = std::fabs(mdot);
        Set::Scalar mbase = std::tanh(4.0 * mdota / (mlocal));
        return mbase * qflux;
    }

    Set::Scalar get_K(const Model::Propellant::PhiArray<n_species> &phi)
    {
        return (k_ap * massfraction + k_htpb * (1.0 - massfraction)) * phi[0] + dispersion1 * (1. - phi[0]);
    }

    Set::Scalar get_rho(const Model::Propellant::PhiArray<n_species> &phi)
    {

        return (rho_ap * massfraction + rho_htpb * (1.0 - massfraction)) * phi[0] + dispersion2 * (1. - phi[0]);
    }

    Set::Scalar get_cp(const Model::Propellant::PhiArray<n_species> &phi)
    {
        return (cp_ap * massfraction + cp_htpb * (1.0 - massfraction)) * phi[0] + dispersion3 * (1. - phi[0]);
    }

    Set::Scalar get_L(const Model::Propellant::PhiArray<n_species> &phi, Set::Scalar T)
    {
        Set::Scalar L = NAN;
        if (mob_ap == 1)
            L = m_ap * P * std::exp(-E_ap / T) * phi[0];
        else
            L = m_ap * std::exp(-E_ap / T) * phi[0];
        L += m_htpb * std::exp(-E_htpb / T) * (1.0 - phi[0]);
        if (T <= bound)
            L = 0.0;
        return L;
    }

    static void Parse(Homogenize &value, IO::ParmParse &pp)
    {

        // K; dispersion variables are use to create an inert region for the void grain case.
        // An inert region is one that dissipates energy fast enough to remove regression and
        // thermal strain effects.
        pp_query("dispersion1", value.dispersion1);
        // rho; dispersion variables are use to create an inert region for the void grain case.
        // An inert region is one that dissipates energy fast enough to remove regression and
        // thermal strain effects.
        pp_query("dispersion2", value.dispersion2);
        // cp; dispersion variables are use to create an inert region for the void grain case.
        // An inert region is one that dissipates energy fast enough to remove regression and
        // thermal strain effects.
        pp_query("dispersion3", value.dispersion3);

        // AP Density
        pp_query_required("rho_ap", value.rho_ap);
        // HTPB Density
        pp_query_required("rho_htpb", value.rho_htpb);
        // AP Thermal Conductivity
        pp_query_required("k_ap", value.k_ap);
        // HTPB Thermal Conductivity
        pp_query_required("k_htpb", value.k_htpb);
        // AP Specific Heat
        pp_query_required("cp_ap", value.cp_ap);
        // HTPB Specific Heat
        pp_query_required("cp_htpb", value.cp_htpb);

        // Surgate heat flux model paramater - Homogenized
        pp_query_default("h1", value.h1, 1.81);
        // Surgate heat flux model paramater - Homogenized
        pp_query_default("h2", value.h2, 1.34);

        // AP/HTPB ratio for homogenized domain
        pp_query_default("massfraction", value.massfraction, 0.8);

        // AP mass flux reference value
        pp_query_default("mlocal_ap", value.mlocal_ap, 0.0);
        value.mlocal_htpb = 685000.0 - 850e3 * value.massfraction;

        // AP Pre-exponential factor for Arrhenius Law
        pp_query_required("m_ap", value.m_ap, Unit::Velocity());
        // HTPB Pre-exponential factor for Arrhenius Law
        pp_query_required("m_htpb", value.m_htpb, Unit::Velocity());
        // AP Activation Energy for Arrhenius Law
        pp_query_required("E_ap", value.E_ap, Unit::Energy());
        // HTPB Activation Energy for Arrhenius Law
        pp_query_required("E_htpb", value.E_htpb, Unit::Energy());
        // Whether to include pressure to the arrhenius law [??]
        pp_query_default("mob_ap", value.mob_ap, false);
        // HTPB Activation Energy for Arrhenius Law
        pp_query_default("bound", value.bound, 0.0);
    }

    // IO variables
    Set::Scalar massfraction = NAN;
    Set::Scalar dispersion1 = NAN, dispersion2 = NAN, dispersion3 = NAN;
    Set::Scalar k_ap = NAN, k_htpb = NAN;
    Set::Scalar rho_ap = NAN, rho_htpb = NAN;
    Set::Scalar cp_ap = NAN, cp_htpb = NAN;
    Set::Scalar h1 = NAN, h2 = NAN;
    Set::Scalar mlocal_ap = NAN;
    // Regression variables
    bool mob_ap = -1;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;

    // Calculated variables
    Set::Scalar k4 = NAN;
    Set::Scalar mlocal_htpb = NAN;
};

}
}

#endif
