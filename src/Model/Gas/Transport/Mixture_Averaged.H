#ifndef MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_
#define MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_

#include <vector>
#include <array>
#include <memory>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Thermo/GrossModel.H"
#include "Model/Gas/Transport/Transport.H"

// Type: 0,1 determines if per species viscosity and thermal conductivity are 
//           constant (type=0) or computed using kinetic theory (type=1)
// Either way, mixture averaging is used to compute the mixture values
//
// For kinetic theory: Chapman-Enskog is used to compute viscosity
//                     and Eucken relation used to compute conductivity


// dynamic viscosity, mu, computed using Chapman-Enskog theory
//
//        5E20 sqrt(pi * m * kb * T)  
//  mu = ---------------------------- 
//              16 pi σ^2 Ω           
//
//            MW
//  where m = -- (mass per particle), MW: molecular weight, NA: Avogadro's number, kb: Boltzman constant
//            NA                                                                  
//
//  σ: Lennard-Jones collision diameter in Angstroms, Ω: collision integral


// thermal conductivity coefficient, k, computed using Eucken relation
//
//  k = mu * (cp + 1.25R), where cp and R are both per kg
//

namespace Model {
namespace Gas {
namespace Transport {

class Mixture_Averaged // public Transport
{
public:
    static constexpr const char* name = "mixture_averaged";
    static constexpr bool common_diffusivity = false;

private:
    int nspecies;
    double Rg;
    std::vector<double> MW;
    Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::GrossModel> thermo;
    std::string type;
    std::vector<double> val1; // constant viscosity mu or Lennard-Jones collision diameter depending on type
    std::vector<double> val2; // constant thermal conductivity or Lennard-Jones potential well depth depending on type
    Set::Scalar d_fact;
    Set::Scalar P_1atm, L_1ang, A_1cmsq;

public:
    Mixture_Averaged()
    {
        Rg = Set::Constant::Rg;
    }

    static void Parse(  Mixture_Averaged & value, IO::ParmParse & pp,
                        int a_nspecies, std::vector<double> &a_MW,
                        Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::GrossModel> a_thermo)
    {
        value.nspecies = a_nspecies;
        value.MW = a_MW;
        value.thermo = a_thermo;

        // Compute individual species transport properties as "constant" or from "LJ" (Lennard-Jones) theory
        pp.query_switch("type",{
                {"constant", [&]() {
                    value.type = "constant";
                    // Dynamic viscosity if type=constant
                    pp.queryarr_required("mu", value.val1, Unit::Pressure() * Unit::Time());
                    // Thermal conductivity if type=constant
                    pp.queryarr_required("k", value.val2, Unit::Power() / Unit::Length() / Unit::Temperature());
                    pp.query_default("d_fact", value.d_fact, "1.1", Unit::Less());
                }},
                {"LJ", [&]() {
                    value.type = "LJ";
                    // Lennard-Jones diameter if type=LJ
                    pp.queryarr_required("LJdiameter", value.val1, Unit::Length());
                    // Lennard-Jones well depth if type=LJ
                    pp.queryarr_required("LJwelldepth", value.val2, Unit::Temperature());
                }}
            });
        pp.query_default("P_1atm", value.P_1atm, "1_atm", Unit::Pressure());
        pp.query_default("L_1ang", value.L_1ang, "1_ang", Unit::Length());
        pp.query_default("A_1cmsq", value.A_1cmsq, "1_cm^2", Unit::Area());

        Util::Assert(INFO, TEST(value.val1.size() == (size_t)value.nspecies));
        Util::Assert(INFO, TEST(value.val2.size() == (size_t)value.nspecies));
    }

private:
    double collision_integral(double Tstar) const {
        double S =  1.16145 * pow(Tstar, -0.14874) 
                    + 0.52487 / exp(0.77320 * Tstar)
                    + 2.16178 / exp(2.43787 * Tstar);
        return S;
    }
    double mu_enskog(double T, int i) const {
        const double sigma = val1[i] / L_1ang;
        return 2.669570e-6 * sqrt(MW[i] * T) / (sigma * sigma * collision_integral(T / val2[i]));
    }
    double k_eucken(double T, int n) const
    {
        return mu_enskog(T, n) * (thermo.cp_mol_species(T, n) + 1.25 * Rg) / MW[n];
    }

public:
    auto DeviceData() const
    {
        return std::tie(type, val1, val2, d_fact, P_1atm, L_1ang, A_1cmsq);
    }

    template<class Composition>
    double dynamic_viscosity(double T, const Composition& X, int i, int j, int k) const
    {
        // Wilke's Average: Dynamic viscosity, Pa-s
        double mu_mix = 0.0;
        if ( type == "constant" )
        {
            // dynamic viscosity, mu, is constant per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = val1[a];
                double MWa = MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = val1[b];
                    double MWb = MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i,j,k,b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                mu_mix += X(i, j, k, a) * mua / phi;
            }
        }
        else if ( type == "LJ" )
        {
            // dynamic viscosity, mu, computed by Chapman-Enskog per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = mu_enskog(T, a);
                double MWa = MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = mu_enskog(T, b);
                    double MWb = MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                mu_mix += X(i, j, k, a) * mua / phi;
            }
        }
        else
        {
            Util::Abort(INFO, "[Transport::Mixture_Averaged] Unknown type.");
        }

        if (!(mu_mix == mu_mix) ) mu_mix = 0.0;
        return mu_mix;
    }

    template<class Composition>
    double thermal_conductivity(double T, const Composition& X, int i, int j, int k) const
    {
        // Wilke's Average: Thermal conductivity coefficient, W/(m-K)
        double k_mix = 0.0;
        if ( type == "constant" )
        {
            // thermal conductivity coefficient, k, is constant per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = val1[a];
                double MWa = MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = val1[b];
                    double MWb = MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                k_mix += X(i, j, k, a) * val2[a] / phi;
            }
        }
        else if ( type == "LJ" )
        {
            // thermal conductivity coefficient, k, computed by Eucken relation per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = mu_enskog(T, a);
                double MWa = MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = mu_enskog(T, b);
                    double MWb = MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                k_mix += X(i, j, k, a) * k_eucken(T, a) / phi;
            }
        }
        else
        {
            Util::Abort(INFO, "[Transport::Mixture_Averaged] Unknown type.");
        }

        return k_mix;
    }

    template<class Composition>
    double diffusion_coefficient(double T, double P, const Composition& X,
                                int i, int j, int k, int a) const
    {
        if (nspecies == 1) return 0.0;

        if (type == "constant")
        {
            double mw_mix = 0.0;
            for (int n = 0; n < nspecies; ++n)
                mw_mix += X(i,j,k,n) * MW[n];
            const double rho_mix = P * mw_mix / (Set::Constant::Rg * T);
            const double D_self_a = d_fact * val1[a] / rho_mix;
            double sumD = 0.0;
            for (int b = 0; b < nspecies; ++b) {
                if (a == b) continue;
                const double D_self_b = d_fact * val1[b] / rho_mix;
                const double D_ab = 2.0 * D_self_a * D_self_b /
                                    (D_self_a + D_self_b);
                sumD += X(i,j,k,b) / D_ab;
            }
            return sumD > 0.0 ? (1.0 - X(i,j,k,a)) / sumD : 0.0;
        }
        if (type == "LJ")
        {
            double sumD = 0.0;
            for (int b = 0; b < nspecies; ++b) {
                if (a == b) continue;

                double sigmaAB = 0.5 * (val1[a] + val1[b]) / L_1ang;
                double epsAB   = sqrt(val2[a] * val2[b]);
                double Tstar   = T / epsAB;

                double omegaAB =
                    1.06036 / pow(Tstar, 0.15610) +
                    0.19300 / exp(0.47635 * Tstar) +
                    1.03587 / exp(1.52996 * Tstar) +
                    1.76474 / exp(3.89411 * Tstar);

                double DAB = 
                    0.0018583 * sqrt(T*T*T * (1.0 / MW[a] + 1.0 / MW[b])) /
                    ((P / P_1atm) * sigmaAB*sigmaAB * omegaAB);
                DAB *= A_1cmsq;

                sumD += X(i, j, k, b) / DAB;
            }
            return sumD > 0.0 ? (1.0 - X(i,j,k,a)) / sumD : 0.0;
        }
        Util::Abort(INFO, "[Transport::Mixture_Averaged] Unknown type.");
        return NAN;
    }

    template<class Composition>
    void diffusion_coeffs(Set::Patch<Set::Scalar> DKM, double T, double P,
                        const Composition& X, int i, int j, int k) const
    {
        for (int a = 0; a < nspecies; ++a)
            DKM(i,j,k,a) = diffusion_coefficient(T, P, X, i, j, k, a);
    }

}; // class Mixture_Averaged

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif
