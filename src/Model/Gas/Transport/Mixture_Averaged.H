#ifndef MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_
#define MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_

#include <vector>
#include <memory>
#include "IO/ParmParse.H"
#include "Model/Gas/Transport/Transport.H"

// Type: 0,1 determines if per species viscosity and thermal conductivity are 
//           constant (type=constant) or computed using kinetic theory (type=LJ)
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

class Mixture_Averaged : public Transport {
public:
    static constexpr const char* name = "mixture_averaged";
    virtual const char* model_name() const override { return name; }

private:
    int nspecies;
    std::vector<double> &MW;
    Thermo::Thermo const * const thermo;
    std::string type;
    std::vector<double> val1; // constant viscosity mu or Lennard-Jones collision diameter depending on type
    std::vector<double> val2; // constant thermal conductivity or Lennard-Jones potential well depth depending on type
    double d_fact; // When type=constant, this is a correction factor such that self diffusion D_n = d_fact * mu_n/rho_mix
    Set::Scalar P_1atm, L_1ang, A_1cmsq;

public:
    Mixture_Averaged() = delete;
    Mixture_Averaged(int a_nspecies, std::vector<double> &a_MW, Thermo::Thermo *a_thermo, IO::ParmParse& pp, std::string name)
        : nspecies(a_nspecies), MW(a_MW), thermo(a_thermo) // thermo is owned by Gas, should only be referenced here
    {
        pp.queryclass(name, *this);
    }

    ~Mixture_Averaged() override = default;

    static void Parse(Mixture_Averaged & value, IO::ParmParse & pp)
    {
        pp.query_validate("type", value.type, {"constant", "LJ"});
        if (value.type == "constant")
        {
            pp.queryarr_required("mu", value.val1, Unit::Pressure() * Unit::Time());
            pp.queryarr_required("k", value.val2, Unit::Power() / Unit::Length() / Unit::Temperature());
            pp.query_default("d_fact", value.d_fact, "1.1", Unit::Less());
        }
        else if (value.type == "LJ")
        {
            pp.queryarr_required("LJdiameter", value.val1, Unit::Length());
            pp.queryarr_required("LJwelldepth", value.val2, Unit::Temperature());
        }
        else
        {
            Util::Abort(INFO, "Unrecognized type for transport model ", value.type);
        }
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
        return 2.669570e-6 * sqrt(MW[i] * T) / (val1[i]/L_1ang * val1[i]/L_1ang * collision_integral(T / val2[i]));
    }
    double k_eucken(double T, int n) const 
    {
        if (!thermo) Util::Abort(INFO, "[Transport::Mixture_Averaged::k_eucken] Thermo model not set.");
        double cp = thermo->cp_mol_species(T, n);
        return mu_enskog(T, n) * (cp + 1.25 * Set::Constant::Rg) / MW[n];
    }

public:
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override 
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

    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override 
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

    void diffusion_coeffs(
            Set::Patch<Set::Scalar>& DKM, double T, double P,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) override
    {
        if ( type == "constant" )
        {
            // Mixture diffusion coefficients based on simple harmonic mean, m^2/s
            double rho_mix = 0.0;
            for (int n=0; n<nspecies; ++n)
                rho_mix += X(i,j,k,n) * MW[n];
            rho_mix *= P/(Set::Constant::Rg*T);
            
            // compute self-diffusion proxies D_n = mu_n/rho_mix
            std::vector<double> D_self(nspecies);
            for (int n=0; n<nspecies; ++n)
                D_self[n] = d_fact * val1[n] / rho_mix;
            
            // compute mixture diffusion coefficients
            for (int a=0; a<nspecies; ++a) {
                double sumD = 0.0;
                for (int b=0; b<nspecies; ++b) {
                    if (a==b) continue;
                    double D_ab = 2.0 * D_self[a] * D_self[b] / (D_self[a] + D_self[b]); // harmonic mean
                    sumD += X(i,j,k,b) / D_ab;
                }
                DKM(i, j, k, a) = (1.0 - X(i, j, k ,a)) / sumD;
                if (!(DKM(i, j, k, a) == DKM(i, j, k, a))) DKM(i, j, k, a) = 0.0; // NaN guard
            }
        }
        else if ( type == "LJ" )
        {
            // Mixture diffusion coefficients (Chapman–Enskog), m^2/s
            for (int a = 0; a < nspecies; ++a) {
                double sumD = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    if (a == b) continue;

                    double sigmaAB = 0.5 * (val1[a] + val1[b]) / L_1ang; // in angstroms
                    double epsAB   = sqrt(val2[a] * val2[b]);
                    double Tstar   = T / epsAB;

                    double omegaAB =
                        1.06036 / pow(Tstar, 0.15610) +
                        0.19300 / exp(0.47635 * Tstar) +
                        1.03587 / exp(1.52996 * Tstar) +
                        1.76474 / exp(3.89411 * Tstar);

                    double DAB = 
                        0.0018583 * sqrt(T*T*T * (1.0 / MW[a] + 1.0 / MW[b])) /
                        ( (P / 101325.0) * sigmaAB*sigmaAB * omegaAB );
                    DAB *= A_1cmsq; // convert cm^2/s to m^2/s

                    sumD += X(i, j, k, b) / DAB;
                }

                DKM(i, j, k, a) = (1.0 - X(i, j, k ,a)) / sumD;
                if (!(DKM(i, j, k, a) == DKM(i, j, k, a))) DKM(i, j, k, a) = 0.0; // NaN guard
            }
        }
    }

}; // class Mixture_Averaged

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif
