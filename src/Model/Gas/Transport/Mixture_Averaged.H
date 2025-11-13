#ifndef MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_
#define MODEL_GAS_TRANSPORT_MIXTURE_AVERAGED_H_

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

class Mixture_Averaged : public virtual Gas {
public:
    int type;
    std::vector<double> val1; // constant viscosity mu or Lennard-Jones collision diameter depending on type
    std::vector<double> val2; // constant thermal conductivity or Lennard-Jones potential well depth depending on type

    Mixture_Averaged(int type_in, std::vector<double>& val1_in, std::vector<double>& val2_in)
        : Gas(0, val1) /*Dummy constructor*/, type(type_in), val1(val1_in), val2(val2_in) {};
    virtual ~Mixture_Averaged() = default;

private:
    double collision_integral(double Tstar) {
        double S =  1.16145 * pow(Tstar, -0.14874) 
                    + 0.52487 / exp(0.77320 * Tstar)
                    + 2.16178 / exp(2.43787 * Tstar);
        return S;
    }
    double mu_enskog(double T, int i) {
        return 2.669570e-6 * sqrt(_MW[i] * T) / (val1[i] * val1[i] * collision_integral(T / val2[i]));
    }
    double k_eucken(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int n) {
        return mu_enskog(T, n) * (X(i, j, k, n) * cp_mole(T, X, i, j, k) + 1.25 * Ru) / _MW[n];
    }

public:
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Wilke's Average: Dynamic viscosity, Pa-s
        switch (type) {
        case 0:
            // dynamic viscosity, mu, is constant per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = val1[a];
                double MWa = _MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = val1[b];
                    double MWb = _MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i,j,k,b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                mu_mix += X(i, j, k, a) * mua / phi;
            }
        case 1:
            // dynamic viscosity, mu, computed by Chapman-Enskog per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = mu_enskog(T, a)
                double MWa = _MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = mu_enskog(T, b);
                    double MWb = _MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                mu_mix += X(i, j, k, a) * mua / phi;
            }
        }
        return mu_mix;
    }

    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Wilke's Average: Thermal conductivity coefficient, W/(m-K)
        double k_mix = 0.0;
        switch (type) {
        case 0:
            // thermal conductivity coefficient, k, is constant per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = val1[a];
                double MWa = _MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = val1[b];
                    double MWb = _MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                k_mix += X(i, j, k, a) * val2[a] / phi;
            }
        case 1:
            // thermal conductivity coefficient, k, computed by Eucken relation per species
            for (int a = 0; a < nspecies; ++a) {
                double mua = mu_enskog(T, a)
                double MWa = _MW[a];
                double phi = 0.0;
                for (int b = 0; b < nspecies; ++b) {
                    double mub = mu_enskog(T, b);
                    double MWb = _MW[b];
                    double ratio = pow(MWb / MWa, 0.25);
                    double mu_ratio = sqrt(mua / mub);
                    double term = sqrt(1.0 + MWa / MWb);
                    phi += X(i, j, k, b) * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
                }
                k_mix += X(i, j, k, a) * k_eucken(T, X, i, j, k, a) / phi;
            }
        }
        return k_mix;
    }

    void diffusion_coeffs(
            Set::Patch<Set::Scalar>& DKM, double T, double P,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) override
    {
        // Mixture diffusion coefficients (Chapman–Enskog), m^2/s
        for (int a = 0; a < nspecies; ++a) {
            double sumD = 0.0;
            for (int b = 0; b < nspecies; ++b) {
                if (a == b) continue;

                double sigmaAB = 0.5 * (val1[a] + val1[b]);
                double epsAB   = sqrt(val2[a] * val2[b]);
                double Tstar   = T / epsAB;

                double omegaAB =
                    1.06036 / pow(Tstar, 0.15610) +
                    0.19300 / exp(0.47635 * Tstar) +
                    1.03587 / exp(1.52996 * Tstar) +
                    1.76474 / exp(3.89411 * Tstar);

                double DAB = 
                    0.0018583 * sqrt(T*T*T * (1.0 / _MW[a] + 1.0 / _MW[b])) /
                    ( (P / 101325.0) * sigmaAB*sigmaAB * omegaAB );
                DAB /= 1.0e4;

                sumD += X(i, j, k, b) / DAB;
            }

            DKM(i, j, k, a) = (1.0 - X(i, j, k ,a)) / sumD;
            if (!(DKM(i, j, k, a) == DKM(i, j, k, a))) DKM(i, j, k, a) = 0.0; // NaN guard
        }

    }

} // class Mixture_Averaged

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif
