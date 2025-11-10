#ifndef MODEL_GAS_MIXTURE_H
#define MODEL_GAS_MIXTURE_H

// This is where we compute mixture dependent values for thermodynamics
// and transport.  Species specific data is stored in Model::Gas::Species

namespace Model {
namespace Gas {

// If Mixture is constructed without passing a transport model,
// then no assume no transport (Euler equations)
struct Mixture {
    // Mixed thermo values
    double density;
    double cp;
    double cv;
    double gamma;
    double MW;
    double R;
    std::vector<double> X;
    std::vector<double> Y;

    // Constructor
    explicit Mixture() { density=0.0, cp=0.0, cv=0.0, gamma=0.0, MW=0.0, R=0.0; }

    void ComputeThermo(const std::vector<double>& rhoY,
                       const Gas::Species& species,
                       const double T) {
        if (T != species.T) species.UpdateTemperature(T);
        int nspecies = rhoY.size();
        X = std::vector<double>(nspecies, 0.0);
        Y = std::vector<double>(nspecies, 0.0);
        density = 0.0;
        for (int n=0; n<nspecies; ++n) density += rhoY[n];
        double moles = 0.0;
        for (int n=0; n<nspecies; ++n) {
            Y[n] = rhoY[n]/density;
            moles += Y[n]/species.MW[n];
        }
        cp = 0.0;
        MW = 0.0;
        for (int n=0; n<nspecies; +=n) {
            X[n] = Y[n] / species.MW[n] / moles;
            cp += Y[n] * species.cp[n];
            MW += X[n] * species.MW[n];
        }
        R = Ru/MW;
        cv = cp - R;
        gamma = cp/cv;
    }
        
}

// If Mixture is constructed with a transport model specified:
// Use these function to compute mixture values of viscosity, thermal conductivity
// and species diffusion.  Species specific values are stored in Model::Gas::Species
struct Mixture {
    // Mixed thermo values
    double density;
    double cp;
    double cv;
    double gamma;
    double MW;
    double R;
    std::vector<double> X;
    std::vector<double> Y;

    // Mixed transport values
    const std::string model;
    double thermal_conductivity;
    double dynamic_viscosity;
    std::vector<double> diffusion_coefficients;

    // Constructor
    explicit Transport(const std::string& modelName) : model(modelName)  {
        if (model != "mixture-averaged" && 
            //model != "multi-component"  &&
            model != "rocfire") {
            Util::Exception(INFO,"Unknown transport model: ", model);
        }
        density=0.0, cp=0.0, cv=0.0, gamma=0.0, MW=0.0, R=0.0;
    }

    void ComputeThermo(const std::vector<double>& rhoY,
                       const Gas::Species& species,
                       const double T) {
        if (T != species.T) species.UpdateTemperature(T);
        int nspecies = rhoY.size();
        X = std::vector<double>(nspecies, 0.0);
        Y = std::vector<double>(nspecies, 0.0);
        density = 0.0;
        for (int n=0; n<nspecies; ++n) density += rhoY[n];
        double moles = 0.0;
        for (int n=0; n<nspecies; ++n) {
            Y[n] = rhoY[n]/density;
            moles += Y[n]/species.MW[n];
        }
        cp = 0.0;
        MW = 0.0;
        for (int n=0; n<nspecies; +=n) {
            X[n] = Y[n] / species.MW[n] / moles;
            cp += Y[n] * species.cp[n];
            MW += X[n] * species.MW[n];
        }
        R = Ru/MW;
        cv = cp - R;
        gamma = cp/cv;
    }
    
    void ComputeTransport(const std::vector<double>& rhoY,
                          const Gas::Species& species,
                          const double T,
                          const double P) {
        if (T != species.T) species.UpdateTemperature(T);
        if (model == "mixture-averaged") {
            // Multispecies effects
            // Mixture viscosity and heat conduction coefficient
            // Multicomponent extension of Chapman-Enskog
            // Transport Phenomena, Revised 2nd Edition / Bird, Stewart, Lightfoot // Ch. 1.4 & 9.3
            //
            // Binary diffusion constants
            // Chapman-Enskog
            // Transport PHenomena, Revised 2nd Edition / Bird, Stewart, Lightfoot // Ch. 17.3
            //
            // Upper triangular matrix for species row and species column since D_AB = D_BA, i.e.
            //
            //                species A | species B | species C
            //              ------------------------------------
            //    species A |    D_AA   |   D_AB    |   D_AC   |
            //              ------------------------------------
            //    species B |           |   D_BB    |   D_BC   |
            //              ------------------------------------
            //    species C |           |           |   D_CC   |
            //              ------------------------------------

            double collision_integral( double nondimT ) const {
                double I = 1.06036/pow(nondimT, 0.15610) + 0.19300/exp(0.47635*nondimT) + 
                           1.03587/exp(1.52996*nondimT) + 1.76474/exp(3.89411*nondimT);
                return I;
            }

            int nspecies = rhoY.size();
            thermal_conductivity = 0.0;
            dynamic_viscosity = 0.0;

            double sigmaAB = 0.0;
            double epsAB = 0.0;
            double nondimT = 0.0;
            double omegaAB = 0.0;
            double DAB = 0.0;
            for (int a=0; a<nspecies; ++a)
            {
                double phi = 0.0;
                diffusion_coefficient[a] = 0.0;
                for (int b=0; b<nspecies; ++b)
                {
                    phi += X[b] * 1.0/sqrt(8.0) *
                           pow(1.0 + species.MW[a]/species.MW[b], -0.5) *
                           pow(1.0 + sqrt(species.dynamic_viscosity[a]/species.dynamic_viscosity[b]) *
                           pow(species.MW[b]/species.MW[a], 0.25), 2.0);
                    if ( b != a )
                    {
                        sigmaAB = 0.5*(species.LJdiameter[a] + species.LJdiameter[b]);
                        epsAB = sqrt(species.LJwelldepth[a] * species.LJwelldepth[b]);
                        nondimT = T/epsAB;
                        omegaAB = collision_integral(nondimT);
                        DAB = 0.0018583*sqrt(pow(T, 3.0) * (1.0/species.MW[a] + 1.0/species.MW[b])) /
                                    (P/101325.0*pow(sigmaAB,2.0)*omegaAB);
                        DAB /= 10000.0; // convert from cm^2/s to m^2/s
                        diffusion_coefficient[a] += X[b]/DAB;
                    }
                }
                diffusion_coefficient[a] = (1.0 - X[a])/diffusion_coefficicent[a];

                // Set to zero if nan or inf (pure species locally)
                if ( diffusion_coefficient[a] != diffusion_coefficient[a] ) {
                    diffusion_coefficient[a] = 0.0;
                }

                dynamic_viscosity += X[a] * species.dynamic_viscosity[a] / phi;
                thermal_conductivity  += X[a] * species.thermal_conductivity[a] / phi;
            }
        } // mixture-averaged model
    } // ComputeTransport
}; // Mixture struct

} // namespace Gas
} // namespace Model

#endif
