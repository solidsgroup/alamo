#ifndef MODEL_GAS_SPECIES_H
#define MODEL_GAS_SPECIES_H

// This is where we compute mixture dependent values for thermodynamics
// and transport.  Species specific data is stored in Model::Gas::Species

namespace Model {
namespace Gas {

// Constructor when passing species, mw, cp, mu, and k
// This assumes that they are constant
struct Species {
    double T;
    const std::vector<std::string> species;
    const std::vector<double> MW;
    std::vector<double> cp;
    std::vector<double> dynamic_viscosity;
    std::vector<double> thermal_conductivity;

    // Constructor
    explicit Species( std::vector<std::string>& Asp,
                      std::vector<double>& Amw,
                      std::vector<double>& Acp,
                      std::vector<double> &Amu,
                      std::vector<double>& Ak ) 
             : species(Asp), MW(Amw), cp(Acp), dynamic_viscosity(Amu), thermal_conductivity(Ak) {}

    void UpdateTemperature(double Tnew) { T = Tnew; };
    
        
}; // Species struct

} // namespace Gas
} // namespace Model

#endif
