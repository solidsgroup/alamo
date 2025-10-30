#include <iostream>
#include <map>
#include "IO/YamlParse.H"
#include "Util/Util.H"

CanteraYAML::CanteraKinetics CanteraYAML::ParseYaml(std::string &filename) {
    CanteraYAML::CanteraKinetics kinetics;
    CanteraYAML::CanteraParser parser;

    std::string yamlfile = filename;
    if (!parser.loadFile(yamlfile)) {
        Util::Message(INFO, "Invalid yaml file ", yamlfile);
        Util::Abort(INFO);
    }

    auto units = parser.getUnits();
#ifdef AMREX_DEBUG
    if (units.size() > 0) {
        for (const auto& pair : units) {
            std::cout << "Units: " << pair.first << ": " << pair.second << "\n";
        }
    }
#endif

    auto phases = parser.getPhases();
    kinetics.phase = phases;
#ifdef AMREX_DEBUG
    std::cout << "Phases: " << phases.size() << "\n";
    for (auto& p : phases) {
        std::cout << " - " << p.name << " (" << p.species.size() << " species)\n";
    }
#endif

    auto species = parser.getSpecies();
    kinetics.species = species;
#ifdef AMREX_DEBUG
    std::cout << "Species: " << species.size() << "\n";
    if (!species.empty()) {
        for (size_t n=0; n<species.size(); ++n) {
            std::cout << "Species[" << n << "]: " << species[n].name << "\t" << species[n].transportLJdepth << "\t" << species[n].transportLJdiameter << "\n";
            for (size_t nt=0; nt<species[n].thermoTemp.size(); ++nt) {
                std::cout << "\t" << species[n].thermoTemp[nt];
            }
            std::cout << "\n";
            for (size_t nt=0; nt<species[n].thermoData[0].size(); ++nt) {
                std::cout << "\t\t" << species[n].thermoData[0][nt] << "\t" << species[n].thermoData[1][nt] << "\n";
            }
            for (const auto& pair : species[n].composition) {
                std::cout << "\t" << pair.first << ": " << pair.second << "\n";
            }
        }
    }
#endif

    auto reactions = parser.getReactions();
    kinetics.reaction = reactions;
#ifdef AMREX_DEBUG
    std::cout << "Reactions: " << reactions.size() << "\n";
    for (size_t n=0; n<reactions.size(); ++n) {
        std::cout << n << " - " << reactions[n].equation << "\n";
        std::cout << "  [A=" << reactions[n].A
                  << ", b=" << reactions[n].b
                  << ", Ea=" << reactions[n].E << "]\n";
        std::cout << "\t Type: " << reactions[n].type << "\n";
        std::cout << "\t Third Body: " << reactions[n].third_body << "\n";
        std::cout << "\t reversible: " << reactions[n].reversible << "\n";
        std::cout << "\t falloff: " << reactions[n].falloff << "\n";
        for (const auto& pair : reactions[n].reactants) {
            std::cout << "\t Reactants: " << pair.first << " " << pair.second << "\n";
        }
        for (const auto& pair : reactions[n].products) {
            std::cout << "\t Products: " << pair.first << " " << pair.second << "\n";
        }
        for (const auto& pair : reactions[n].efficiencies) {
            std::cout << "\t Efficiencies: " << pair.first << " " << pair.second << "\n";
        }
        for (const auto& pair : reactions[n].Troe) {
            std::cout << "\t Troe: " << pair.first << " " << pair.second << "\n";
        }
        for (const auto& pair : reactions[n].orders) {
            std::cout << "\t Orders: " << pair.first << " " << pair.second << "\n";
        }
    }
#endif

    return kinetics;
}
