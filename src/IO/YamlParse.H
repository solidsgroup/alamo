#ifndef IO_YAMLPARSE
#define IO_YAMLPARSE

#include "Util/Util.H"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <regex>
#include <yaml-cpp/yaml.h>

namespace CanteraYAML {

// ---------------------------
// Data structures
// ---------------------------

struct Phase {
    std::string name;
    std::string thermoModel;
    std::vector<std::string> species;
};

struct Species {
    std::string name;
    std::map<std::string, double> composition;
    std::string thermoModel;
    std::vector<double> thermoTemp;
    // Thermo: nested vectors for NASA coefficients
    std::vector<std::vector<double>> thermoData;
    double transportLJdepth = 0.0;
    double transportLJdiameter = 0.0;
};

struct Reaction {
    std::string equation;
    std::string type = "arrhenius";
    bool third_body = false;
    bool reversible = false;
    bool falloff = false;
    std::map<std::string, double> reactants;
    std::map<std::string, double> products;
    std::map<std::string, double> orders;
    // Arrhenius constants (high values if falloff)
    double A = 0.0;
    double b = 0.0;
    double E = 0.0;
    // Arrhenius constants for low if falloff
    double A2 = 0.0;
    double b2 = 0.0;
    double E2 = 0.0;
    double F = 1.0;
    std::map<std::string, double> efficiencies;
    std::map<std::string, double> Troe;
};

struct CanteraKinetics {
    std::vector<Phase> phase;
    std::vector<Species> species;
    std::vector<Reaction> reaction;
};

// ---------------------------
// Utility functions
// ---------------------------

// Convert a YAML node to double safely
inline double asDoubleSafe(const YAML::Node& node) {
    if (!node) return 0.0;
    return node.as<double>();
}

// Convert a YAML node to string safely
inline std::string asStringSafe(const YAML::Node& node) {
    if (!node) return "";
    return node.as<std::string>();
}

// Parse inline map {key: value(double), ...}
inline std::map<std::string, double> parseMap(const YAML::Node& node) {
    std::map<std::string, double> result;
    for (auto it = node.begin(); it != node.end(); ++it) {
        result[it->first.as<std::string>()] = double(it->second.as<double>());
    }
    return result;
}

// Parse inline map {key: value(string), ...}
inline std::map<std::string, std::string> parseMapString(const YAML::Node& node) {
    std::map<std::string, std::string> result;
    for (auto it = node.begin(); it != node.end(); ++it) {
        result[it->first.as<std::string>()] = it->second.as<std::string>();
    }
    return result;
}

// Parse a vector of doubles from a YAML sequence or inline list
inline std::vector<double> parseDoubleList(const YAML::Node& node) {
    std::vector<double> result;
    for (const auto& n : node) {
        result.push_back(n.as<double>());
    }
    return result;
}

// Helper: trim whitespace
inline std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\n\r");
    auto end   = s.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// Helper: split string by delimiter
inline std::vector<std::string> split(const std::string& s, const std::string& delim) {
    std::vector<std::string> elems;
    if (delim.empty()) {  // handle degenerate case
        elems.push_back(trim(s));
        return elems;
    }

    size_t start = 0;
    size_t end = 0;

    while ((end = s.find(delim, start)) != std::string::npos) {
        std::string token = s.substr(start, end - start);
        elems.push_back(trim(token));
        start = end + delim.length();
    }

    // Add final segment
    if (start < s.size())
        elems.push_back(trim(s.substr(start)));
    else if (!s.empty() && s.rfind(delim) == s.size() - delim.size())
        elems.push_back("");  // optional: preserve trailing empty element
    return elems;
}

// Replace all occurrences of "(+Species)" with " + +Species"
inline std::string replaceThirdBody(const std::string& eq) {
    // Regex explanation:
    //   \(           : literal '('
    //   \s*\+        : optional whitespace, then a literal '+'
    //   ([A-Za-z0-9_]+) : capture species name (letters, digits, underscore)
    //   \s*\)        : optional whitespace, then ')'
    //
    // Example matches: "(+M)", "( +Ar )", "(+CO2)", "(+X_123)"
    static const std::regex pattern(R"(\(\s*\+([A-Za-z0-9_]+)\s*\))");

    // Replace with "+ +\1" where \1 is the captured species name
    return std::regex_replace(eq, pattern, "+ +$1");
}

// Parse string to get number of species
auto parseSide = [](const std::string& side) -> std::map<std::string, double> {
    std::map<std::string, double> species;
    auto terms = split(side, " + ");
    for (auto& term : terms) {
        if (term.empty()) continue;
        std::stringstream ss(term);
        double coeff;
        std::string name;

        // Try reading a number
        if ((ss >> coeff)) {
            if (!(ss >> name)) {
                // If we didn’t get a name, maybe the coefficient was actually the species name
                name = std::to_string(coeff);
                coeff = 1.0;
            }
        } else {
            // No coefficient found
            coeff = 1.0;
            name = term;
        }

        name = trim(name);
        if (name.empty()) continue;
        species[name] += coeff; // sum duplicates
    }
    return species;
};

// ---------------------------
// Main parser class
// ---------------------------

class CanteraParser {
public:
    YAML::Node root;

    bool loadFile(const std::string& filename) {
        try {
            root = YAML::LoadFile(filename);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Failed to load YAML file: " << e.what() << "\n";
            return false;
        }
    }
    // Extract units
    // TODO: convert provided units to SI
    std::map<std::string, std::string> getUnits() const {
        std::map<std::string, std::string> result;
        if (!root["units"]) return result;
        else return parseMapString(root["units"]);
    }        

    // Extract phases
    std::vector<Phase> getPhases() const {
        std::vector<Phase> result;
        if (!root["phases"]) return result;
        for (const auto& node : root["phases"]) {
            Phase p;
            p.name = asStringSafe(node["name"]);
            p.thermoModel = asStringSafe(node["thermo"]);
            if (node["species"]) {
                for (const auto& s : node["species"]) {
                    p.species.push_back(s.as<std::string>());
                }
            }
            result.push_back(p);
        }
        return result;
    }

    // Extract species
    std::vector<Species> getSpecies() const {
        auto units = getUnits();
        std::vector<Species> result;
        if (!root["species"]) return result;
        for (const auto& node : root["species"]) {
            Species s;
            s.name = asStringSafe(node["name"]);
            if (node["composition"]) s.composition = parseMap(node["composition"]);
            if (node["thermo"] && node["thermo"]["data"]) {
                s.thermoModel = asStringSafe(node["thermo"]["model"]);
                s.thermoTemp = parseDoubleList(node["thermo"]["temperature-ranges"]);
                for (const auto& arr : node["thermo"]["data"]) {
                    s.thermoData.push_back(parseDoubleList(arr));
                }
            }
            if (node["transport"]) {
                //boltzmann = 1.380649e-23;
                // We want well-depth in E/boltzmann = T (units Kelvin)
                // diameter in Angstroms
                // These are the required units for the mixture models in Hydro
                s.transportLJdepth = asDoubleSafe(node["transport"]["well-depth"]);///boltzmann;
                s.transportLJdiameter = asDoubleSafe(node["transport"]["diameter"]);
            }
            result.push_back(s);
        }
        return result;
    }

    // Extract reactions
    std::vector<Reaction> getReactions() const {
        auto units = getUnits();
        std::vector<Reaction> result;
        if (!root["reactions"]) return result;
        for (const auto& node : root["reactions"]) {
            Reaction r;
            r.equation = asStringSafe(node["equation"]);
            if (node["rate-constant"]) {
                const auto& rc = node["rate-constant"];
                r.A = asDoubleSafe(rc["A"]);
                r.b = asDoubleSafe(rc["b"]);
                r.E = asDoubleSafe(rc["Ea"]);
            }
            if (node["high-P-rate-constant"]) {
                const auto& rc = node["high-P-rate-constant"];
                r.A = asDoubleSafe(rc["A"]);
                r.b = asDoubleSafe(rc["b"]);
                r.E = asDoubleSafe(rc["Ea"]);
            }
            if (node["low-P-rate-constant"]) {
                const auto& rc = node["low-P-rate-constant"];
                r.A2 = asDoubleSafe(rc["A"]);
                r.b2 = asDoubleSafe(rc["b"]);
                r.E2 = asDoubleSafe(rc["Ea"]);
            }
            if (node["type"]) r.type = asStringSafe(node["type"]);
            if (r.type == "three-body") r.third_body = true;
            else if (r.type == "falloff") {
                r.falloff = true;
                r.third_body = true;
            }
            std::string lhs, rhs;
            std::string eq_mod = replaceThirdBody(r.equation);
            size_t found_forward = eq_mod.find("=>");
            size_t found_reversible = eq_mod.find("<=>");
            if (found_reversible != std::string::npos) {
                r.reversible = true;
                lhs = eq_mod.substr(0, found_reversible);
                rhs = eq_mod.substr(found_reversible + 3);
            } else if (found_forward != std::string::npos) {
                lhs = eq_mod.substr(0, found_forward);
                rhs = eq_mod.substr(found_forward + 2);
            } else {
                Util::Message(INFO, "No valid arrow (=> or <=>) found in reaction.");
                Util::Abort(INFO);
            }
            r.reactants = parseSide(lhs);
            r.products = parseSide(rhs);
            if (node["orders"]) r.orders = parseMap(node["orders"]);
            if (node["efficiencies"]) r.efficiencies = parseMap(node["efficiencies"]);
            if (node["Troe"]) {
                r.Troe = parseMap(node["Troe"]);
                r.F = 1.0; // TODO: update this based on Troe Fcent
            }

            // Convert to proper units
            double order = 0.0;
            for (const auto& pair : r.reactants) {
                if (node["orders"] and r.orders[pair.first]) order += r.orders[pair.first]; 
                else order += pair.second;
            }
            order = (order - 1.0);

            if (units["activation-energy"] == "cal/mol") {
                r.E *= 4184.0; // convert to J/kmol
                r.E2 *= 4184.0; // convert to J/kmol
            }

            if (units["length"] == "cm") {
                // cm^3 to m^3
                r.A *= pow(1e-2, 3*order);
                r.A2 *= pow(1e-2, 3*order);

            }
            if (units["quantity"] == "mol") {
                // 1/mol to 1/kmo1
                r.A *= pow(1.0e3, order);
                r.A2 *= pow(1.0e3, order);
            }
                  
            result.push_back(r);
        }
        return result;
    }
};

CanteraKinetics ParseYaml(std::string &filename);

} // namespace CanteraYAML

#endif
