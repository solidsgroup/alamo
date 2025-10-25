#ifndef IO_YAMLPARSE
#define IO_YAMLPARSE

#include "Util/Util.H"
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace CanteraYAML {

// ---------------------------
// Data structures
// ---------------------------

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

struct Phase {
    std::string name;
    std::string thermoModel;
    std::vector<std::string> species;
};

struct Reaction {
    std::string equation;
    std::string type = "arrhenius";
    bool third_body = false;
    bool reversible = false;
    std::map<std::string, double> reactants;
    std::map<std::string, double> products;
    double A = 0.0;
    double b = 0.0;
    double Ea = 0.0;
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
inline std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) elems.push_back(trim(item));
    return elems;
}

// Parse string to get number of species
auto parseSide = [](const std::string& side) -> std::map<std::string, double> {
    std::map<std::string, double> species;
    auto terms = split(side, '+');
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
                s.transportLJdepth = asDoubleSafe(node["transport"]["well-depth"]);
                s.transportLJdiameter = asDoubleSafe(node["transport"]["diameter"]);
            }
            result.push_back(s);
        }
        return result;
    }

    // Extract reactions
    std::vector<Reaction> getReactions() const {
        std::vector<Reaction> result;
        if (!root["reactions"]) return result;
        for (const auto& node : root["reactions"]) {
            Reaction r;
            r.equation = asStringSafe(node["equation"]);
            if (node["rate-constant"]) {
                const auto& rc = node["rate-constant"];
                r.A = asDoubleSafe(rc["A"]);
                r.b = asDoubleSafe(rc["b"]);
                r.Ea = asDoubleSafe(rc["Ea"]);
            }
            if (node["type"]) r.type = asStringSafe(node["type"]);
            if (r.type == "three-body") r.third_body = true;
            size_t found_forward = r.equation.find("=>");
            size_t found_reversible = r.equation.find("<=>");
            std::string lhs, rhs;
            if (found_reversible != r.equation.npos) {
                r.reversible = true;
                lhs = r.equation.substr(0, found_reversible);
                rhs = r.equation.substr(found_reversible + 3);
            } else if (found_forward != r.equation.npos) {
                lhs = r.equation.substr(0, found_forward);
                rhs = r.equation.substr(found_forward + 2);
            } else {
                Util::Message(INFO, "No valid arrow (=> or <=>) found in reaction.");
                Util::Abort(INFO);
            }
            r.reactants = parseSide(lhs);
            r.products = parseSide(rhs);
            result.push_back(r);
        }
        return result;
    }
};

class CanteraKinetics
{
public:
    std::vector<Phase> phase;
    std::vector<Species> species;
    std::vector<Reaction> reaction;
};

CanteraKinetics ParseYaml(std::string &filename);

} // namespace CanteraYAML

#endif
