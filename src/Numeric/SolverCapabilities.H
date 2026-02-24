// SolverCapabilities.H
#ifndef NUMERIC_SOLVER_CAPABILITIES_H
#define NUMERIC_SOLVER_CAPABILITIES_H

#include "NumericTypes.H"
#include <memory>
#include <string>
#include <vector>

namespace Numeric {

// Forward declaration for GenericVariableAccessor
class GenericVariableAccessor;

class SolverCapabilities {
public:

        struct MethodSupport {
            bool supported = false;
            std::vector<std::string> reasons;
            std::vector<std::string> suggestions;

            operator bool() const { return supported; }

            static MethodSupport Supported() {
                return {true, {}, {}};
            }

            static MethodSupport Unsupported(
                const std::vector<std::string>& reasons = {},
                const std::vector<std::string>& suggestions = {}) {
                return {false, reasons, suggestions};
            }
        };

        struct MethodValidationResult {
            bool isValid = false;
            std::vector<std::string> warnings;
            std::vector<std::string> errors;
        };
        
        struct DefaultConfiguration {
            FluxReconstructionType fluxReconstruction;
            FluxScheme fluxScheme;
            TimeSteppingSchemeType timeSteppingScheme;
            ReconstructionMode reconstructionMode;
            WenoVariant wenoVariant;
        };

    // Virtual Destructor
    virtual ~SolverCapabilities() = default;

    // Core Solver Identification
    virtual std::string getIdentifier() const = 0;
    virtual std::string getDescription() const = 0;

    // Method Support Queries with Detailed Feedback
    virtual MethodSupport supportsFluxReconstruction(
        FluxReconstructionType method) const = 0;

    virtual MethodSupport supportsFluxScheme(
        FluxScheme scheme) const = 0;

    virtual MethodSupport supportsTimeSteppingScheme(
        TimeSteppingSchemeType scheme) const = 0;

    virtual MethodSupport supportsReconstructionMode(
        ReconstructionMode mode) const = 0;

    virtual MethodSupport supportsWenoVariant(
        WenoVariant variant) const = 0;

    // Comprehensive Method Combination Validation
    virtual MethodValidationResult validateMethodCombination(
        FluxReconstructionType fluxReconstruction,
        FluxScheme fluxScheme,
        TimeSteppingSchemeType timeSteppingScheme,
        ReconstructionMode reconstructionMode,
        WenoVariant wenoVariant = WenoVariant::WENOJS5
    ) const = 0;

    // Default Configuration Recommendations
    virtual DefaultConfiguration getDefaultConfiguration() const = 0;

    // Variable Accessor Creation
    virtual std::shared_ptr<GenericVariableAccessor> createVariableAccessor(
        ReconstructionMode mode,
        int numGhostCells = 2
    ) const = 0;
};

// Registry for managing solver capabilities
class SolverCapabilitiesRegistry {
public:
    // Singleton Access
    static SolverCapabilitiesRegistry& getInstance() {
        static SolverCapabilitiesRegistry instance;
        return instance;
    }

    // Register Solver Capabilities
    void registerCapabilities(std::shared_ptr<SolverCapabilities> caps) {
        solverCapabilities[caps->getIdentifier()] = std::move(caps);
    }

    // Retrieve Solver Capabilities
    std::shared_ptr<SolverCapabilities> getCapabilities(
        const std::string& identifier) const {
        auto it = solverCapabilities.find(identifier);
        if (it == solverCapabilities.end()) {
            throw std::runtime_error("No capabilities found for solver: " + identifier);
        }
        return it->second;
    }

private:
    std::map<std::string, std::shared_ptr<SolverCapabilities>> solverCapabilities;
};

} // namespace Numeric

#endif // NUMERIC_SOLVER_CAPABILITIES_H
