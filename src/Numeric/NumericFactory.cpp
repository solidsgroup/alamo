// NumericFactory.cpp
#include "NumericFactory.H"
#include "WENOReconstruction.H"
#include "TimeStepper.H"
#include "FluxHandler.H"
#include "IntegratorVariableAccessLayer.H"

namespace Numeric {
    // Singleton Implementation
    SolverMetadataRegistry& SolverMetadataRegistry::getInstance() {
        static SolverMetadataRegistry instance;
        return instance;
    }

    // String Conversion Implementations
    std::string NumericFactory::toString(FluxReconstruction method) {
        switch(method) {
            case FluxReconstructionType::FirstOrder: return "FirstOrder";
            case FluxReconstructionType::WENO: return "WENO";
            default: throw std::runtime_error("Unknown FluxReconstruction");
        }
    }

    std::string NumericFactory::toString(FluxScheme scheme) {
        switch(scheme) {
            case FluxScheme::LocalLaxFriedrichs: return "LocalLaxFriedrichs";
            case FluxScheme::HLLC: return "HLLC";
            default: throw std::runtime_error("Unknown FluxScheme");
        }
    }

    std::string NumericFactory::toString(TimeSteppingScheme scheme) {
        switch(scheme) {
            case TimeSteppingSchemeType::ForwardEuler: return "ForwardEuler";
            case TimeSteppingSchemeType::RK3: return "RK3";
            default: throw std::runtime_error("Unknown TimeSteppingScheme");
        }
    }

    std::string NumericFactory::toString(ReconstructionMode mode) {
        switch(mode) {
            case ReconstructionMode::Primitive: return "Primitive";
            case ReconstructionMode::Conservative: return "Conservative";
            case ReconstructionMode::Characteristic: return "Characteristic";
            default: throw std::runtime_error("Unknown ReconstructionMode");
        }
    }

    // Parsing Implementations
    FluxReconstruction NumericFactory::parseFluxReconstruction(const std::string& str) {
        if (str == "FirstOrder") return FluxReconstructionType::FirstOrder;
        if (str == "WENO") return FluxReconstructionType::WENO;
        throw std::runtime_error("Invalid FluxReconstruction: " + str);
    }

    FluxScheme NumericFactory::parseFluxScheme(const std::string& str) {
        if (str == "LocalLaxFriedrichs") return FluxScheme::LocalLaxFriedrichs;
        if (str == "HLLC") return FluxScheme::HLLC;
        throw std::runtime_error("Invalid FluxScheme: " + str);
    }

    FluxMethod NumericFactory::parseTimeSteppingScheme(const std::string& str) {
        if (str == "ForwardEuler") return TimeSteppingSchemeType::ForwardEuler;
        if (str == "RK3") return TimeSteppingSchemeType::RK3;
        throw std::runtime_error("Invalid TimeSteppingScheme: " + str);
    }

    ReconstructionMode NumericFactory::parseReconstructionMode(const std::string& str) {
        if (str == "Primitive") return ReconstructionMode::Primitive;
        if (str == "Conservative") return ReconstructionMode::Conservative;
        if (str == "Characteristic") return ReconstructionMode::Characteristic;
        throw std::runtime_error("Invalid ReconstructionMode: " + str);
    }

std::shared_ptr<GenericVariableAccessor> NumericFactory::createVariableAccessor(
    const std::string& solverTypeName,
    ReconstructionMode mode,
    int numGhostCells) {
    
    // Use a map to look up the correct factory function
    static std::map<std::string, std::function<std::shared_ptr<GenericVariableAccessor>(ReconstructionMode, int)>> factories = {
        {"CompressibleEuler", [](ReconstructionMode mode, int ghosts) {
            return std::make_shared<CompressibleEuler::CompressibleEulerVariableAccessor>(ghosts, mode);
        }},
        // Add other solver types as needed
    };
    
    auto it = factories.find(solverTypeName);
    if (it != factories.end()) {
        return it->second(mode, numGhostCells);
    }
    
    // Default to CompressibleEuler if not found
    Util::Warning(INFO, "Unknown solver type, defaulting to CompressibleEuler");
    return std::make_shared<CompressibleEuler::CompressibleEulerVariableAccessor>(numGhostCells, mode);
}    
    // Implementations of Factory Methods will go here
    // These will be template method implementations
    // that use the SolverMetadataRegistry to validate methods
}
