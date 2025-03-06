// NumericFactory.H
#ifndef NUMERIC_NUMERIC_FACTORY_H
#define NUMERIC_NUMERIC_FACTORY_H

#include <memory>
#include <string>
#include <map>
#include <set>
#include <typeinfo>
#include <stdexcept>
#include <vector>
#include "Numeric/WENOReconstruction.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/SolverCapabilities.H"

namespace Numeric {

/**
 * @brief Factory for creating numeric methods with comprehensive validation
 */
class NumericFactory {
public:
    // Flux Reconstruction Factory Method 
    template<typename SolverType>
    static std::shared_ptr<Numeric::FluxReconstruction<SolverType>> 
    createReconstructionMethod(
        FluxReconstructionType method = FluxReconstructionType::WENO,
        WenoVariant variant = WenoVariant::WENOJS5) {
        
        // Get solver capabilities
        auto caps = SolverCapabilitiesRegistry::getInstance()
            .getCapabilities(getSolverIdentifier<SolverType>());
        
        // Validate method support
        auto reconstructionSupport = caps->supportsFluxReconstruction(method);
        if (!reconstructionSupport) {
            throw std::runtime_error(
                "Unsupported Flux Reconstruction: " + 
                toString(method) + "\n" +
                "Reasons: " + 
                (reconstructionSupport.reasons.empty() ? 
                    "No specific reasons provided" : 
                    reconstructionSupport.reasons[0])
            );
        }

        // Create appropriate reconstruction method
        switch(method) {
            case FluxReconstructionType::WENO: {
                auto wenoSupport = caps->supportsWenoVariant(variant);
                if (!wenoSupport) {
                    // Fallback to default WENO variant
                    variant = WenoVariant::WENOJS5;
                }
                
                switch(variant) {
                    case WenoVariant::WENOJS5:
                        return std::make_shared<Numeric::WENOJS5<SolverType>>();
                    case WenoVariant::WENOZ5:
                        return std::make_shared<Numeric::WENOZ5<SolverType>>();
                    default:
                        throw std::runtime_error("Unsupported WENO variant");
                }
            }
            case FluxReconstructionType::FirstOrder:
                return std::make_shared<Numeric::FirstOrderReconstruction<SolverType>>();
            
            default:
                throw std::runtime_error("Unsupported reconstruction method");
        }
    }

    // Flux Method Factory Method
    template<typename SolverType>
    static std::shared_ptr<Numeric::FluxMethod<SolverType>> 
    createFluxMethod(
        FluxScheme scheme = FluxScheme::LocalLaxFriedrichs) {
        
        // Get solver capabilities
        auto caps = SolverCapabilitiesRegistry::getInstance()
            .getCapabilities(getSolverIdentifier<SolverType>());
        
        // Validate method support
        auto fluxSupport = caps->supportsFluxScheme(scheme);
        if (!fluxSupport) {
            throw std::runtime_error(
                "Unsupported Flux Scheme: " + 
                toString(scheme) + "\n" +
                "Reasons: " + 
                (fluxSupport.reasons.empty() ? 
                    "No specific reasons provided" : 
                    fluxSupport.reasons[0])
            );
        }

        // Create appropriate flux method
        switch(scheme) {
            case FluxScheme::LocalLaxFriedrichs:
                return std::make_shared<Numeric::LocalLaxFriedrichsMethod<SolverType>>();
            case FluxScheme::HLLC:
                // TODO: Implement HLLC method
                throw std::runtime_error("HLLC method not yet implemented");
            default:
                throw std::runtime_error("Unsupported flux method");
        }
    }

    // Time Stepping Scheme Factory Method
    template<typename SolverType>
    static std::shared_ptr<Numeric::TimeSteppingScheme<SolverType>> 
    createTimeSteppingScheme(
        TimeSteppingSchemeType scheme = TimeSteppingSchemeType::RK3) {
        
        // Get solver capabilities
        auto caps = SolverCapabilitiesRegistry::getInstance()
            .getCapabilities(getSolverIdentifier<SolverType>());
        
        // Validate method support
        auto timeStepSupport = caps->supportsTimeSteppingScheme(scheme);
        if (!timeStepSupport) {
            throw std::runtime_error(
                "Unsupported Time Stepping Scheme: " + 
                toString(scheme) + "\n" +
                "Reasons: " + 
                (timeStepSupport.reasons.empty() ? 
                    "No specific reasons provided" : 
                    timeStepSupport.reasons[0])
            );
        }

        // Create appropriate time stepping scheme
        switch(scheme) {
            case TimeSteppingSchemeType::ForwardEuler:
                return std::make_shared<Numeric::EulerForwardScheme<SolverType>>();
            case TimeSteppingSchemeType::RK3:
                return std::make_shared<Numeric::RK3Scheme<SolverType>>();
            default:
                throw std::runtime_error("Unsupported time stepping scheme");
        }
    }

    // Variable Accessor Factory Method
    template<typename SolverType>
    static std::shared_ptr<Numeric::GenericVariableAccessor>
    createVariableAccessor(
        ReconstructionMode mode = ReconstructionMode::Primitive,
        int numGhostCells = 2) {
        
        // Get solver capabilities
        auto caps = SolverCapabilitiesRegistry::getInstance()
            .getCapabilities(getSolverIdentifier<SolverType>());
        
        // Validate reconstruction mode support
        auto modeSupport = caps->supportsReconstructionMode(mode);
        if (!modeSupport) {
            throw std::runtime_error(
                "Unsupported Reconstruction Mode: " + 
                toString(mode) + "\n" +
                "Reasons: " + 
                (modeSupport.reasons.empty() ? 
                    "No specific reasons provided" : 
                    modeSupport.reasons[0])
            );
        }

        // Create variable accessor
        return caps->createVariableAccessor(mode, numGhostCells);
    }

    // String Conversion Utilities
    static std::string toString(FluxReconstructionType method);
    static std::string toString(FluxScheme scheme);
    static std::string toString(TimeSteppingSchemeType scheme);
    static std::string toString(ReconstructionMode mode);
    static std::string toString(WenoVariant variant);

    // Enum Parsing Utilities
    static FluxReconstructionType parseFluxReconstruction(const std::string& str);
    static FluxScheme parseFluxScheme(const std::string& str);
    static TimeSteppingSchemeType parseTimeSteppingScheme(const std::string& str);
    static ReconstructionMode parseReconstructionMode(const std::string& str);
    static WenoVariant parseWenoVariant(const std::string& str);
};

// Specialization for ScimitarX defined in NumericFactory.cpp
template<>
inline std::string getSolverIdentifier<Integrator::ScimitarX>() {
    return "CompressibleEuler";
}

} // namespace Numeric

#endif // NUMERIC_NUMERIC_FACTORY_H
