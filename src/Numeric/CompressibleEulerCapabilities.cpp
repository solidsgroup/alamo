// CompressibleEulerCapabilities.cpp
#include "Numeric/IntegratorVariableAccessLayer.H"

namespace Numeric {
namespace CompressibleEuler {

SolverCapabilities::MethodSupport 
CompressibleEulerCapabilities::supportsFluxReconstruction(FluxReconstructionType method) const {
    switch(method) {
        case FluxReconstructionType::FirstOrder:
            return MethodSupport::Supported();
        case FluxReconstructionType::WENO:
            return MethodSupport::Supported();
        default:
            return MethodSupport::Unsupported(
                {"Unsupported reconstruction method"},
                {"Use FirstOrder or WENO reconstruction"}
            );
    }
}

SolverCapabilities::MethodSupport 
CompressibleEulerCapabilities::supportsFluxScheme(FluxScheme scheme) const {
    switch(scheme) {
        case FluxScheme::LocalLaxFriedrichs:
            return MethodSupport::Supported();
        case FluxScheme::HLLC:
            return MethodSupport::Supported();
        default:
            return MethodSupport::Unsupported(
                {"Unknown flux scheme"},
                {"Use standard flux methods"}
            );
    }
}

SolverCapabilities::MethodSupport 
CompressibleEulerCapabilities::supportsTimeSteppingScheme(TimeSteppingSchemeType scheme) const {
    switch(scheme) {
        case TimeSteppingSchemeType::ForwardEuler:
        case TimeSteppingSchemeType::RK3:
            return MethodSupport::Supported();
        default:
            return MethodSupport::Unsupported(
                {"Unsupported time stepping scheme"},
                {"Use ForwardEuler or RK3"}
            );
    }
}

SolverCapabilities::MethodSupport 
CompressibleEulerCapabilities::supportsReconstructionMode(ReconstructionMode mode) const {
    switch(mode) {
        case ReconstructionMode::Primitive:
        case ReconstructionMode::Conservative:
        case ReconstructionMode::Characteristic:
            return MethodSupport::Supported();
        default:
            return MethodSupport::Unsupported(
                {"Unsupported reconstruction mode"},
                {"Use Primitive, Conservative, or Characteristic"}
            );
    }
}

SolverCapabilities::MethodSupport 
CompressibleEulerCapabilities::supportsWenoVariant(WenoVariant variant) const {
    switch(variant) {
        case WenoVariant::WENOJS5:
            return MethodSupport::Supported();
        case WenoVariant::WENOZ5:
            return MethodSupport::Supported();
        case WenoVariant::WENOJS3:
            return MethodSupport::Supported();
        default:
            return MethodSupport::Unsupported(
                {"Unsupported WENO variant"},
                {"Use WENOJS5 or WENOZ5"}
            );
    }
}

SolverCapabilities::MethodValidationResult 
CompressibleEulerCapabilities::validateMethodCombination(
    FluxReconstructionType fluxReconstruction,
    FluxScheme fluxScheme,
    TimeSteppingSchemeType timeSteppingScheme,
    ReconstructionMode reconstructionMode,
    WenoVariant wenoVariant) const 
{
    MethodValidationResult result;
    result.isValid = true;

    // Add specific validation logic here
   /* if (fluxReconstruction == FluxReconstructionType::WENO && 
        fluxScheme == FluxScheme::HLLC) {
        result.isValid = false;
        result.warnings.push_back(
            "WENO reconstruction with HLLC flux might be numerically unstable"
        );
    }*/

    return result;
}

SolverCapabilities::DefaultConfiguration 
CompressibleEulerCapabilities::getDefaultConfiguration() const {
    return {
        FluxReconstructionType::WENO,
        FluxScheme::LocalLaxFriedrichs,
        TimeSteppingSchemeType::RK3,
        ReconstructionMode::Characteristic,
        WenoVariant::WENOJS5
    };
}

std::shared_ptr<GenericVariableAccessor>
CompressibleEulerCapabilities::createVariableAccessor(
    ReconstructionMode mode,
    int numGhostCells) const 
{
    return std::make_shared<CompressibleEulerVariableAccessor>(
        numGhostCells, mode);
}

} // namespace CompressibleEuler
} // namespace Numeric
