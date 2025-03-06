// NumericFactory.cpp
#include "Numeric/NumericFactory.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/WENOReconstruction.H"
#include "Numeric/TimeStepper.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/IntegratorVariableAccessLayer.H"

namespace Numeric {
    // String Conversion Implementations
    std::string NumericFactory::toString(FluxReconstructionType method) {
        switch(method) {
            case FluxReconstructionType::FirstOrder: return "FirstOrder";
            case FluxReconstructionType::WENO: return "WENO";
            default: throw std::runtime_error("Unknown FluxReconstruction");
        }
    }

    std::string NumericFactory::toString(WenoVariant variant) {
        switch(variant) {
            case WenoVariant::WENOJS5: return "WENOJS5";
            case WenoVariant::WENOZ5: return "WENOZ5";
            default: throw std::runtime_error("Unknown WENO Variant");
        }
    }

    std::string NumericFactory::toString(FluxScheme scheme) {
        switch(scheme) {
            case FluxScheme::LocalLaxFriedrichs: return "LocalLaxFriedrichs";
            case FluxScheme::HLLC: return "HLLC";
            default: throw std::runtime_error("Unknown FluxScheme");
        }
    }

    std::string NumericFactory::toString(TimeSteppingSchemeType scheme) {
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
    FluxReconstructionType NumericFactory::parseFluxReconstruction(const std::string& str) {
        if (str == "FirstOrder") return FluxReconstructionType::FirstOrder;
        if (str == "WENO") return FluxReconstructionType::WENO;
        throw std::runtime_error("Invalid FluxReconstruction: " + str);
    }

    WenoVariant NumericFactory::parseWenoVariant(const std::string& str) {
        if (str == "WENOJS5") return WenoVariant::WENOJS5;
        if (str == "WENOZ5") return WenoVariant::WENOZ5;
        throw std::runtime_error("Invalid WENO Variant: " + str);
    }

    FluxScheme NumericFactory::parseFluxScheme(const std::string& str) {
        if (str == "LocalLaxFriedrichs") return FluxScheme::LocalLaxFriedrichs;
        if (str == "HLLC") return FluxScheme::HLLC;
        throw std::runtime_error("Invalid FluxScheme: " + str);
    }

    TimeSteppingSchemeType NumericFactory::parseTimeSteppingScheme(const std::string& str) {
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
}
