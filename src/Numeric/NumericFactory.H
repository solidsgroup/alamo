// NumericFactory.H
#ifndef NUMERIC_NUMERIC_FACTORY_H
#define NUMERIC_NUMERIC_FACTORY_H

#include <memory>
#include <string>
#include <map>
#include <set>
#include <typeinfo>
#include <stdexcept>

// Forward declarations
namespace Numeric {
    template <typename IntegratorType> class FluxReconstruction;
    template <typename IntegratorType> class FluxMethod;
    template <typename IntegratorType> class TimeSteppingScheme;
    class GenericVariableAccessor;
}

namespace Numeric {
    // Enum Classes for Numeric Methods
    enum class FluxReconstructionType {
        FirstOrder,    // Basic first-order upwind reconstruction
        WENO,          // Weighted Essentially Non-Oscillatory
    };

    enum class WenoVariant {
        WENOJS5,           // WENO Jiang-Shu 5th Order
        WENOZ5,            // WENO Z 5th Order
    };

    enum class FluxScheme {
        LocalLaxFriedrichs,  // Local Lax-Friedrichs (Rusanov) flux
        HLLC,                // Harten-Lax-van Leer Contact flux
    };

    enum class TimeSteppingSchemeType {
        ForwardEuler,  // First-order explicit Euler method
        RK3,           // Third-order Runge-Kutta
    };

    enum class ReconstructionMode {
        Primitive,       // Reconstruct in primitive variable space
        Conservative,    // Reconstruct in conservative variable space
        Characteristic   // Reconstruct in characteristic variable space
    };

    // Abstract Base Class for Solver Metadata
    class SolverMetadata {
    public:
        virtual ~SolverMetadata() = default;
        
        // Solver Identification
        virtual std::string getName() const = 0;
        
        // Compatibility Checks
        virtual bool supportsFluxReconstruction(FluxReconstructionType method) const = 0;
        virtual bool supportsFluxScheme(FluxScheme scheme) const = 0;
        virtual bool supportsTimeSteppingScheme(TimeSteppingSchemeType scheme) const = 0;
        virtual bool supportsReconstructionMode(ReconstructionMode mode) const = 0;
        
        // Default Method Recommendations
        virtual FluxReconstructionType getDefaultFluxReconstruction() const = 0;
        virtual FluxScheme getDefaultFluxScheme() const = 0;
        virtual TimeSteppingSchemeType getDefaultTimeSteppingScheme() const = 0;
        virtual ReconstructionMode getDefaultReconstructionMode() const = 0;
    };

    // Solver Metadata Registry
    class SolverMetadataRegistry {
    public:
        // Singleton Access
        static SolverMetadataRegistry& getInstance();

        // Register Solver Metadata
        template<typename SolverType>
        void registerSolverMetadata(std::unique_ptr<SolverMetadata> metadata) {
            solverMetadata[typeid(SolverType).name()] = std::move(metadata);
        }

        // Get Metadata for Specific Solver Type
        template<typename SolverType>
        SolverMetadata* getSolverMetadata() {
            auto it = solverMetadata.find(typeid(SolverType).name());
            return it != solverMetadata.end() ? it->second.get() : nullptr;
        }

    private:
        std::map<std::string, std::unique_ptr<SolverMetadata>> solverMetadata;
    };

    // Numeric Method Factory
    class NumericFactory {
    public:
        // Flux Reconstruction Factory Method
        template<typename SolverType>
        static std::shared_ptr<FluxReconstruction<SolverType>> 
        createReconstructionMethod(
            FluxReconstructionType method = FluxReconstructionType::WENO,
            WenoVariant variant = WenoVariant::WENOJS5);

        // Flux Method Factory Method
        template<typename SolverType>
        static std::shared_ptr<FluxMethod<SolverType>> 
        createFluxMethod(
            FluxScheme scheme = FluxScheme::LocalLaxFriedrichs);

        // Time Stepping Scheme Factory Method
        template<typename SolverType>
        static std::shared_ptr<TimeSteppingScheme<SolverType>> 
        createTimeSteppingScheme(
            TimeSteppingSchemeType scheme = TimeSteppingSchemeType::RK3);

        // Variable Accessor Factory Method
        template<typename SolverType>
        static std::shared_ptr<GenericVariableAccessor>
        createVariableAccessor(
            ReconstructionMode mode = ReconstructionMode::Primitive,
            int numGhostCells = 2);

        // String Conversion Utilities
        static std::string toString(FluxReconstructionType method);
        static std::string toString(FluxScheme scheme);
        static std::string toString(TimeSteppingSchemeType scheme);
        static std::string toString(ReconstructionMode mode);

        // Enum Parsing Utilities
        static FluxReconstructionType parseFluxReconstruction(const std::string& str);
        static FluxScheme parseFluxScheme(const std::string& str);
        static TimeSteppingSchemeType parseTimeSteppingScheme(const std::string& str);
        static ReconstructionMode parseReconstructionMode(const std::string& str);
    };
}

#endif // NUMERIC_NUMERIC_FACTORY_H
