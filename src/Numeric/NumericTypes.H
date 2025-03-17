// NumericTypes.H
#ifndef NUMERIC_TYPES_H
#define NUMERIC_TYPES_H

#include <vector>
#include <string>
#include <memory>

// Explicit specialization for Integrator::ScimitarX
 namespace Integrator {
        class ScimitarX;
    }

namespace Numeric {
    // Forward declarations for template classes
    template <typename IntegratorType> class FluxReconstruction;
    template <typename IntegratorType> class FluxMethod;
    template <typename IntegratorType> class TimeSteppingScheme;
    template <typename IntegratorType> class EulerForwardScheme;
    template <typename IntegratorType> class RK3Scheme;
    template <typename IntegratorType> class FirstOrderReconstruction;
    template <typename IntegratorType, typename VariantTraits> class WENOReconstruction;
    template <typename IntegratorType> class LocalLaxFriedrichsMethod;
    
    // Forward declaration of non-template classes
    class GenericVariableAccessor;
    class SolverCapabilities;
    
    // Enum Classes
    enum class FluxReconstructionType {
        FirstOrder,    // Basic first-order upwind reconstruction
        WENO          // Weighted Essentially Non-Oscillatory
    };

    enum class WenoVariant {
        WENOJS5,       // WENO Jiang-Shu 5th Order
        WENOZ5,        // WENO Z 5th Order
        WENOJS3        // WENO JS 3rd Order
    };

    enum class FluxScheme {
        LocalLaxFriedrichs,  // Local Lax-Friedrichs (Rusanov) flux
        HLLC                // Harten-Lax-van Leer Contact flux
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
    
    
    // Template specialization forward declaration for getSolverIdentifier
    template<typename SolverType>
    std::string getSolverIdentifier();


    template<>
    std::string getSolverIdentifier<Integrator::ScimitarX>();
}

#endif // NUMERIC_TYPES_H
