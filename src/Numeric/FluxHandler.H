#ifndef NUMERIC_FLUXHANDLER_H
#define NUMERIC_FLUXHANDLER_H

#include "Set/Set.H"
#include "Integrator/ScimitarX.H"
#include <memory>

namespace Numeric {

// Template-based base class for flux reconstruction
template <typename IntegratorDerivedClass>
class FluxReconstruction {
public:
    virtual ~FluxReconstruction() = default;

    /**
     * @brief Perform flux reconstruction in the specified direction.
     * 
     * @param direction The spatial direction (0 for x, 1 for y, 2 for z).
     * @param lev The AMR level.
     * @param solver A pointer to the derived solver class.
     */
    virtual void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver) const = 0;
};

// Template-based base class for flux methods (e.g., Local Lax-Friedrichs, HLLC)
template <typename IntegratorDerivedClass>
class FluxMethod {
public:
    virtual ~FluxMethod() = default;

    /**
     * @brief Compute the numerical flux in the specified direction.
     * 
     * @param direction The spatial direction (0 for x, 1 for y, 2 for z).
     * @param lev The AMR level.
     * @param solver A pointer to the derived solver class.
     */
    virtual void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver) const = 0;
};

// FluxHandler class to manage both reconstruction and flux computation
template <typename IntegratorDerivedClass>
class FluxHandler {
public:
    FluxHandler() {}

    /**
     * @brief Set the flux reconstruction method.
     * 
     * @param recon Shared pointer to a flux reconstruction object.
     */
    void SetReconstruction(std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> recon) {
        reconstruction = std::move(recon);
    }

    /**
     * @brief Set the flux computation method.
     * 
     * @param method Shared pointer to a flux method object.
     */
    void SetFluxMethod(std::shared_ptr<FluxMethod<IntegratorDerivedClass>> method) {
        flux_method = std::move(method);
    }

    /**
     * @brief Construct the fluxes in all spatial directions.
     * 
     * @param lev The AMR level.
     * @param solver A pointer to the derived solver class.
     */
    void ConstructFluxes(int lev, IntegratorDerivedClass* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
            return;
        }

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            reconstruction->ReconstructFlux(dir, lev, solver);
            flux_method->ComputeFlux(dir, lev, solver);
        }
    }

private:
    std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> reconstruction;  // Flux reconstruction method
    std::shared_ptr<FluxMethod<IntegratorDerivedClass>> flux_method;             // Flux computation method
};

// Concrete implementation of first-order flux reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver) const override {
        Util::Message(INFO, "Performing first-order flux reconstruction.");
        // Implement first-order reconstruction logic here.
    }
};

// Concrete implementation of Local Lax-Friedrichs flux method
template <typename IntegratorDerivedClass>
class LocalLaxFriedrichsMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver) const override {
        Util::Message(INFO, "Computing Local Lax-Friedrichs flux.");
        // Implement LLF flux computation here.
    }
};

// Concrete implementation of HLLC flux method
template <typename IntegratorDerivedClass>
class HLLCFluxMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver) const override {
        Util::Message(INFO, "Computing HLLC flux.");
        // Implement HLLC flux computation here.
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H

