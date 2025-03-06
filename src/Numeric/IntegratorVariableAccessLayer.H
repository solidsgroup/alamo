// IntegratorVariableAccessLayer.H
#ifndef NUMERIC_VARIABLE_ACCESS_LAYER_H
#define NUMERIC_VARIABLE_ACCESS_LAYER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <vector>
#include <map>
#include <memory>
#include "Set/Set.H"
#include "Util/Util.H"
#include "Numeric/Stencil.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/SolverCapabilities.H"

namespace Integrator {
    class ScimitarX;
}

namespace Numeric {

class GenericVariableAccessor {
public:
    virtual ~GenericVariableAccessor() = default;

    // Enhanced method to get solver capabilities
    virtual std::shared_ptr<SolverCapabilities> getSolverCapabilities() const = 0;

    // Create working buffer with more flexible configuration
    virtual amrex::MultiFab CreateWorkingBuffer(
        const amrex::BoxArray& baseGrids, 
        const amrex::DistributionMapping& dm, 
        int num_components, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const = 0;

    // Enhanced variable transformation methods
    virtual void TransformVariables(
        int lev, 
        void* solver, 
        amrex::MultiFab& working_buffer, 
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult,
        int direction = 0
    ) const = 0;

    virtual void ReverseTransform(
        int lev, 
        void* solver, 
        amrex::MultiFab& QL_stencil,
        amrex::MultiFab& QR_stencil, 
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult,
        int direction = 0
    ) const = 0;

    // More comprehensive ghost cell management
    virtual int getRequiredGhostCells(
        ReconstructionMode mode, 
        FluxReconstructionType reconstructionType
    ) const = 0;

    // Current reconstruction mode accessor
    virtual ReconstructionMode getCurrentReconstructionMode() const = 0;
};

// Specific implementation for Compressible Euler flows
namespace CompressibleEuler {

class CompressibleEulerCapabilities : public SolverCapabilities {
public:
    std::string getIdentifier() const override { 
        return "CompressibleEuler"; 
    }

    std::string getDescription() const override { 
        return "Compressible Euler Equations Solver"; 
    }

    MethodSupport supportsFluxReconstruction(FluxReconstructionType method) const override;
    MethodSupport supportsFluxScheme(FluxScheme scheme) const override;
    MethodSupport supportsTimeSteppingScheme(TimeSteppingSchemeType scheme) const override;
    MethodSupport supportsReconstructionMode(ReconstructionMode mode) const override;
    MethodSupport supportsWenoVariant(WenoVariant variant) const override;
    MethodValidationResult validateMethodCombination(
        FluxReconstructionType fluxReconstruction,
        FluxScheme fluxScheme,
        TimeSteppingSchemeType timeSteppingScheme,
        ReconstructionMode reconstructionMode,
        WenoVariant wenoVariant = WenoVariant::WENOJS5
    ) const override;
    DefaultConfiguration getDefaultConfiguration() const override;
    std::shared_ptr<GenericVariableAccessor> createVariableAccessor(
        ReconstructionMode mode,
        int numGhostCells = 2
    ) const override;
};

class CompressibleEulerVariableAccessor : public GenericVariableAccessor {
public:
    CompressibleEulerVariableAccessor(
        int total_ghosts = 2,
        ReconstructionMode mode = ReconstructionMode::Primitive
    );

    // Implement the new virtual methods from base class
    std::shared_ptr<SolverCapabilities> getSolverCapabilities() const override;

    amrex::MultiFab CreateWorkingBuffer(
        const amrex::BoxArray& baseGrids, 
        const amrex::DistributionMapping& dm, 
        int num_components, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    void TransformVariables(
        int lev, 
        void* solver_void,
        amrex::MultiFab& working_buffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult,
        int direction = 0
    ) const override;

    void ReverseTransform(
        int lev,
        void* solver_void,
        amrex::MultiFab& QL_stencil,
        amrex::MultiFab& QR_stencil, 
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult,
        int direction = 0
    ) const override;

    int getRequiredGhostCells(
        ReconstructionMode mode, 
        FluxReconstructionType reconstructionType
    ) const override;

    ReconstructionMode getCurrentReconstructionMode() const override {
        return current_mode;
    }

private:
    // Private helper methods
    void CopyPrimitiveVariablesToWorkingBuffer(
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& working_buffer
    ) const;

    void CopyConservativeVariablesToWorkingBuffer(
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& working_buffer
    ) const;

    void ToCharacteristic(
        int direction, 
        int lev, 
        const Integrator::ScimitarX* solver, 
        amrex::MultiFab& input_mf
    ) const;

    void FromCharacteristic(
        int direction, 
        int lev, 
        const Integrator::ScimitarX* solver,
        amrex::MultiFab& QL_stencil, 
        amrex::MultiFab& QR_stencil
    ) const;

    // Companion method for matrix transformations
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::Matrix ComputeRightEigenvectorMatrix(
        const Set::Vector& W, 
        int direction
    );

    ReconstructionMode current_mode;
    int total_ghost_cells;
};

} // namespace CompressibleEuler

// Placeholder namespaces for other solver-specific variable accessors
namespace ElastoPlastic {
    class ElastoPlasticVariableAccessor;
}

namespace FiveEquationModel {
    class FiveEquationVariableAccessor;
}

} // namespace Numeric

#endif // NUMERIC_VARIABLE_ACCESS_LAYER_H
