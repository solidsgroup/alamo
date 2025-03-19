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
        const int num_components,
        const int ghost_cells         
    ) const = 0;

    // Enhanced variable transformation methods
    virtual void TransformVariables(
        int direction, 
        int lev, 
        void* solver, 
        amrex::MultiFab& VariableBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const = 0;

    // Enhanced variable transformation methods
    virtual void TransformFluxes(
        int direction,
        int lev, 
        void* solver, 
        amrex::MultiFab& CellFluxBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const = 0;

    virtual void ReverseTransformFluxes(
        int direction,
        int lev, 
        void* solver, 
        amrex::MultiFab& SummedFlux,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
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
        const int num_components,
        const int ghost_cells 
    ) const override;

    void TransformVariables(
        int direction, 
        int lev, 
        void* solver_void,
        amrex::MultiFab& VariableBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    void TransformFluxes(
        int direction,
        int lev, 
        void* solver_void,
        amrex::MultiFab& CellFluxBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    void ReverseTransformFluxes(
        int direction,
        int lev,
        void* solver_void,
        amrex::MultiFab& SummedFlux,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    int getRequiredGhostCells(
        ReconstructionMode mode, 
        FluxReconstructionType reconstructionType
    ) const override;

    ReconstructionMode getCurrentReconstructionMode() const override {
        return current_mode;
    }

    // Computing Roe Averages
    static Set::MultiVector ComputeRoeAverages(
        const Set::MultiVector& WL, 
        const Set::MultiVector& WR,
        int num_components
    );

private:
    // Private helper methods
    void CopyPrimitiveVariablesToVariableBuffer(
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& VariableBuffer
    ) const;

    void CopyConservativeVariablesToVariableBuffer(
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& VariableBuffer
    ) const;

    void CopyPrimitiveFluxesToCellFluxBuffer(
        int direction,    
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& CellFluxBuffer
    ) const;

    void CopyConservativeFluxesToCellFluxBuffer(
        int direction,    
        int lev,
        Integrator::ScimitarX* solver,
        amrex::MultiFab& CellFluxBuffer
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
        amrex::MultiFab& SummedFlux 
    ) const;

    void FromConservative(
        int direction, 
        int lev, 
        const Integrator::ScimitarX* solver,
        amrex::MultiFab& SummedFlux 
    ) const;

    // Companion method for matrix transformations
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::MultiMatrix ComputeRightEigenvectorMatrix(
        const Set::MultiVector& W, 
        int direction, int num_components
    );

    // Companion method for matrix transformations
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::MultiMatrix ComputeLeftEigenvectorMatrix(
        const Set::MultiVector& W, 
        int direction, int num_components
    );


    ReconstructionMode current_mode;
    int total_ghost_cells;

    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

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
