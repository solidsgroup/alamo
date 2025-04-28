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

    // Copy of getVariableIndex struct from ScimitarX_Util.H
    struct VariableIndices {
        std::map<int, int> variableIndexMap; // Mapping from variable enums to indices
        int NVAR_MAX = 0;                    // Total number of variables
        int DENS = -1;                       // Index for Density
        int UVEL = -1;                       // Index for U-velocity
        int VVEL = -1;                       // Index for V-velocity
        int WVEL = -1;                       // Index for W-velocity
        int IE = -1;                         // Index for Internal Energy
        
        // ElastoPlastic specific variables
        int SXX = -1;                    // Index for Deviatoric stress XX
        int SYY = -1;                    // Index for Deviatoric stress YY
        int SZZ = -1;                    // Index for Deviatoric stress ZZ
        int SXY = -1;                    // Index for Deviatoric stress XY
        int SXZ = -1;                    // Index for Deviatoric stress XZ
        int SYZ = -1;                    // Index for Deviatoric stress YZ
        int EPSBAR = -1;                 // Index for Effective plastic strain
        int IE_ELASTIC = -1;             // Index for Elastic Internal Energy  
    };
    
    // Add a static instance of the struct (similar to ScimitarX's static member)
    static VariableIndices variableIndices;
    
    // Add method to initialize indices
    virtual void initializeIndices(const VariableIndices& indices) {
        // Copy only the fields that have valid values (not -1)
        if (indices.NVAR_MAX > 0) variableIndices.NVAR_MAX = indices.NVAR_MAX;
        if (indices.DENS >= 0) variableIndices.DENS = indices.DENS;
        if (indices.UVEL >= 0) variableIndices.UVEL = indices.UVEL;
        if (indices.VVEL >= 0) variableIndices.VVEL = indices.VVEL;
        if (indices.WVEL >= 0) variableIndices.WVEL = indices.WVEL;
        if (indices.IE >= 0) variableIndices.IE = indices.IE;
        
        // ElastoPlastic specific variables - only copy if they're valid
        if (indices.SXX >= 0) variableIndices.SXX = indices.SXX;
        if (indices.SYY >= 0) variableIndices.SYY = indices.SYY;
        if (indices.SZZ >= 0) variableIndices.SZZ = indices.SZZ;
        if (indices.SXY >= 0) variableIndices.SXY = indices.SXY;
        if (indices.SXZ >= 0) variableIndices.SXZ = indices.SXZ;
        if (indices.SYZ >= 0) variableIndices.SYZ = indices.SYZ;
        if (indices.EPSBAR >= 0) variableIndices.EPSBAR = indices.EPSBAR;
        if (indices.IE_ELASTIC >= 0) variableIndices.IE_ELASTIC = indices.IE_ELASTIC;
        
        // Copy only the entries that exist in the provided map
        for (const auto& [varEnum, index] : indices.variableIndexMap) {
            variableIndices.variableIndexMap[varEnum] = index;
        }
    }
    
    // Get accessor for the indices
    virtual const VariableIndices& getIndices() const {
        return variableIndices;
    }

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
    virtual void CopyVariables(
        int direction, 
        int lev, 
        void* solver, 
        amrex::MultiFab& VariableBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const = 0;

    // Enhanced variable transformation methods
    virtual void CopyFluxes(
        int direction,
        int lev, 
        void* solver, 
        amrex::MultiFab& CellFluxBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const = 0;

    // Populate Average States
    virtual void PopulateAverageStates(
        int direction,
        int lev, 
        void* solver,
        amrex::MultiFab& AverageStateBuffer
    ) const = 0;

    // Transforming Stencils to Characteristic Space
    virtual Set::MultiMatrix TransformStencilToCharacteristic(
        int i, int j, int k,    
        const Set::MultiMatrix& stencil_matrix,
        int direction,
        const Set::MultiVector& avg_state
    ) const = 0;
    
    //Reverse Transforming From Characteristic Space
    virtual Set::MultiVector TransformFromCharacteristic(
        const Set::MultiVector& reconstructed_value,
        int direction,
        const Set::MultiVector& avg_state
    ) const = 0;
    
    // Store directional flux to appropriate arrays
    virtual void StoreDirectionalFlux(
        int direction,
        int lev, 
        void* solver, 
        amrex::MultiFab& SummedFlux
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

    void CopyVariables(
        int direction, 
        int lev, 
        void* solver_void,
        amrex::MultiFab& VariableBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    void CopyFluxes(
        int direction,
        int lev, 
        void* solver_void,
        amrex::MultiFab& CellFluxBuffer,
        ReconstructionMode mode, 
        const SolverCapabilities::MethodValidationResult& validationResult
    ) const override;

    // Populate Average States
    void PopulateAverageStates(
        int direction,
        int lev, 
        void* solver,
        amrex::MultiFab& AverageStateBuffer
    ) const override;

    // Transforming Stencils to Characteristic Space
    Set::MultiMatrix TransformStencilToCharacteristic(
        int i, int j, int k,    
        const Set::MultiMatrix& stencil_matrix,
        int direction,
        const Set::MultiVector& avg_state
    ) const override;
    
    //Reverse Transforming From Characteristic Space
    Set::MultiVector TransformFromCharacteristic(
        const Set::MultiVector& reconstructed_value,
        int direction,
        const Set::MultiVector& avg_state
    ) const override;

    void StoreDirectionalFlux(
        int direction,
        int lev, 
        void* solver_void,
        amrex::MultiFab& SummedFlux
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

    // Companion method for matrix transformations
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::MultiMatrix ComputeRightEigenvectorMatrix(
        const Set::MultiVector& W, 
        int direction, int num_components
    );

    // Companion method for matrix transformations
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::MultiMatrix ComputeLeftEigenvectorMatrix(
        int i, int j, int k,    
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
