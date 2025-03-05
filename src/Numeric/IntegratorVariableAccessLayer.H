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
#include "Numeric/NumericFactory.H"

namespace Integrator {
    class ScimitarX;
}

namespace Numeric {

// Base class for generic variable accessor
class GenericVariableAccessor {
public:
    virtual ~GenericVariableAccessor() = default;

    virtual amrex::MultiFab CreateWorkingBuffer(
        const amrex::BoxArray& baseGrids, 
        const amrex::DistributionMapping& dm, 
        int num_components, 
        int num_ghost_cells = 0
    ) const = 0;

    virtual void TransformVariables(
        int lev, 
        void* solver, 
        amrex::MultiFab& working_buffer, 
        Numeric::ReconstructionMode mode, 
        int direction = 0
    ) const = 0;

    virtual void ReverseTransform(
        int lev, 
        void* solver, 
        amrex::MultiFab& QL_stencil,
        amrex::MultiFab& QR_stencil, 
        Numeric::ReconstructionMode mode, 
        int direction = 0
    ) const = 0;

    virtual int GetTotalGhostCells() const = 0;

  // Add this declaration
    virtual Numeric::ReconstructionMode GetReconstructionMode() const = 0;
protected:

    Numeric::ReconstructionMode current_mode;
    int total_ghost_cells;
};

// Placeholder namespaces for other solver-specific variable accessors
namespace ElastoPlastic {
    class ElastoPlasticVariableAccessor;
}

namespace FiveEquationModel {
    class FiveEquationVariableAccessor;
}

namespace CompressibleEuler {

class CompressibleEulerVariableAccessor : public GenericVariableAccessor {
public:
        
    // Constants specific to Compressible Euler solver
    static constexpr int NumComponents = AMREX_SPACEDIM + 2;  // ρ, ρu, ρv, [ρw], E
    static constexpr Set::Scalar Gamma = 1.4;  // Ratio of specific heats

    CompressibleEulerVariableAccessor(
        int total_ghosts = 2,
        Numeric::ReconstructionMode mode = Numeric::ReconstructionMode::Primitive
    ) {
        total_ghost_cells = total_ghosts;
        current_mode = mode;
    }

    // Implement GetReconstructionMode method
    Numeric::ReconstructionMode GetReconstructionMode() const override {
        return current_mode;
    }

    amrex::MultiFab CreateWorkingBuffer(
        const amrex::BoxArray& baseGrids, 
        const amrex::DistributionMapping& dm, 
        int num_components, 
        int num_ghost_cells = 0
    ) const override {
        return amrex::MultiFab(
            baseGrids, 
            dm, 
            num_components, 
            std::max(total_ghost_cells, num_ghost_cells)
        );
    }

    void TransformVariables(
        int lev,
        void* solver_void,
        amrex::MultiFab& working_buffer,
        Numeric::ReconstructionMode mode,
        int direction = 0
    ) const override {
        Integrator::ScimitarX* solver = static_cast<Integrator::ScimitarX*>(solver_void);

        switch (mode) {
            case Numeric::ReconstructionMode::Primitive:
                CopyPrimitiveVariablesToWorkingBuffer(lev, solver, working_buffer);
                break;
            case Numeric::ReconstructionMode::Conservative:
                CopyConservativeVariablesToWorkingBuffer(lev, solver, working_buffer);
                break;
            case Numeric::ReconstructionMode::Characteristic:
                CopyConservativeVariablesToWorkingBuffer(lev, solver, working_buffer);
                ToCharacteristic(direction, lev, solver, working_buffer);
                break;
            default:
                Util::Abort(INFO, "Unknown reconstruction mode");
        }
    }

    void ReverseTransform(
        int lev,
        void* solver_void,
        amrex::MultiFab& QL_stencil,
        amrex::MultiFab& QR_stencil,
        Numeric::ReconstructionMode mode,
        int direction = 0
    ) const override {
        if (mode == Numeric::ReconstructionMode::Characteristic) {
            Integrator::ScimitarX* solver = static_cast<Integrator::ScimitarX*>(solver_void);
            FromCharacteristic(direction, lev, solver, QL_stencil, QR_stencil);
        }
    }

    int GetTotalGhostCells() const override {
        return total_ghost_cells;
    }

private:
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::Matrix ComputeRightEigenvectorMatrix(
        const Set::Vector& W, 
        int direction
    );
};

/*// Forward declarations for CompressibleEuler's Inline Methods
inline void CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToWorkingBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& working_buffer
) const;

inline void CompressibleEulerVariableAccessor::CopyConservativeVariablesToWorkingBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& working_buffer
) const;

inline void CompressibleEulerVariableAccessor::ToCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver, 
    amrex::MultiFab& input_mf
) const;

inline void CompressibleEulerVariableAccessor::FromCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& QL_stencil, 
    amrex::MultiFab& QR_stencil
) const;*/

}

/*// Placeholder Variable Accessors for other Solver Types
namespace ElastoPlastic {
    class ElastoPlasticVariableAccessor : public GenericVariableAccessor {
    public:
        // Implement base class methods with placeholders
        ElastoPlasticVariableAccessor(
            int total_ghosts = 2,
            ReconstructionMode mode = ReconstructionMode::Primitive
        ) {
            total_ghost_cells = total_ghosts;
            current_mode = mode;
        }

        amrex::MultiFab CreateWorkingBuffer(
            const amrex::BoxArray& baseGrids, 
            const amrex::DistributionMapping& dm, 
            int num_components, 
            int num_ghost_cells = 0
        ) const override {
            Util::Abort(INFO, "ElastoPlastic CreateWorkingBuffer not implemented");
            return amrex::MultiFab();
        }

        void TransformVariables(
            int lev, 
            void* solver, 
            amrex::MultiFab& working_buffer, 
            ReconstructionMode mode, 
            int direction = 0
        ) const override {
            Util::Abort(INFO, "ElastoPlastic TransformVariables not implemented");
        }

        void ReverseTransform(
            int lev, 
            void* solver, 
            amrex::MultiFab& QL_stencil,
            amrex::MultiFab& QR_stencil, 
            ReconstructionMode mode, 
            int direction = 0
        ) const override {
            Util::Abort(INFO, "ElastoPlastic ReverseTransform not implemented");
        }

        int GetTotalGhostCells() const override {
            return total_ghost_cells;
        }

        ReconstructionMode GetReconstructionMode() const override {
            return current_mode;
        }
    };
}

namespace FiveEquationModel {
    class FiveEquationVariableAccessor : public GenericVariableAccessor {
    public:
        // Implement base class methods with placeholders
        FiveEquationVariableAccessor(
            int total_ghosts = 2,
            ReconstructionMode mode = ReconstructionMode::Primitive
        ) {
            total_ghost_cells = total_ghosts;
            current_mode = mode;
        }

        amrex::MultiFab CreateWorkingBuffer(
            const amrex::BoxArray& baseGrids, 
            const amrex::DistributionMapping& dm, 
            int num_components, 
            int num_ghost_cells = 0
        ) const override {
            Util::Abort(INFO, "FiveEquationModel CreateWorkingBuffer not implemented");
            return amrex::MultiFab();
        }

        void TransformVariables(
            int lev, 
            void* solver, 
            amrex::MultiFab& working_buffer, 
            ReconstructionMode mode, 
            int direction = 0
        ) const override {
            Util::Abort(INFO, "FiveEquationModel TransformVariables not implemented");
        }

        void ReverseTransform(
            int lev, 
            void* solver, 
            amrex::MultiFab& QL_stencil,
            amrex::MultiFab& QR_stencil, 
            ReconstructionMode mode, 
            int direction = 0
        ) const override {
            Util::Abort(INFO, "FiveEquationModel ReverseTransform not implemented");
        }

        int GetTotalGhostCells() const override {
            return total_ghost_cells;
        }

        ReconstructionMode GetReconstructionMode() const override {
            return current_mode;
        }
    };
}*/

} // namespace Numeric

#endif // NUMERIC_VARIABLE_ACCESS_LAYER_H
