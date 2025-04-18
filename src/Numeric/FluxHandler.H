#ifndef NUMERIC_FLUXHANDLER_H
#define NUMERIC_FLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/NarrowBandLevelset_Backup.H"
#include "Integrator/NarrowBandLevelset.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H"  // Include fluid model helper functions
#include "Numeric/Stencil.H"                                
#include <memory>
#include <cmath>
#include <algorithm>

// Include foward declare for ScimitarX.H
namespace Integrator {
class ScimitarX;
class NarrowBandLevelset_Backup;
class NarrowBandLevelset;
}

namespace Numeric {

/*// Clamp all indices
AMREX_GPU_DEVICE
inline void ClampIndices(int indices[3], const int lower_bounds[3], const int upper_bounds[3]) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        indices[d] = amrex::max(lower_bounds[d], amrex::min(indices[d], upper_bounds[d]));
    }
}

// Clamp a single index (optimized)
AMREX_GPU_DEVICE
inline void ShiftAndClampIndex(int& index, int offset, int lower_bound, int upper_bound) {
    index = amrex::max(lower_bound, amrex::min(index + offset, upper_bound));
}*/

// ClampIndices: Clamps indices to valid bounds (branch-free logic)
AMREX_GPU_DEVICE
inline void ClampIndices(int indices[3], const int lower_bounds[3], const int upper_bounds[3]) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        indices[d] = lower_bounds[d] + std::min(std::max(indices[d] - lower_bounds[d], 0), upper_bounds[d] - lower_bounds[d]);
    }
}

// ShiftAndClampIndices: Shifts and clamps indices (branch-free logic)
AMREX_GPU_DEVICE
inline void ShiftAndClampIndices(int indices[3], const int offsets[3],
    const int lower_bounds[3], const int upper_bounds[3], int direction) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        int shifted_index = indices[d] + (d == direction ? offsets[d] : 0);
        indices[d] = lower_bounds[d] + std::min(std::max(shifted_index - lower_bounds[d], 0), upper_bounds[d] - lower_bounds[d]);
    }
}

// Base class for flux reconstruction
template <typename IntegratorDerivedClass>
class FluxReconstruction {
public:
    virtual ~FluxReconstruction() = default;

    /**
    * @brief Perform flux reconstruction for the specified direction.
    */
    virtual void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const = 0;
};

// Base class for flux method (e.g., Local Lax-Friedrichs, HLLC)
template <typename IntegratorDerivedClass>
class FluxMethod {
public:
    virtual ~FluxMethod() = default;

    /**
    * @brief Compute the flux for the specified direction.
    */
    virtual void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil) const = 0;
};

// NarrowBand specific flux class
template <>
class FluxMethod<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxMethod() = default;

    /**
    * @brief Compute the flux for the specified direction.
    */
    virtual void ComputeFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil, amrex::MultiFab& flux) const = 0;
};

// FluxHandler class for managing reconstruction and flux computation
template <typename IntegratorDerivedClass>
class FluxHandler {
public:
    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

    FluxHandler() {}

    void SetReconstruction(std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> recon) {
        reconstruction = std::move(recon);
    }

    void SetFluxMethod(std::shared_ptr<FluxMethod<IntegratorDerivedClass>> method) {
        flux_method = std::move(method);
    }

    void ConstructFluxes(int lev, IntegratorDerivedClass* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
        }

        //  Util::Message(INFO, "Constructing fluxes for level " + std::to_string(lev));

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            amrex::BoxArray face_ba = solver->grids[lev];
            amrex::IntVect face_type = amrex::IntVect::TheZeroVector();  // `{0, 0, 0}`
            face_type[dir] = 1;  // Make it face-centered in `dir`
            face_ba.convert(face_type);

            // Create temporary QL and QR MultiFabs at face-centered positions
            amrex::MultiFab QL_stencil(face_ba, solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);
            amrex::MultiFab QR_stencil(face_ba, solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);

            reconstruction->ReconstructFlux(dir, lev, solver, QL_stencil, QR_stencil);
            flux_method->ComputeFlux(dir, lev, solver, QL_stencil, QR_stencil);
        }

    }

private:
    std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<IntegratorDerivedClass>> flux_method;             // Flux method

};

inline void ConvertToFaceCentered(int dir, int lev, Set::Field<Set::Scalar>& flux_mf_wrapper) {
    // Obtain the BoxArray of XFlux_mf/YFlux_mf
    amrex::BoxArray flux_ba = flux_mf_wrapper[lev]->boxArray();  // Assuming we're working with level 0

    // Create a zero vector for the face_type
    amrex::IntVect face_type = amrex::IntVect::TheZeroVector();  // {0, 0, 0}
    
    // Set the corresponding direction to 1 for face-centered
    face_type[dir] = 1;  // 0: cell-centered, 1: face-centered in the 'dir' direction

    // Convert the BoxArray to face-centered for the given direction
    flux_ba.convert(face_type);

    // Define a new MultiFab with the updated BoxArray and existing distribution mapping
    const amrex::DistributionMapping& dmap = flux_mf_wrapper[lev]->DistributionMap();
    flux_mf_wrapper.Define(lev, flux_ba, dmap, flux_mf_wrapper[0]->nComp(), flux_mf_wrapper[0]->nGrow());
}

// FluxHandler class for Integrator::NarrowBandLevelset to correctly assess distribution 
// mapping for narrowband
template <>
class FluxHandler<Integrator::NarrowBandLevelset> {
public:
    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

    FluxHandler() {}

    void SetReconstruction(std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> recon) {
        reconstruction = std::move(recon);
    }

    void SetFluxMethod(std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> method) {
        flux_method = std::move(method);
    }

    void ConstructFluxes(int lev, Integrator::NarrowBandLevelset* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
        }

        // Get number of levelset fields
        int ncomp = solver->number_of_components;
        int nghosts = solver->number_of_ghost_cells;

        // Loop through all levelset fields
        for (int ils = 0; ils < ncomp; ++ils) {
            // Get proper distribution mapping for narrowband
            const amrex::BoxArray& narrowband_ba = solver->level_sets[ils].narrowband_ba;
            const amrex::DistributionMapping& narrowband_dm = solver->level_sets[ils].narrowband_dm;

            // Loop through all dimensions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Convert narrowband_ba to face-centered direction dir
                amrex::BoxArray face_ba(narrowband_ba);  // make a copy to modify
                amrex::IntVect face_type = amrex::IntVect::TheZeroVector(); // {0, 0, 0}
                face_type[dir] = 1;
                face_ba.convert(face_type);  // now face_ba is face-centered in `dir`

                // Get the corresponding full flux MultiFab and its index type
                auto& flux_mf =
                    (dir == 0) ? solver->XFlux_mf:
                    (dir == 1) ? solver->YFlux_mf:
                                 solver->ZFlux_mf;
                
                // Convert multifab index type if index type does not match correct direction
                if (face_ba.ixType() != flux_mf[lev]->boxArray().ixType()){
                    ConvertToFaceCentered(dir, lev, flux_mf);
                }      

                // Create temporary QL and QR MultiFabs at face-centered positions
                amrex::MultiFab QL_stencil(face_ba, narrowband_dm, ncomp, nghosts);
                amrex::MultiFab QR_stencil(face_ba, narrowband_dm, ncomp, nghosts);

                // Temporary flux MultiFab with same indexing
                amrex::MultiFab temp_flux(face_ba, narrowband_dm, ncomp, nghosts);

                // Perform reconstruction and flux computation
                reconstruction->ReconstructFlux(dir, lev, solver, QL_stencil, QR_stencil);
                flux_method->ComputeFlux(dir, lev, solver, QL_stencil, QR_stencil, temp_flux);
            }
        }
    }

private:
    std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> flux_method;             // Flux method
};

/*// FluxHandler class for Integrator::NarrowBandLevelset to correctly assess distribution
// mapping for narrowband
template <>
class FluxHandler<Integrator::NarrowBandLevelset> {
public:
    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

    FluxHandler() {}

    void SetReconstruction(std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> recon) {
        reconstruction = std::move(recon);
    }

    void SetFluxMethod(std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> method) {
        flux_method = std::move(method);
    }

    void ConstructFluxes(int lev, Integrator::NarrowBandLevelset* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
        }

        // Get number of levelset fields
        int ncomp = solver->number_of_components;
        int nghosts = solver->number_of_ghost_cells;

        // Loop through all levelset fields
        for (int ils=0; ils<ncomp; ils++){
            // Get proper distribution mapping for narrowband
            const amrex::BoxArray narrowband_ba = solver->level_sets[ils].narrowband_ba;
            const amrex::DistributionMapping narrowband_dm = solver->level_sets[ils].narrowband_dm;

            // Loop through all dimensions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Convert narrowband_ba to face-centered direction dir
                amrex::BoxArray face_ba(narrowband_ba);  // make a copy to modify
                amrex::IntVect face_type = amrex::IntVect::TheZeroVector(); // {0, 0, 0}
                face_type[dir] = 1;
                face_ba.convert(face_type);  // now face_ba is face-centered in `dir`

                // Create temporary QL and QR MultiFabs at face-centered positions
                amrex::MultiFab QL_stencil(face_ba, narrowband_dm, ncomp, nghosts);
                amrex::MultiFab QR_stencil(face_ba, narrowband_dm, ncomp, nghosts);

                // Get the proper index information for flux multifab
                auto const& flux_mf = (dir == 0) ? solver->XFlux_mf[lev] :
                (dir == 1) ? solver->YFlux_mf[lev] :
                solver->ZFlux_mf[lev];

                // Create temporary flux multifab to save computations and parallel copy back
                amrex::MultiFab temp_flux(face_ba, narrowband_dm, ncomp, nghosts);

                // Perform reconstruction and flux computation
                reconstruction->ReconstructFlux(dir, lev, solver, QL_stencil, QR_stencil);
                flux_method->ComputeFlux(dir, lev, solver, QL_stencil, QR_stencil, temp_flux);
            }
        }
    }

private:
    std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> flux_method;             // Flux method

};*/

// First-order upwind reconstruction base class
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    // ReconstructFlux: First-order upwind flux reconstruction with branch-free framework
    void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
        //Util::Message(INFO, "Performing first-order upwind flux reconstruction for direction: " + std::to_string(direction));
        //Util::ScimitarX_Util::Debug debug;
        //debug.SetTargetDebugLocationIndices(101, 0, 0);  // Debug target location
        //debug.SetTargetDebugLocation(0.49, 0.12, 0.0, solver->geom[lev]);  // Debug target location
        const int nghosts = solver->number_of_ghost_cells;
        // Loop through face-centered indices for flux reconstruction
        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            auto const& QL_arr = QL_stencil.array(mfi);  // Left reconstructed state
            auto const& QR_arr = QR_stencil.array(mfi);  // Right reconstructed state
            auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);  // Conservative variables

            const int num_components = solver->number_of_components;  // Number of components (e.g., rho, u, v, E)

            /*// Define velocity components based on direction
            int velocity_component[3] = {
                solver->variableIndex.UVEL,  // X-direction -> UVEL
                solver->variableIndex.VVEL,  // Y-direction -> VVEL
                solver->variableIndex.WVEL   // Z-direction -> WVEL
            };*/

            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

                // Construct index arrays with correct number of components
                int index[3];
                int left_index[3];
                int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];

#if AMREX_SPACEDIM == 2
                index[0] = iface; index[1] = jface; index[2] = kface;
                left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1; upper_bounds[2] = 0;
#else
                index[0] = iface; index[1] = jface; index[2] = kface;
                left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1; upper_bounds[2] = face_bx_with_ghosts.bigEnd(2) - 1;
#endif

                // Clamp the main index
                ClampIndices(index, lower_bounds, upper_bounds);

                // Shift and clamp for left and right indices
                int left_offset[3] = { 0, 0, 0 };
                int right_offset[3] = { 0, 0, 0 };
                left_offset[direction] = 0;  // Shift left
                right_offset[direction] = 1;  // Shift right

                ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                //debug.DebugIndicesAndBounds(iface, jface, kface, index, left_index, right_index, lower_bounds, upper_bounds, "QL and QR Construction");

                    //const Set::Scalar vel = p_arr(index[0], index[1], index[2], velocity_component[direction]);  // Normal velocity
                    //const Set::Scalar sigma = std::copysign(1.0, vel);  // 1 if vel > 0, 0 if vel <= 0

                for (int n = 0; n < num_components; ++n) {
                    // Branch-free computation of QL and QR
                    //const Set::Scalar QL = 0.5 * ((1 + sigma)* p_arr(index[0], index[1], index[2], n) + (1.0 - sigma) * p_arr(left_index[0], left_index[1], left_index[2], n));
                    //const Set::Scalar QR = 0.5 * ((1 + sigma)* p_arr(right_index[0], right_index[1], right_index[2], n) + (1.0 - sigma)*p_arr(index[0], index[1], index[2],n));  // Since we are using first-order upwinding  

                    const Set::Scalar QL = p_arr(left_index[0], left_index[1], left_index[2], n);
                    const Set::Scalar QR = p_arr(right_index[0], right_index[1], right_index[2], n);  // Since we are using first-order upwinding  
                    // Assign values to the stencil arrays
                    QL_arr(index[0], index[1], index[2], n) = QL;
                    QR_arr(index[0], index[1], index[2], n) = QR;

                    // Debugging output
                    //debug.DebugValuesIfTarget(iface, jface, kface, QL, "QL_arr", "Reconstructed Left State", false, n);
                    //debug.DebugValuesIfTarget(iface, jface, kface, QR, "QR_arr", "Reconstructed Right State", false, n);
                }
            });
        }

    }

};

inline void FillBoundaryField(Set::Field<Set::Scalar>& field, int level) {
    if (field.size() > level && field[level]) {
        field[level]->FillBoundary();
    }
}

// First-order specific to NarrowBandLevelset integrator class
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset> : public FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    void ReconstructFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver,
        amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
        // Define multifab constants 
        const int nghosts = solver->number_of_ghost_cells;
        const int ncomp = solver->number_of_components;

        // Define OutsideNarrowband
        const int OutsideNarrowband = solver->GetOutsideNarrowBandPos();

        // Loop through all the levelset fields
        for (int ils = 0; ils < ncomp; ils++) {
            // Get the narrowband box array and distibution mapping for temporary multifabs
            amrex::BoxArray narrowband_ba = solver->level_sets[ils].narrowband_ba;
            amrex::DistributionMapping narrowband_dm = solver->level_sets[ils].narrowband_dm;

            // Create temporary multifabs
            amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, nghosts);
            amrex::MultiFab narrowband_LS(narrowband_ba, narrowband_dm, 1, nghosts);
            amrex::MultiFab narrowband_Vel(narrowband_ba, narrowband_dm, AMREX_SPACEDIM, nghosts);

            // ParallelCopy
            auto const& Tube = solver->level_sets[ils].Tube_imf;
            auto const& velocity = solver->level_sets[ils].velocity_mf;
            narrowband.ParallelCopy(*Tube, 0, 0, 1, nghosts, nghosts);
            narrowband_LS.ParallelCopy(*solver->ls_mf[lev], ils, 0, 1, nghosts, nghosts);
            narrowband_Vel.ParallelCopy(*velocity[lev], 0, 0, AMREX_SPACEDIM, nghosts, nghosts);

            // Fill Boundary for ghost cells
            narrowband.FillBoundary();
            narrowband_LS.FillBoundary();
            narrowband_Vel.FillBoundary();

            // Loop through the face centered stencil
            for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
                // Extend the box to access ghost cells
                const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

                // Compute valid index bounds using compiler
                int lower_bounds[3] = { 0, 0, 0 };
                int upper_bounds[3] = { 0, 0, 0 };

#if AMREX_SPACEDIM >= 1
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0);
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
#endif

#if AMREX_SPACEDIM >= 2
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1);
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
#endif

#if AMREX_SPACEDIM == 3
                lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
                upper_bounds[2] = face_bx_with_ghosts.bigEnd(2) - 1;
#endif

                // Define multifab arrays
                auto const& QL_arr = QL_stencil.array(mfi);
                auto const& QR_arr = QR_stencil.array(mfi);
                auto const& nb_arr = narrowband.array(mfi);
                auto const& ls = narrowband_LS.array(mfi);
                auto const& vel_arr = narrowband_Vel.array(mfi);

                // Loop through all face cell indices
                amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                    // Construct cell-centered indices
                    int cell_centered_index[3] = { iface, jface, kface };
                    int left_cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_index[3] = { iface, jface, kface };
                    int left_cell_centered_offset[3] = { 0, 0, 0 };
                    int right_cell_centered_offset[3] = { 0, 0, 0 };

                    // Define offsets (modular for future higher-order schemes)
                    left_cell_centered_offset[direction] = -1;  // Left offset for first-order upwind scheme
                    right_cell_centered_offset[direction] = 1;  // Right offset for first-order upwind scheme

                    // Clamp cell centered indices
                    ClampIndices(cell_centered_index, lower_bounds, upper_bounds);
                    ShiftAndClampIndices(left_cell_centered_index, left_cell_centered_offset, lower_bounds, upper_bounds, direction);
                    ShiftAndClampIndices(right_cell_centered_index, right_cell_centered_offset, lower_bounds, upper_bounds, direction);

                    // Create amrex::IntVect to store coords
                    const amrex::IntVect cell_center_IV(AMREX_D_DECL(cell_centered_index[0], 
                                                                     cell_centered_index[1], 
                                                                     cell_centered_index[2])); 
                    const amrex::IntVect left_cell_center_IV(AMREX_D_DECL(left_cell_centered_index[0], 
                                                                          left_cell_centered_index[1], 
                                                                          left_cell_centered_index[2]));
                    const amrex::IntVect right_cell_center_IV(AMREX_D_DECL(right_cell_centered_index[0], 
                                                                           right_cell_centered_index[1], 
                                                                           right_cell_centered_index[2]));

                    // Exit early if not narrowband cell
                    int nb_val = nb_arr(cell_center_IV);
                    //amrex::Print() << "nb_arr" << cell_center_IV << " : " << nb_arr(cell_center_IV) << std::endl;
                    if (std::abs(nb_val) == OutsideNarrowband) return;

                    // Clamp face centered indices
                    int face_centered_index[3] = { iface, jface, kface };
                    ClampIndices(face_centered_index, lower_bounds, upper_bounds);
                    const amrex::IntVect face_centered_IV(AMREX_D_DECL(face_centered_index[0], 
                                                                       face_centered_index[1], 
                                                                       face_centered_index[2]));

                    // Compute velocity and its sign (branch-free)
                    Set::Scalar vel = vel_arr(cell_center_IV, direction);
                    int sigma = vel >= 0 ? 1 : 0;

                    // First-order upwind reconstruction
                    Set::Scalar QL = sigma * ls(cell_center_IV) + (1.0 - sigma) * ls(left_cell_center_IV);
                    Set::Scalar QR = sigma * ls(right_cell_center_IV) + (1.0 - sigma) * ls(cell_center_IV);

                    // Store results
                    QL_arr(face_centered_IV, ils) = QL;
                    QR_arr(face_centered_IV, ils) = QR;
                });
            }
        }
    }
};

// First-order specific to NarrowBandLevelset integrator class
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset_Backup> : public FluxReconstruction<Integrator::NarrowBandLevelset_Backup> {
public:
    void ReconstructFlux(int direction, int lev, Integrator::NarrowBandLevelset_Backup* solver,
        amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
        const int nghosts = solver->number_of_ghost_cells;
        //printf("direction: %d\n", direction);

        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);
            auto const& ls = solver->ls_mf.Patch(lev, mfi);            // Level-set field
            auto const& vel_arr = solver->velocity_mf.Patch(lev, mfi); // Velocity field

            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays with correct number of components
                int index[3] = { iface, jface, kface };
                int left_index[3] = { iface, jface, kface };
                int right_index[3] = { iface, jface, kface };

                // Compute valid index bounds
                int lower_bounds[3] = {
                    face_bx_with_ghosts.smallEnd(0),
                    (AMREX_SPACEDIM > 1) ? face_bx_with_ghosts.smallEnd(1) : 0,
                    (AMREX_SPACEDIM > 2) ? face_bx_with_ghosts.smallEnd(2) : 0
                };

                int upper_bounds[3] = {
                    face_bx_with_ghosts.bigEnd(0) - 1,
                    (AMREX_SPACEDIM > 1) ? face_bx_with_ghosts.bigEnd(1) - 1 : 0,
                    (AMREX_SPACEDIM > 2) ? face_bx_with_ghosts.bigEnd(2) - 1 : 0
                };

                // Clamp the main index
                ClampIndices(index, lower_bounds, upper_bounds);

                // Define offsets (modular for future higher-order schemes)
                int left_offset[3] = { 0, 0, 0 };
                int right_offset[3] = { 0, 0, 0 };
                left_offset[direction] = -1;  // Left offset for first-order upwind scheme
                right_offset[direction] = 1;  // Right offset for first-order upwind scheme

                ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                // Compute velocity and its sign (branch-free)
                //Set::Scalar vel = vel_arr(index[0], index[1], index[2])[direction];
                Set::Scalar vel = vel_arr(index[0], index[1], index[2], direction);
                int sigma = vel >= 0 ? 1 : 0;

                // First-order upwind reconstruction
                Set::Scalar QL = sigma * ls(index[0], index[1], index[2], 0) +
                    (1.0 - sigma) * ls(left_index[0], left_index[1], left_index[2], 0);
                Set::Scalar QR = sigma * ls(right_index[0], right_index[1], right_index[2], 0) +
                    (1.0 - sigma) * ls(index[0], index[1], index[2], 0);

                // Store results
                QL_arr(index[0], index[1], index[2], 0) = QL;
                QR_arr(index[0], index[1], index[2], 0) = QR;
            });
        }
    }
};

// Local Lax-Friedrichs flux method
template <typename IntegratorDerivedClass>
class LocalLaxFriedrichsMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil) const override {
        //Util::Message(INFO, "Computing Local Lax-Friedrichs flux for direction: " + std::to_string(direction));

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats
        //Util::ScimitarX_Util::Debug debug;
        //debug.SetTargetDebugLocationIndices(101, 0, 0);  // Debug target location
        //debug.SetTargetDebugLocation(0.49, 0.12, 0.0,solver->geom[lev]);  // Debug target location

        const int nghosts = solver->number_of_ghost_cells;

        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {

            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts); //Includes two ghost layers

            auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                solver->ZFlux_mf.Patch(lev, mfi);
            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);

            // **Direction to velocity and energy mapping**
            int normal = direction;                 // Normal velocity component
            int trans1 = (direction + 1) % AMREX_SPACEDIM;  // First transverse velocity
#if AMREX_SPACEDIM == 3                                                              
            int trans2 = (direction + 2) % AMREX_SPACEDIM;  // Second transverse velocity (3D only)
#endif
            // **Direction to velocity mapping**
            int velocity_component[3] = {
                solver->variableIndex.UVEL,  // X-direction -> UVEL
                solver->variableIndex.VVEL,  // Y-direction -> VVEL
                solver->variableIndex.WVEL   // Z-direction -> WVEL
            };

            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {


                // Construct index arrays with correct number of components
                int index[3];
                // int left_index[3];
                // int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];

#if AMREX_SPACEDIM == 2
                index[0] = iface; index[1] = jface; index[2] = kface;
                //left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                //right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1; upper_bounds[2] = 0;
#else
                index[0] = iface; index[1] = jface; index[2] = kface;
                //left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                //right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1; upper_bounds[2] = face_bx_with_ghosts.bigEnd(2);
#endif

                // Clamp the main index
                ClampIndices(index, lower_bounds, upper_bounds);

                // Shift and clamp for left and right indices
                //int left_offset[3] = {0, 0, 0};
                //int right_offset[3] = {0, 0, 0};
                //left_offset[direction] = 0;  // No Left offset. For Flux left has the same indices as the cell center indices.
                //right_offset[direction] = 1;  // Shift right

                //ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                //ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                // **Left and right primitive variables (QL and QR)**
                Set::Scalar rho_L = QL_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
                Set::Scalar rho_R = QR_arr(index[0], index[1], index[2], solver->variableIndex.DENS);

                rho_L = std::max(rho_L, 1e-8);  // Avoid zero density
                rho_R = std::max(rho_R, 1e-8);  // Avoid zero density

                Set::Scalar e_L = std::max(QL_arr(index[0], index[1], index[2], solver->variableIndex.IE), 1e-8);
                Set::Scalar e_R = std::max(QR_arr(index[0], index[1], index[2], solver->variableIndex.IE), 1e-8);

                Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
                Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;

                Set::Scalar vel_L[3] = { 0.0, 0.0, 0.0 };
                Set::Scalar vel_R[3] = { 0.0, 0.0, 0.0 };

                vel_L[normal] = QL_arr(index[0], index[1], index[2], velocity_component[normal]);
                vel_R[normal] = QR_arr(index[0], index[1], index[2], velocity_component[normal]);

#if AMREX_SPACEDIM >= 2    
                vel_L[trans1] = QL_arr(index[0], index[1], index[2], velocity_component[trans1]);
                vel_R[trans1] = QR_arr(index[0], index[1], index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_L[trans2] = QL_arr(index[0], index[1], index[2], velocity_component[trans2]);
                vel_R[trans2] = QR_arr(index[0], index[1], index[2], velocity_component[trans2]);
#endif            


                // **Compute Kinetic Energy for Left and Right States**
                Set::Scalar KE_L = 0.5 * (vel_L[0] * vel_L[0] + vel_L[1] * vel_L[1] + vel_L[2] * vel_L[2]);  // Kinetic energy for left state
                Set::Scalar KE_R = 0.5 * (vel_R[0] * vel_R[0] + vel_R[1] * vel_R[1] + vel_R[2] * vel_R[2]);  // Kinetic energy for right state

                // **Total Internal Energy (E = e + KE)**
                Set::Scalar E_L = e_L + KE_L;  // Total internal energy for left state
                Set::Scalar E_R = e_R + KE_R;  // Total internal energy for right state

                // **Speed of sound**
                Set::Scalar c_L = std::sqrt(gamma * std::max(pressure_L, 1e-8) / rho_L);
                Set::Scalar c_R = std::sqrt(gamma * std::max(pressure_R, 1e-8) / rho_R);

                // **Maximum wave speed (lambda_max)**
                Set::Scalar lambda_max = std::max(std::abs(vel_L[normal]) + c_L, std::abs(vel_R[normal]) + c_R);


                //debug.CheckNegativePressureAndAbort(iface, jface, kface, pressure_L, "Pressure Field Left", "LLF Step", true);
                //debug.CheckNegativePressureAndAbort(iface, jface, kface, pressure_R, "Pressure Field Right", "LLF Step", true);

            // **Flux arrays**
                int num_components = solver->number_of_components;
                std::vector<Set::Scalar> F_L(num_components, 0.0);
                std::vector<Set::Scalar> F_R(num_components, 0.0);

                F_L[solver->variableIndex.DENS] = rho_L * vel_L[direction];  // Mass flux
                F_R[solver->variableIndex.DENS] = rho_R * vel_R[direction];

                F_L[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal] * vel_L[normal] + pressure_L;
                F_R[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal] * vel_R[normal] + pressure_R;

#if AMREX_SPACEDIM >= 2
                F_L[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[normal] * vel_L[trans1];
                F_R[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[normal] * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
                F_L[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[normal] * vel_L[trans2];
                F_R[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[normal] * vel_R[trans2];
#endif

                F_L[solver->variableIndex.IE] = (rho_L * E_L + pressure_L) * vel_L[normal];
                F_R[solver->variableIndex.IE] = (rho_R * E_R + pressure_R) * vel_R[normal];


                // **Conservative state vectors (left and right)**
                std::vector<Set::Scalar> QL_cons(num_components, 0.0);
                std::vector<Set::Scalar> QR_cons(num_components, 0.0);

                // **Convert primitive variables to conservative variables**
                QL_cons[solver->variableIndex.DENS] = rho_L;  // Density
                QR_cons[solver->variableIndex.DENS] = rho_R;

                QL_cons[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal];  // Momentum in normal direction
                QR_cons[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal];

#if AMREX_SPACEDIM >= 2
                QL_cons[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[trans1];  // Momentum in first transverse direction
                QR_cons[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[trans1];
#endif

#if AMREX_SPACEDIM == 3
                QL_cons[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[trans2];  // Momentum in second transverse direction
                QR_cons[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[trans2];
#endif

                QL_cons[solver->variableIndex.IE] = rho_L * E_L;  // Total energy
                QR_cons[solver->variableIndex.IE] = rho_R * E_R;


                // **Local Lax-Friedrichs flux calculation**
                for (int n = 0; n < num_components; ++n) {

                    flux_arr(index[0], index[1], index[2], n) = 0.5 * (F_L[n] + F_R[n]) - 0.5 * lambda_max * (QR_cons[n] - QL_cons[n]);

                }

                // Debug check for NaNs
                    // Debug
                    /*debug.DebugLocalLaxFriedrichsFlux(iface, jface, kface,num_components,
                                                    rho_L, vel_L[0], vel_L[1], vel_L[2], e_L, pressure_L,
                                                    rho_R, vel_R[0], vel_R[1], vel_R[2], e_R, pressure_R,
                                                    c_L, c_R, lambda_max,
                                                    F_L.data(), F_R.data(),
                                                    "ComputeFlux", true);

                debug.CheckNaNAndAbort(iface, jface, kface, lambda_max, "Lambda Max", "wave speed", true);*/
            });


        }

        //Util::Message(INFO, "Local Lax-Friedrichs flux computed for direction: " + std::to_string(direction));
    }
};

template <>
class LocalLaxFriedrichsMethod<Integrator::NarrowBandLevelset> : public FluxMethod<Integrator::NarrowBandLevelset> {
public:
    void ComputeFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver,
        const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil,
        amrex::MultiFab& temp_flux) const override {
        const int nghosts = solver->number_of_ghost_cells;  // Get number of ghost cells
        const int ncomp = solver->number_of_components;

        // Define OutsideNarrowband
        const int OutsideNarrowband = solver->GetOutsideNarrowBandPos();

        // Loop through all the levelset fields
        for (int ils = 0; ils < ncomp; ils++) {
            // Get the narrowband box array and distibution mapping for temporary multifabs
            amrex::BoxArray narrowband_ba = solver->level_sets[ils].narrowband_ba;
            amrex::DistributionMapping narrowband_dm = solver->level_sets[ils].narrowband_dm;

            // Create temporary multifabs
            amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, nghosts);
            amrex::MultiFab narrowband_Vel(narrowband_ba, narrowband_dm, AMREX_SPACEDIM, nghosts);

            // ParallelCopy & initialize
            temp_flux.setVal(0.0, nghosts); // Initialize because each ls has unique band
            auto const& Tube = solver->level_sets[ils].Tube_imf;
            auto const& velocity = solver->level_sets[ils].velocity_mf;
            narrowband.ParallelCopy(*Tube, 0, 0, 1, nghosts, nghosts);
            narrowband_Vel.ParallelCopy(*velocity[lev], 0, 0, AMREX_SPACEDIM, nghosts, nghosts);

            // Fill Boundary for ghost cells
            narrowband.FillBoundary();
            narrowband_Vel.FillBoundary();
            temp_flux.FillBoundary();

            // Loop through face-centered cells
            for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
                // Use a face-centered box to compute fluxes
                const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

                // Compute valid index bounds using compiler
                int lower_bounds[3] = { 0, 0, 0 };
                int upper_bounds[3] = { 0, 0, 0 };

#if AMREX_SPACEDIM >= 1
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0);
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
#endif

#if AMREX_SPACEDIM >= 2
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1);
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
#endif

#if AMREX_SPACEDIM == 3
                lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
                upper_bounds[2] = face_bx_with_ghosts.bigEnd(2) - 1;
#endif

                // Get multifab arrays
                auto const& QL_arr = QL_stencil.array(mfi);
                auto const& QR_arr = QR_stencil.array(mfi);
                auto const& nb_arr = narrowband.array(mfi);
                auto const& velocity_arr = narrowband_Vel.array(mfi);
                auto const& flux_arr = temp_flux.array(mfi);

                // Loop through face indices
                amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                    // Define cell centered indices/offsets
                    int cell_centered_index[3] = { iface, jface, kface };
                    int left_cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_index[3] = { iface, jface, kface };
                    int left_cell_centered_offset[3] = { 0, 0, 0 };
                    int right_cell_centered_offset[3] = { 0, 0, 0 };

                    // Define offsets (modular for future higher-order schemes)
                    left_cell_centered_offset[direction] = -1;  // Left offset for first-order upwind scheme
                    right_cell_centered_offset[direction] = 1;  // Right offset for first-order upwind scheme

                    // Clamp cell centered indices
                    ClampIndices(cell_centered_index, lower_bounds, upper_bounds);
                    ShiftAndClampIndices(left_cell_centered_index, left_cell_centered_offset, lower_bounds, upper_bounds, direction);
                    ShiftAndClampIndices(right_cell_centered_index, right_cell_centered_offset, lower_bounds, upper_bounds, direction);

                    // Create amrex::IntVect to store coords
                    const amrex::IntVect cell_center_IV(AMREX_D_DECL(cell_centered_index[0], 
                                                                     cell_centered_index[1], 
                                                                     cell_centered_index[2])); 
                    const amrex::IntVect left_cell_center_IV(AMREX_D_DECL(left_cell_centered_index[0], 
                                                                          left_cell_centered_index[1], 
                                                                          left_cell_centered_index[2]));
                    const amrex::IntVect right_cell_center_IV(AMREX_D_DECL(right_cell_centered_index[0], 
                                                                           right_cell_centered_index[1], 
                                                                           right_cell_centered_index[2]));

                    // Exit early if not narrowband cell
                    int nb_val = nb_arr(cell_center_IV);
                    if (std::abs(nb_val) == OutsideNarrowband) return;

                    // Clamp face centered indices
                    int face_centered_index[3] = { iface, jface, kface };
                    ClampIndices(face_centered_index, lower_bounds, upper_bounds);

                    const amrex::IntVect face_centered_IV(AMREX_D_DECL(face_centered_index[0], 
                                                                       face_centered_index[1], 
                                                                       face_centered_index[2]));

                    // Get velocities at left and right faces
                    Set::Scalar vel_L = velocity_arr(left_cell_center_IV, direction);
                    Set::Scalar vel_R = velocity_arr(right_cell_center_IV, direction);

                    // Compute max wave speed λ for LF flux
                    Set::Scalar lambda_max = std::max(std::abs(vel_L), std::abs(vel_R));

                    // Fetch reconstructed level-set values
                    Set::Scalar QL = QL_arr(face_centered_IV, ils);
                    Set::Scalar QR = QR_arr(face_centered_IV, ils);

                    // Compute Local Lax-Friedrichs numerical flux
                    Set::Scalar flux = 0.5 * (vel_L * QL + vel_R * QR) - 0.5 * lambda_max * (QR - QL);

                    // Store computed flux
                    flux_arr(face_centered_IV, ils) = flux;
                });
            }

            // Fill boundary for temp_flux
            temp_flux.FillBoundary();

            // Copy the temp_flux to full grid multifab
            auto& full_grid_flux = (direction == 0) ? solver->XFlux_mf[lev]:
                                   (direction == 1) ? solver->YFlux_mf[lev]:
                                    solver->ZFlux_mf[lev];
            full_grid_flux->ParallelCopy(temp_flux,0,0,ncomp,nghosts,nghosts);
        }
    }
};

template <>
class LocalLaxFriedrichsMethod<Integrator::NarrowBandLevelset_Backup> : public FluxMethod<Integrator::NarrowBandLevelset_Backup> {
public:
    void ComputeFlux(int direction, int lev, Integrator::NarrowBandLevelset_Backup* solver,
        const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil) const override {

        const int nghosts = solver->number_of_ghost_cells;  // Get number of ghost cells

        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            // Use a face-centered box to compute fluxes
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            // Select the correct flux MultiFab based on direction
            auto const& flux_arr = (direction == 0) ? solver->XFlux_mf.Patch(lev, mfi) :
                (direction == 1) ? solver->YFlux_mf.Patch(lev, mfi) :
                solver->ZFlux_mf.Patch(lev, mfi);

            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);
            auto const& velocity_arr = solver->velocity_mf.Patch(lev, mfi);

            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                int index[3] = { iface, jface, kface };
                int left_index[3] = { iface, jface, kface };
                int right_index[3] = { iface, jface, kface };

                // Compute valid index bounds
                int lower_bounds[3] = {
                    face_bx_with_ghosts.smallEnd(0),
                    (AMREX_SPACEDIM > 1) ? face_bx_with_ghosts.smallEnd(1) : 0,
                    (AMREX_SPACEDIM > 2) ? face_bx_with_ghosts.smallEnd(2) : 0
                };

                int upper_bounds[3] = {
                    face_bx_with_ghosts.bigEnd(0) - 1,
                    (AMREX_SPACEDIM > 1) ? face_bx_with_ghosts.bigEnd(1) - 1 : 0,
                    (AMREX_SPACEDIM > 2) ? face_bx_with_ghosts.bigEnd(2) - 1 : 0
                };

                // Clamp the main index to valid bounds
                ClampIndices(index, lower_bounds, upper_bounds);

                // Define offsets for left and right indices
                int left_offset[3] = { 0, 0, 0 };
                int right_offset[3] = { 0, 0, 0 };
                left_offset[direction] = -1;
                right_offset[direction] = 1;

                // Shift and clamp left/right indices
                ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                // Get velocities at left and right faces
                //Set::Scalar vel_L = velocity_arr(left_index[0], left_index[1], left_index[2])[direction];
                //Set::Scalar vel_R = velocity_arr(right_index[0], right_index[1], right_index[2])[direction];
                Set::Scalar vel_L = velocity_arr(left_index[0], left_index[1], left_index[2], direction);
                Set::Scalar vel_R = velocity_arr(right_index[0], right_index[1], right_index[2], direction);

                // Compute max wave speed λ for LF flux
                Set::Scalar lambda_max = std::max(std::abs(vel_L), std::abs(vel_R));

                // Fetch reconstructed level-set values
                Set::Scalar QL = QL_arr(index[0], index[1], index[2], 0);
                Set::Scalar QR = QR_arr(index[0], index[1], index[2], 0);

                // Compute Local Lax-Friedrichs numerical flux
                Set::Scalar flux = 0.5 * (vel_L * QL + vel_R * QR) - 0.5 * lambda_max * (QR - QL);

                // Store computed flux
                flux_arr(index[0], index[1], index[2], 0) = flux;
            });
        }
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H
