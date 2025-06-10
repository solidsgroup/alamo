#ifndef NUMERIC_NARROWBANDFLUXHANDLER_H
#define NUMERIC_NARROWBANDFLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/NarrowBandLevelset.H"
#include "FluxHandler.H"
#include "Util/Util.H"                              
#include <memory>
#include <cmath>
#include <algorithm>

// Include foward declare for ScimitarX.H
namespace Integrator {
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

// Base class for flux reconstruction
template <>
class FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxReconstruction() = default;

    /**
    * @brief Perform flux reconstruction for the specified direction.
    */
    virtual void ReconstructFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband,
        const amrex::MultiFab& narrowband_Vel, const amrex::MultiFab& narrowband_LS, 
        amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const = 0;
};

// NarrowBand specific flux class
template <>
class FluxMethod<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxMethod() = default;

    /**
    * @brief Compute the flux for the specified direction.
    */
    virtual void ComputeFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband, 
        const amrex::MultiFab& narrowband_Vel, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& temp_flux_mf, Set::Field<Set::Scalar>& full_flux_mf) const = 0;
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
    flux_mf_wrapper.Define(lev, flux_ba, dmap, flux_mf_wrapper[lev]->nComp(), flux_mf_wrapper[lev]->nGrow());
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
            // Create all temporary multifabs here so no need to perform multiple copies
            auto const& ls_data = solver->level_sets[ils];

            // Get proper distribution mapping for narrowband
            const amrex::BoxArray& narrowband_ba = ls_data.narrowband_ba;
            const amrex::DistributionMapping& narrowband_dm = ls_data.narrowband_dm;

            // Define cell centered multifabs from solver
            amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, nghosts);
            amrex::MultiFab narrowband_Vel(narrowband_ba, narrowband_dm, AMREX_SPACEDIM, nghosts);
            amrex::MultiFab narrowband_LS(narrowband_ba, narrowband_dm, 1, nghosts);
            
            // ParallelCopy
            auto const& Tube = ls_data.Tube_imf;
            auto const& velocity = ls_data.velocity_mf;
            narrowband.ParallelCopy(*Tube, 0, 0, 1, nghosts, nghosts);
            narrowband_Vel.ParallelCopy(*velocity[lev], 0, 0, AMREX_SPACEDIM, nghosts, nghosts);
            narrowband_LS.ParallelCopy(*solver->ls_mf[lev], ils, 0, 1, nghosts, nghosts);

            // Loop through all dimensions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Convert narrowband_ba to face-centered direction dir
                amrex::BoxArray narrowband_face_ba(narrowband_ba);  
                amrex::IntVect narrowband_face_type = amrex::IntVect::TheZeroVector(); 
                narrowband_face_type[dir] = 1;
                narrowband_face_ba.convert(narrowband_face_type);

                // Create temporary QL and QR MultiFabs at face-centered positions
                amrex::MultiFab QL_stencil(narrowband_face_ba, narrowband_dm, ncomp, nghosts);
                amrex::MultiFab QR_stencil(narrowband_face_ba, narrowband_dm, ncomp, nghosts);
                
                // Create temporary flux to store flux values for narrowband cells
                amrex::MultiFab temp_flux_mf(narrowband_face_ba, narrowband_dm, ncomp, nghosts);
                temp_flux_mf.setVal(0.0, nghosts);

                // Get the corresponding full flux MultiFab and its index type
                auto& full_flux_mf =
                    (dir == 0) ? solver->XFlux_mf:
                    (dir == 1) ? solver->YFlux_mf:
                                solver->ZFlux_mf;

                // Get the index types for temp_flux_mf and full_flux_mf
                const amrex::IndexType full_flux_ixType = full_flux_mf[lev]->boxArray().ixType();
                const amrex::IndexType temp_flux_ixType = temp_flux_mf.boxArray().ixType();

                // Check if index types are same and if not, convert
                // THERE IS BUG IN HANDLING OF FACE_CENTERED MFS - NEED TO FIND SOURCE
                // THIS IS ONLY TEMPORARY FIX
                if (full_flux_ixType != temp_flux_ixType) ConvertToFaceCentered(dir, lev, full_flux_mf);

                // Zero out the full domain flux so that no remnant fluxes are present
                full_flux_mf[lev]->setVal(0.0, nghosts);

                // Perform reconstruction and flux computation
                reconstruction->ReconstructFlux(dir, lev, solver, narrowband,
                    narrowband_Vel, narrowband_LS, QL_stencil, QR_stencil);
                flux_method->ComputeFlux(dir, lev, solver, narrowband, narrowband_Vel, QL_stencil, QR_stencil, temp_flux_mf, full_flux_mf);
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

// First-order specific to NarrowBandLevelset integrator class
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset> : public FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    void ReconstructFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband,
        const amrex::MultiFab& narrowband_Vel, const amrex::MultiFab& narrowband_LS, 
        amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
        // Define multifab constants 
        const int nghosts = solver->number_of_ghost_cells;
        const int ncomp = solver->number_of_components;

        // Define OutsideNarrowband
        const int OutsideNarrowband = solver->GetOutsideNarrowBandPos();

        const amrex::IntVect debug(AMREX_D_DECL(12, 12, 0));

        // Loop through all the levelset fields
        for (int ils = 0; ils < ncomp; ils++) {
            // Fill Boundary
            QL_stencil.FillBoundary(ils, 1); // Fill only the current component
            QR_stencil.FillBoundary(ils, 1);

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
                auto const& nb_arr = narrowband.const_array(mfi);
                auto const& vel_arr = narrowband_Vel.const_array(mfi);
                auto const& ls = narrowband_LS.const_array(mfi);
                auto const& QL_arr = QL_stencil.array(mfi);
                auto const& QR_arr = QR_stencil.array(mfi);
               
                // Loop through all face cell indices
                amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                    // Construct cell-centered indices
                    int cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_offset[3] = { 0, 0, 0 };

                    // Define offsets (modular for future higher-order schemes)
                    right_cell_centered_offset[direction] = 1;  // Right offset for first-order upwind scheme

                    // Clamp cell centered indices
                    ClampIndices(cell_centered_index, lower_bounds, upper_bounds);
                    ShiftAndClampIndices(right_cell_centered_index, right_cell_centered_offset, lower_bounds, upper_bounds, direction);

                    // Create amrex::IntVect to store coords
                    const amrex::IntVect cell_center_IV(AMREX_D_DECL(cell_centered_index[0], 
                                                                    cell_centered_index[1], 
                                                                    cell_centered_index[2])); 
                    const amrex::IntVect right_cell_center_IV(AMREX_D_DECL(right_cell_centered_index[0], 
                                                                        right_cell_centered_index[1], 
                                                                        right_cell_centered_index[2]));

                    // Debug print
                    //if (cell_center_IV == debug) amrex::Print() << "debugging point " << debug << " for direction " << direction << std::endl;

                    // Exit early if not narrowband cell
                    int nb_val = nb_arr(cell_center_IV);
                    //if (cell_center_IV == debug) amrex::Print() << "Narrowband: " << nb_val << std::endl;
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
                    Set::Scalar QL = sigma * ls(cell_center_IV) + (1.0 - sigma) * ls(right_cell_center_IV);
                    Set::Scalar QR = sigma * ls(right_cell_center_IV) + (1.0 - sigma) * ls(cell_center_IV);

                    // Store results
                    QL_arr(face_centered_IV, ils) = QL;
                    QR_arr(face_centered_IV, ils) = QR;

                    // Debug
                    /*if (cell_center_IV == debug){
                        amrex::Print() << "Velocity: " << vel << std::endl;
                        amrex::Print() << "LS: " << ls(cell_center_IV) << std::endl;
                        amrex::Print() << "LS Right: " << ls(right_cell_center_IV) << std::endl;
                        amrex::Print() << "QL: " << QL << std::endl;
                        amrex::Print() << "QR: " << QR << std::endl;
                    }*/
                });
            }

            // Fill Boundary after updating array
            QL_stencil.FillBoundary(ils, 1); // Fill only the current component
            QR_stencil.FillBoundary(ils, 1);
        }
    }
};

template <>
class LocalLaxFriedrichsMethod<Integrator::NarrowBandLevelset> : public FluxMethod<Integrator::NarrowBandLevelset> {
public:
    void ComputeFlux(int direction, int lev, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband,
        const amrex::MultiFab& narrowband_Vel, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& temp_flux_mf, Set::Field<Set::Scalar>& full_flux_mf) const override {
        const int nghosts = solver->number_of_ghost_cells;  // Get number of ghost cells
        const int ncomp = solver->number_of_components;

        // Define OutsideNarrowband
        const int OutsideNarrowband = solver->GetOutsideNarrowBandPos();

        const amrex::IntVect debug(AMREX_D_DECL(12, 12, 0));

        // Loop through all the levelset fields
        for (int ils = 0; ils < ncomp; ils++) {
            // Fill Boundary for ghost cells
            temp_flux_mf.FillBoundary(ils, 1);

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
                auto const& nb_arr = narrowband.const_array(mfi);
                auto const& velocity_arr = narrowband_Vel.const_array(mfi);
                auto const& QL_arr = QL_stencil.const_array(mfi);
                auto const& QR_arr = QR_stencil.const_array(mfi);
                auto const& flux_arr = temp_flux_mf.array(mfi);

                // Loop through face indices
                amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                    // Define cell centered indices/offsets
                    int cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_index[3] = { iface, jface, kface };
                    int right_cell_centered_offset[3] = { 0, 0, 0 };

                    // Define offsets (modular for future higher-order schemes)
                    right_cell_centered_offset[direction] = 1;  // Left offset for first-order upwind scheme

                    // Clamp cell centered indices
                    ClampIndices(cell_centered_index, lower_bounds, upper_bounds);
                    ShiftAndClampIndices(right_cell_centered_index, right_cell_centered_offset, lower_bounds, upper_bounds, direction);

                    // Create amrex::IntVect to store coords
                    const amrex::IntVect cell_center_IV(AMREX_D_DECL(cell_centered_index[0], 
                                                                    cell_centered_index[1], 
                                                                    cell_centered_index[2])); 
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
                    Set::Scalar VL = velocity_arr(cell_center_IV, direction);
                    Set::Scalar VR = velocity_arr(right_cell_center_IV, direction);
                    Set::Scalar VFace = 0.5 * (VL + VR); // Linear appriximation for cell face

                    // Get lambda_max
                    Set::Scalar lambda_max = std::max(std::abs(VL), std::abs(VR));

                    // Fetch reconstructed level-set values
                    Set::Scalar QL = QL_arr(face_centered_IV, ils);
                    Set::Scalar QR = QR_arr(face_centered_IV, ils);

                    // Compute Local Lax-Friedrichs numerical flux
                    Set::Scalar flux = VL * QL; //0.5 * VFace * (QL + QR) - 0.5 * lambda_max * (QR - QL);

                    // Store computed flux
                    flux_arr(face_centered_IV, ils) = flux;
                    /*if (cell_center_IV == debug){
                        amrex::Print() << "Debugging point " << debug << " for direction " << direction << " in compute flux" << std::endl;
                        amrex::Print() << "VL: " << VL << std::endl;
                        amrex::Print() << "QL: " << QL << std::endl;
                        amrex::Print() << "temp_flux_mf: " << flux << std::endl;
                    }*/
                });
            }

            // Fill boundary before copy
            temp_flux_mf.FillBoundary(ils, 1);

            // Copy the temp_flux to full grid multifab
            full_flux_mf[lev]->ParallelCopy(temp_flux_mf, 0, 0, ncomp, nghosts, nghosts);
        }
    }
};

}  // namespace Numeric

#endif  // NUMERIC_NARROWBANDFLUXHANDLER_H
