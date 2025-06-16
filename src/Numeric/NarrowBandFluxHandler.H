#ifndef NUMERIC_NARROWBANDFLUXHANDLER_H
#define NUMERIC_NARROWBANDFLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/NarrowBandLevelset.H"
#include "FluxHandler.H"
#include "Util/Util.H"
#include "Numeric/Stencil.H"                                
#include <memory>
#include <cmath>
#include <algorithm>

// Include foward declare for NarrowBandLevelset.H
namespace Integrator {
    class NarrowBandLevelset;
}

namespace Numeric {
// Base class for flux reconstruction
template <>
class FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxReconstruction() = default;

    /**
    * @brief Perform flux reconstruction for the specified direction.
    */
    virtual void ReconstructFlux(int direction, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband,
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
    virtual void ComputeFlux(int direction, int lev, int ils, Integrator::NarrowBandLevelset* solver, const amrex::iMultiFab& narrowband, 
        const amrex::MultiFab& narrowband_Vel, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& temp_flux_mf) const = 0;
};

/*inline void ConvertToFaceCentered(int dir, int lev, Set::Field<Set::Scalar>& flux_mf_wrapper) {
    amrex::Print() << "In Convert Function for dir " << dir << std::endl;

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
}*/

/*void ApplyDirectionalNeumannBC(
    amrex::MultiFab& narrowband_LS,
    const amrex::iMultiFab& narrowband,
    int nghosts,
    int dir,  // 0 = x, 1 = y, 2 = z
    int OuterGhost,
    int OuterTube)
{
    AMREX_ALWAYS_ASSERT(narrowband_LS.nGrow() >= nghosts);
    AMREX_ALWAYS_ASSERT(narrowband.nGrow() >= nghosts);

    for (amrex::MFIter mfi(narrowband_LS, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& nb_arr = narrowband.const_array(mfi);
        auto const& phi_arr = narrowband_LS.array(mfi);
        const amrex::Box& bx = mfi.growntilebox(nghosts);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            int tag = nb_arr(iv);

            // Only apply to ghost cells
            if (std::abs(tag) != OuterGhost) return;

            // Check both directions along the specified axis
            for (int s = -1; s <= 1; s += 2) {
                for (int offset = 1; offset <= nghosts; ++offset) {
                    amrex::IntVect shift(AMREX_D_DECL(0, 0, 0));
                    shift[dir] = -s * offset;
                    amrex::IntVect src_iv = iv + shift;

                    if (!bx.contains(src_iv)) continue;

                    int src_tag = nb_arr(src_iv);
                    if (std::abs(src_tag) <= OuterTube) {
                        phi_arr(iv) = phi_arr(src_iv);
                        return;
                    }
                }
            }
        });
    }

    // Final fill of ghost regions from updated data
    narrowband_LS.FillBoundary();
}*/

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

        // Get solver variables
        const int ncomp = solver->number_of_components;
        const int nghosts = solver->number_of_ghost_cells;

        // Loop through all levelset fields
        for (int ils = 0; ils < ncomp; ++ils) {
            // Create all temporary multifabs here so no need to perform multiple copies
            auto const& ls_data = solver->level_sets[lev][ils];

            // Get proper distribution mapping for narrowband
            const amrex::BoxArray& narrowband_ba = ls_data.narrowband_ba;
            const amrex::DistributionMapping& narrowband_dm = ls_data.narrowband_dm;

            // Define narrowband variables
            const auto& narrowband = *ls_data.Tube;
            const auto& ls = *ls_data.Tube_ls;
            const auto& velocity = *ls_data.Tube_velocity;

            // Loop through all dimensions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Convert narrowband_ba to face-centered direction dir
                amrex::BoxArray narrowband_face_ba(narrowband_ba);  
                amrex::IntVect narrowband_face_type = amrex::IntVect::TheCellVector(); 
                narrowband_face_type[dir] = 1;
                narrowband_face_ba.convert(narrowband_face_type);

                // Create temporary QL and QR MultiFabs at face-centered positions
                amrex::MultiFab QL_stencil(narrowband_face_ba, narrowband_dm, 1, nghosts);
                amrex::MultiFab QR_stencil(narrowband_face_ba, narrowband_dm, 1, nghosts);
                
                // Create temporary flux to store flux values for narrowband cells
                amrex::MultiFab temp_flux_mf(narrowband_face_ba, narrowband_dm, 1, nghosts);
                temp_flux_mf.setVal(0.0);

                // Perform reconstruction and flux computation
                reconstruction->ReconstructFlux(dir, solver, narrowband,
                    velocity, ls, QL_stencil, QR_stencil);
                flux_method->ComputeFlux(dir, lev, ils, solver, narrowband, velocity, QL_stencil, QR_stencil, temp_flux_mf);

            }
        }
    }

private:
    std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> flux_method;             // Flux method
};

// First-order specific to NarrowBandLevelset integrator class
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset> : public FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    void ReconstructFlux(int direction, Integrator::NarrowBandLevelset* solver, 
        const amrex::iMultiFab& narrowband, 
        const amrex::MultiFab& narrowband_Vel, 
        const amrex::MultiFab& narrowband_LS, 
        amrex::MultiFab& QL_stencil, 
        amrex::MultiFab& QR_stencil) const override {
        // Define multifab constants 
        const int nghosts = 1;

        // Define OutsideNarrowband
        const int OutsideTube = solver->GetOutsideTube();

        const amrex::IntVect debug(AMREX_D_DECL(19, 49, 0));
        const amrex::IntVect debug2(AMREX_D_DECL(20, 49, 0));

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
                int nb_L = std::abs(nb_arr(cell_center_IV));
                int nb_R = std::abs(nb_arr(right_cell_center_IV));
                if (nb_L == OutsideTube && nb_R == OutsideTube) return;

                //int nb_val = nb_arr(cell_center_IV);
                //if (cell_center_IV == debug) amrex::Print() << "Narrowband: " << nb_val << std::endl;
                //if (std::abs(nb_val) == OutsideNarrowband) return;

                // Clamp face centered indices
                int face_centered_index[3] = { iface, jface, kface };
                ClampIndices(face_centered_index, lower_bounds, upper_bounds);
                const amrex::IntVect face_centered_IV(AMREX_D_DECL(face_centered_index[0], 
                                                                face_centered_index[1], 
                                                                face_centered_index[2]));

                // Compute velocity at face and its sign (branch-free)
                Set::Scalar V = vel_arr(cell_center_IV, direction);
                int sigma = V > 0 ? 1 : 0;

                // Get left and right LS values
                Set::Scalar LS = ls(cell_center_IV);
                Set::Scalar LS_R = ls(right_cell_center_IV);

                // Perform reconstruction
                Set::Scalar QL = sigma * LS + (1.0 - sigma) * LS_R;
                Set::Scalar QR = sigma * LS_R + (1.0 - sigma) * LS;

                // Store results
                QL_arr(face_centered_IV, 0) = QL;
                QR_arr(face_centered_IV, 0) = QR;

                // Debug
                /*if ((cell_center_IV == debug  || cell_center_IV == debug2) && direction == 0 ){
                    amrex::Print() << "Debugging point " << cell_center_IV<< " for direction " << direction << " in reconstuction" << std::endl;
                    amrex::Print() << "nb_left: " << nb_L << std::endl;
                    amrex::Print() << "nb_right: " << nb_R << std::endl;
                    amrex::Print() << "Velocity: " << V << std::endl;
                    amrex::Print() << "LS: " << LS << std::endl;
                    amrex::Print() << "LS Right: " << LS_R << std::endl;
                    amrex::Print() << "QL: " << QL << std::endl;
                    amrex::Print() << "QR: " << QR << std::endl;
                }*/
            });
        }

        // Fill Boundary after updating arrays
        QL_stencil.FillBoundary();
        QR_stencil.FillBoundary();
    }
};

template <>
class LocalLaxFriedrichsMethod<Integrator::NarrowBandLevelset> : public FluxMethod<Integrator::NarrowBandLevelset> {
public:
    void ComputeFlux(int direction, int lev, int ils, Integrator::NarrowBandLevelset* solver, 
        const amrex::iMultiFab& narrowband,
        const amrex::MultiFab& narrowband_Vel, 
        const amrex::MultiFab& QL_stencil, 
        const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& temp_flux_mf) const override {
        const int nghosts = 1; 

        // Define OutsideNarrowband
        const int OutsideTube = solver->GetOutsideTube();

        const amrex::IntVect debug(AMREX_D_DECL(19, 49, 0));
        const amrex::IntVect debug2(AMREX_D_DECL(20, 49, 0));

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

            // Define multifab arrays
            auto const& nb_arr = narrowband.const_array(mfi);
            auto const& vel_arr = narrowband_Vel.const_array(mfi);
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
                int nb_L = nb_arr(cell_center_IV);
                int nb_R = nb_arr(right_cell_center_IV);
                if (std::abs(nb_L) == OutsideTube && std::abs(nb_R) == OutsideTube) return;

                // Clamp face centered indices
                int face_centered_index[3] = { iface, jface, kface };
                ClampIndices(face_centered_index, lower_bounds, upper_bounds);

                const amrex::IntVect face_centered_IV(AMREX_D_DECL(face_centered_index[0], 
                                                                    face_centered_index[1], 
                                                                    face_centered_index[2]));

                // Get velocity of cell center
                Set::Scalar V = vel_arr(cell_center_IV, direction);

                // Fetch reconstructed level-set values
                Set::Scalar QL = QL_arr(face_centered_IV, 0);
                Set::Scalar QR = QR_arr(face_centered_IV, 0);

                // Compute first order numerical flux
                Set::Scalar flux = V * QL;

                // Store computed flux
                flux_arr(face_centered_IV, 0) = flux;
                /*if ((cell_center_IV == debug || cell_center_IV == debug2) && direction  == 0){
                    amrex::Print() << "Debugging point " << cell_center_IV << " for direction " << direction << " in compute flux" << std::endl;
                    amrex::Print() << "V: " << V << std::endl;
                    amrex::Print() << "QL: " << QL << std::endl;
                    amrex::Print() << "QR: " << QR << std::endl;
                    amrex::Print() << "temp_flux_mf: " << flux << std::endl;
                }*/
            });
        }

        // Fill ghost cells
        temp_flux_mf.FillBoundary();

        // Copy the temp_flux to full grid multifab
        auto& full_flux_mf =
            (direction == 0) ? solver->XFlux_mf:
            (direction == 1) ? solver->YFlux_mf:
                                solver->ZFlux_mf;
        full_flux_mf[lev]->ParallelCopy(temp_flux_mf, 0, ils, 1, nghosts, nghosts);
    }
};

}  // namespace Numeric

#endif  // NUMERIC_NARROWBANDFLUXHANDLER_H
