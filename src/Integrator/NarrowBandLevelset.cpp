// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "Integrator/NarrowBandLevelset.H"

// BC
#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Expression.H"

// IC
#include "IC/Constant.H"
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"
#include "NarrowBandLevelset.H"

#include "Numeric/FluxHandler.H"
#include "Numeric/TimeStepper.H"

namespace Integrator
{

// Define constructor functions
// Empty constructor
NarrowBandLevelset::NarrowBandLevelset(int a_nghost) 
    : Integrator(), 
    number_of_ghost_cells(a_nghost) 
{}
    
// Constructor that triggers Parse
NarrowBandLevelset::NarrowBandLevelset(IO::ParmParse& pp) : NarrowBandLevelset() // Call default constructor
{
    fluxHandler = std::make_shared<Numeric::FluxHandler<NarrowBandLevelset>>();
    timeStepper = std::make_shared<Numeric::TimeStepper<NarrowBandLevelset>>();

    pp.queryclass(*this); // Call the static Parse function
}

// Define Parse function
void NarrowBandLevelset::Parse(NarrowBandLevelset& value, IO::ParmParse& pp){
    {// Define initial and boundary conditions
    
    // Query the IC assuming either LS::Sphere or LS::Zalesak
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak>("ic.ls",value.ic_ls,value.geom);
    
    // Assume Neumann BC for levelset field
    value.bc_ls = new BC::Constant(value.number_of_components, pp, "bc.ls");
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    }

    // Define levelset data structure
    value.level_sets.resize(value.number_of_components);
    
    {// Initialize constant velocity. Will be removed once integrated with ScimitarX integrator
    // Define constant velocity vector
    pp_queryarr("ic.velocity.value", value.constant_velocity);

    // Define velocity, normal, and curvature multifabs
    for (int ils=0; ils < value.number_of_components; ils++){
        value.RegisterGeneralFab(value.level_sets[ils].velocity_mf, AMREX_SPACEDIM, value.number_of_ghost_cells, "velocity", false); 
        value.RegisterGeneralFab(value.level_sets[ils].normal_mf, AMREX_SPACEDIM, value.number_of_ghost_cells, "normal", false);
        value.RegisterNewFab(value.level_sets[ils].curvature_mf, value.bc_ls, 1, value.number_of_ghost_cells, "curvature", false);   
    }
    }


    {// Register face-centered flux fields
    value.RegisterFaceFab<0>(value.XFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "xflux", false);
    #if AMREX_SPACEDIM >= 2
    value.RegisterFaceFab<1>(value.YFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "yflux", false);
    #endif
    #if AMREX_SPACEDIM == 3
    value.RegisterFaceFab<2>(value.ZFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "zflux", false);
    #endif
    }

    {// Access the maps from the FeatureMaps singleton
    auto& fluxReconstructionMap = getFeatureMaps().getFluxReconstructionMap();
    auto& fluxSchemeMap = getFeatureMaps().getFluxSchemeMap();
    auto& timeSteppingSchemeMap = getFeatureMaps().getTimeSteppingSchemeMap();
    
    // Flux Reconstruction parsing
    std::string fluxReconstructionStr;
    if (pp.query("FluxReconstruction", fluxReconstructionStr)) {
        auto it = fluxReconstructionMap.find(fluxReconstructionStr);
        if (it != fluxReconstructionMap.end()) {
            value.reconstruction_method = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid FluxReconstruction value: " + fluxReconstructionStr);
        }
    }
    
    // Flux Scheme parsing
    std::string fluxSchemeStr;
    if (pp.query("FluxScheme", fluxSchemeStr)) {
        auto it = fluxSchemeMap.find(fluxSchemeStr);
        if (it != fluxSchemeMap.end()) {
            value.flux_scheme = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid FluxScheme value: " + fluxSchemeStr);
        }
    }
    
    // Time-Stepping Scheme parsing
    std::string timeSteppingSchemeStr;
    if (pp.query("TimeSteppingScheme", timeSteppingSchemeStr)) {
        auto it = timeSteppingSchemeMap.find(timeSteppingSchemeStr);
        if (it != timeSteppingSchemeMap.end()) {
            value.temporal_scheme = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid TimeSteppingScheme value: " + timeSteppingSchemeStr);
        }
    }
    }

    // Read CFL number and initial time step
    pp.query_required("cflNumber", value.cflNumber);
}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
         
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    ic_ls->Initialize(lev, ls_mf);

    Zerols_imf.reset(new amrex::iMultiFab(grids[lev], dmap[lev], 1, number_of_ghost_cells));

    // Loop through all levelsets
    for (int ils=0; ils < number_of_components; ils++){
        // Get structure id number
        level_sets[ils].id = ils;

        // Define velocity_mf containing velocity vector
        level_sets[ils].velocity_mf[lev] -> setVal(constant_velocity);

        // Define initial narrowband boxarray and distribution mapping
        UpdateNarrowband(lev, ils);

        // Define Geometric quantities - must be after narrowband to update narrowband boxes!
        ComputeGeometryQuantities(lev, ils);
    } 
}

void NarrowBandLevelset::ClearZerols(){
    Zerols_imf.reset();
}

void NarrowBandLevelset::ComputeNarrowBandBox(int lev, int ls_id){
    // Access structure
    auto & ls_data = level_sets[ls_id];

    // Define narrowband variables
    amrex::BoxList narrow_band_boxes;

    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tube_width = narrow_band_width * min_DX;
    const Set::Scalar threshold = (narrow_band_width + 1) * min_DX;

    // Loop through boxes to find narrow band regions
    for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& tilebox = mfi.tilebox();
        auto const& phi_arr = ls_mf.Patch(lev, mfi);

        // Initialize min/max bounds to track narrow band region
        amrex::IntVect min_idx = tilebox.bigEnd();
        amrex::IntVect max_idx = tilebox.smallEnd();
        bool has_narrow_band = false;

        // Scan through the box to find narrow band region
        amrex::ParallelFor(tilebox, [=, &has_narrow_band, &min_idx, &max_idx] 
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Clamp the values of phi
            Set::Scalar phi = phi_arr(i, j, k, ls_id);
            phi = std::clamp(phi, -threshold, threshold);
            phi_arr(i, j, k, ls_id) = phi;

            if (std::abs(phi) <= tube_width) {
                has_narrow_band = true;
                min_idx[0] = std::min(min_idx[0], i);
                max_idx[0] = std::max(max_idx[0], i);
#if (AMREX_SPACEDIM >= 2)
                min_idx[1] = std::min(min_idx[1], j);
                max_idx[1] = std::max(max_idx[1], j);
#endif
#if (AMREX_SPACEDIM == 3)
                min_idx[2] = std::min(min_idx[2], k);
                max_idx[2] = std::max(max_idx[2], k);
#endif
            }
        });

        // If narrow band cells found, create a box for them
        if (has_narrow_band) {
            amrex::Box nb_box(min_idx, max_idx);
            nb_box = nb_box.grow(1); // Add small buffer
            narrow_band_boxes.push_back(nb_box);
        }
    }
    
    // Update the level set's narrow band boxes
    ls_data.narrowband_boxes = amrex::BoxArray(narrow_band_boxes); 
}

void NarrowBandLevelset::ComputeNarrowBandMapping(int ls_id){
    auto& ls_data = level_sets[ls_id];

    // Check if there are any narrowband boxes
    if (ls_data.narrowband_boxes.size() > 0) {
        ls_data.narrowband_dm = amrex::DistributionMapping(ls_data.narrowband_boxes);
        ls_data.has_narrowband = true;  // Flag indicating valid narrowband
    } else {
        ls_data.has_narrowband = false; // No narrowband region found
    }
}

bool NarrowBandLevelset::BoxIntersectsNarrowBand(const amrex::Box& box, int ls_id, int lev) const
{
    const auto& nb_boxes = level_sets[ls_id].narrowband_boxes;
    
    // Fast rejection test - does this box intersect any narrow band box?
    for (int i = 0; i < nb_boxes.size(); ++i) {
        if (box.intersects(nb_boxes[i])) {
            return true;
        }
    }
    
    return false;
}

void NarrowBandLevelset::UpdateNarrowband(int lev, int ls_id){
    // Compute zero levelset object tagging

    // Define narrowband box list 
    ComputeNarrowBandBox(lev, ls_id);

    // Get the distribution mapping for processors
    ComputeNarrowBandMapping(ls_id);
}

void NarrowBandLevelset::ComputeNormal(int lev, int ls_id){
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();

    // Update boundaries for ghost cells **before** iterating over tileboxes
    ls_mf[lev]->FillBoundary();

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Add a ghost layer for central differencing scheme
        const amrex::Box& grown_bx = mfi.growntilebox(1);
        auto const& phi_arr = ls_mf.Patch(lev, mfi);
        auto const& normal = level_sets[ls_id].normal_mf.Patch(lev, mfi); 

        // Only perform computations in narrowband
        if (!BoxIntersectsNarrowBand(grown_bx, ls_id, lev)) {
            continue;
        }

        amrex::ParallelFor(grown_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Utilize Numeric::Gradient for simplified computations
            Set::Vector gradient = Numeric::Gradient(phi_arr, i, j, k, ls_id, DX);
            Set::Scalar grad_mag = 1.0; // FIX THIS
            normal(i, j, k) = gradient / grad_mag;
        });
    }
}

void NarrowBandLevelset::ComputeCurvature(int lev, int ls_id){
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();

    // Update boundaries for ghost cells **before** iterating over tileboxes
    level_sets[ls_id].normal_mf[lev]->FillBoundary();
    
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Add a ghost layer for central differencing scheme
        const amrex::Box& grown_bx = mfi.growntilebox(1);
        auto const& normal = level_sets[ls_id].normal_mf.Patch(lev, mfi); 
        auto const& curvature = level_sets[ls_id].curvature_mf.Patch(lev, mfi);

        // Only perform computations in narrowband
        if (!BoxIntersectsNarrowBand(grown_bx, ls_id, lev)) {
            continue;
        }

        amrex::ParallelFor(grown_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Utilize Numeric::Gradient for simplified computations
            Set::Matrix gradient_matrix = Numeric::Gradient(normal, i, j, k, DX);

            // Get diagonals of matrix
            for (int dim=0; dim < AMREX_SPACEDIM; dim++){
                curvature(i, j, k) += gradient_matrix(dim, dim);
            }
        });
    }  
}

void NarrowBandLevelset::ComputeGeometryQuantities(int lev, int ls_id){
    ComputeNormal(lev, ls_id);
    ComputeCurvature(lev, ls_id);  
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) { 
    ComputeAndSetNewTimeStep();
}

void NarrowBandLevelset::ComputeAndSetNewTimeStep() {
    // Compute the minimum time step using the CFL condition
    Set::Scalar finest_dt = GetTimeStep();  // GetTimeStep() already applies the CFL number

    // Start with the finest-level time step
    Set::Scalar coarsest_dt = finest_dt;

    // Adjust time step for coarser levels based on refinement ratios
    for (int lev = finest_level; lev > 0; --lev) {
        int refinement_factor = refRatio(lev - 1)[0];  // Assume isotropic refinement
        coarsest_dt *= refinement_factor;  // Scale conservatively for refinement
    }

    // Set the time step for all levels
    Integrator::SetTimestep(coarsest_dt);
}

Set::Scalar NarrowBandLevelset::GetTimeStep() {
    Set::Scalar minDt = std::numeric_limits<Set::Scalar>::max();  // Start with a large value 
    
    // Loop through all levelsets
    for (int ils=0; ils < number_of_components; ils++){
        for (int lev = 0; lev <= maxLevel(); ++lev) {  // Iterate over AMR levels
            const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`
            const Set::Scalar min_DX = *std::min_element(dx, dx + AMREX_SPACEDIM);
    
            for (amrex::MFIter mfi(*level_sets[ils].velocity_mf[lev], false); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
                auto const& velocity_arr = level_sets[ils].velocity_mf.Patch(lev, mfi);  // Velocity field
    
                Set::Scalar minDt_local = std::numeric_limits<Set::Scalar>::max();  // Thread-local minDt
    
                amrex::ParallelFor(bx, [=, &minDt_local](int i, int j, int k) noexcept {
                    // Extract velocity components
                    Set::Scalar u = velocity_arr(i, j, k)[0];
    #if (AMREX_SPACEDIM >= 2)
                    Set::Scalar v = velocity_arr(i, j, k)[1];
    #endif
    #if (AMREX_SPACEDIM == 3)
                    Set::Scalar w = velocity_arr(i, j, k)[2];
    #endif
    
                    // Compute the maximum velocity magnitude
                    Set::Scalar maxSpeed = std::abs(u);
    #if (AMREX_SPACEDIM >= 2)
                    maxSpeed = std::max(maxSpeed, std::abs(v));
    #endif
    #if (AMREX_SPACEDIM == 3)
                    maxSpeed = std::max(maxSpeed, std::abs(w));
    #endif
    
                    // Compute local CFL time step restriction
                    if (maxSpeed > 1e-8) {  // Avoid division by zero
                        Set::Scalar dtLocal = dx[0] / maxSpeed;
    #if (AMREX_SPACEDIM >= 2)
                        dtLocal = std::min(dtLocal, dx[1] / maxSpeed);
    #endif
    #if (AMREX_SPACEDIM == 3)
                        dtLocal = std::min(dtLocal, dx[2] / maxSpeed);
    #endif
                        minDt_local = std::min(minDt_local, dtLocal);
                    }
                    else {
                         minDt_local = min_DX;
                    }
                });
    
                // Update the global minDt
                minDt = std::min(minDt, minDt_local);
            }
        }
    }

    // Reduce across MPI processes to find the global minimum timestep
    amrex::ParallelDescriptor::ReduceRealMin(minDt);
    
    return cflNumber * minDt;  // Return CFL-adjusted time step
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) 
{
    // Update current timestep (for reinitialize)
    current_timestep ++;

    // Advect
    Advect(lev, time, dt);

    // Update the narrowband
    for (int ls_id = 0; ls_id < number_of_components; ls_id++){
        UpdateNarrowband(lev, ls_id);
    }

    // Reinitialize
    if(current_timestep % 5 == 0) {
        Reinitialize(lev);

        // Update the narrowband
        for (int ls_id = 0; ls_id < number_of_components; ls_id++){
            UpdateNarrowband(lev, ls_id);
        }
    }
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar time, Set::Scalar dt){
    // Swap the levelsets
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);

    // Update the interface velocity
    UpdateInterfaceVelocity(lev);

    // Compute the flux
    switch (temporal_scheme) {
        case TimeSteppingScheme::ForwardEuler: {
            
            fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<NarrowBandLevelset>>());
            fluxHandler->SetFluxMethod(std::make_shared<Numeric::LocalLaxFriedrichsMethod<NarrowBandLevelset>>());

            timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<NarrowBandLevelset>>());

            int numStages = timeStepper->GetNumberOfStages();
            // One-stage loop for Forward Euler
            for (int stage = 0; stage < numStages; ++stage) {

                // 2. Perform flux reconstruction and compute fluxes in all directions
                fluxHandler->ConstructFluxes(lev, this);

                // 3. Compute sub-step using the chosen time-stepping scheme
                timeStepper->ComputeSubStep(lev, dt, stage, this);
            }
            break;
        }
    } 
}

void NarrowBandLevelset::UpdateInterfaceVelocity(int lev){
    for (int ils=0; ils < number_of_components; ils++){
        for (amrex::MFIter mfi(*level_sets[ils].velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const& vel_arr = level_sets[ils].velocity_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector velocity = vel_arr(i,j,k);
                vel_arr(i,j,k) = velocity;
            });
        }
    }
}

void NarrowBandLevelset::Reinitialize(int lev) {
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50; 

    // Grid spacing
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);

    // Time constraint for PDE update
    const Set::Scalar tau = cflNumber * min_DX;

    // Loop through all components
    for (int ils=0; ils < number_of_components; ils++){
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Use integer flag for convergence: 0 (not converged), 1 (converged)
            amrex::Gpu::DeviceScalar<int> d_converged(1);  // Start as converged (1)
            int* d_converged_ptr = d_converged.dataPtr();
            
            // Update boundaries for ghost cells **before** iterating over tileboxes
            ls_mf[lev]->FillBoundary();
    
            for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& valid_bx = mfi.tilebox();
                auto const& ls_arr = ls_mf.Patch(lev, mfi);

                // Test if box intersects with narrowband. If not, skip
                if (!BoxIntersectsNarrowBand(valid_bx, ils, lev)) {
                    continue;
                }
    
                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // **Compute sign of phi (avoid 0 cases)**
                    Set::Scalar phi_val = ls_arr(i, j, k, ils);
                    Set::Scalar sign_phi = (phi_val > 0) ? 1.0 : -1.0;
                    
                    // Compute first-order upwind gradients
                    Set::Scalar a_p = std::max(sign_phi * (phi_val - ls_arr(i - 1, j, k, ils)) / DX[0], 0.0);
                    Set::Scalar a_m = std::min(sign_phi * (ls_arr(i + 1, j, k, ils) - phi_val) / DX[0], 0.0);
    
                    Set::Scalar b_p = 0.0, b_m = 0.0, c_p = 0.0, c_m = 0.0;
                    #if AMREX_SPACEDIM >= 2
                    b_p = std::max(sign_phi * (phi_val - ls_arr(i, j - 1, k, ils)) / DX[1], 0.0);
                    b_m = std::min(sign_phi * (ls_arr(i, j + 1, k, ils) - phi_val) / DX[1], 0.0);
                    #endif
                    #if AMREX_SPACEDIM == 3
                    c_p = std::max(sign_phi * (phi_val - ls_arr(i, j, k - 1, ils)) / DX[2], 0.0);
                    c_m = std::min(sign_phi * (ls_arr(i, j, k + 1, ils) - phi_val) / DX[2], 0.0);
                    #endif
    
                    // **Compute gradient magnitude**
                    Set::Scalar nablaG = std::sqrt(
                        std::max(a_p * a_p, a_m * a_m) +
                        std::max(b_p * b_p, b_m * b_m) +
                        std::max(c_p * c_p, c_m * c_m)
                    );
    
                    // **Update level set function**
                    Set::Scalar phi_new = phi_val - tau * sign_phi * (nablaG - 1.0);
    
                    // **Compute local convergence check**
                    int local_converged = (std::abs(phi_new - phi_val) <= reinit_tolerance) ? 1 : 0;
    
                    // **Atomic update: Track if any thread detects non-convergence**
                    amrex::Gpu::Atomic::Max(d_converged_ptr, local_converged);
    
                    ls_arr(i, j, k, ils) = phi_new;
                });
            }

            // Check convergence across processors
            int converged = d_converged.dataValue();
            if (converged) break;  // Stop early if converged
        }
    }
}

void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow){
}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
}

} // namespace Integrator
