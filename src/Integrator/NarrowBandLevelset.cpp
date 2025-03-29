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
//#include "Util/NarrowBandLevelset_Util.H"

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

    // Initialize the zero levelset imf
    //Zerols_imf.reset(new amrex::iMultiFab(grids[lev], dmap[lev], 1, number_of_ghost_cells));
    Zerols_imf.reset(new amrex::iMultiFab(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, number_of_ghost_cells));
    ClearZerols();

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

void NarrowBandLevelset::ClearZerols(){;
    Zerols_imf->setVal(-1);
}

void NarrowBandLevelset::UpdateZerols(int lev, int ls_id) {
    // Define grid spacing min_DX to find zero-crossing cells
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar threshold = (narrow_band_width + 1) * min_DX;

    // Define neighbor offsets based on dimension
    #if AMREX_SPACEDIM == 1
    constexpr int num_neighbors = 2;
    const int dx[] = {-1, 1};
    #endif

    #if AMREX_SPACEDIM == 2
    constexpr int num_neighbors = 4;
    const int dx[] = {-1, 1,  0, 0};
    const int dy[] = { 0, 0, -1, 1};
    #endif

    #if AMREX_SPACEDIM == 3
    constexpr int num_neighbors = 6;
    const int dx[] = {-1, 1,  0, 0,  0, 0};
    const int dy[] = { 0, 0, -1, 1,  0, 0};
    const int dz[] = { 0, 0,  0, 0, -1, 1};
    #endif

    // Update ghost cells before iterating over tileboxes
    ls_mf[lev]->FillBoundary();

    for (amrex::MFIter mfi(*Zerols_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& domain_box = geom[lev].Domain();

        // Get arrays for MultiFabs
        auto const& ls_arr = ls_mf[lev]->array(mfi);
        auto const& zerols_arr = Zerols_imf->array(mfi);  // FIXED ACCESS

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Clamp phi values
            Set::Scalar phi = ls_arr(i, j, k, ls_id);
            phi = std::clamp(phi, -threshold, threshold);
            ls_arr(i, j, k, ls_id) = phi;

            // Compute |phi|
            Set::Scalar abs_phi = std::abs(phi);

            // Early exit: If already a zero-cell, store ls_id and return
            if (abs_phi <= min_DX) {
                zerols_arr(i, j, k) = ls_id;
                return;
            }

            // Check neighbors for sign change (zero-crossing)
            for (int n = 0; n < num_neighbors; ++n) {
                int ni = i + dx[n];
                int nj = j;
                int nk = k;

                #if AMREX_SPACEDIM >= 2
                nj = j + dy[n];
                #endif
                #if AMREX_SPACEDIM == 3
                nk = k + dz[n];
                #endif

                // Ensure neighbor is inside domain
                if (!domain_box.contains(ni, nj, nk)) {
                    continue;
                }

                Set::Scalar phi_neighbor = ls_arr(ni, nj, nk, ls_id);
                if (phi * phi_neighbor < 0) {
                    zerols_arr(i, j, k) = ls_id;
                    return;  // No need to check more neighbors
                }
            }
        });
    }
}

void NarrowBandLevelset::ComputeNarrowBandBoxList(int lev, int ls_id) {
    // Access structure
    auto & ls_data = level_sets[ls_id];

    // Define narrowband variables
    amrex::BoxList narrow_band_boxes;

    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    amrex::Gpu::DeviceVector<int> narrow_band_flags(ls_mf[lev]->local_size(), 0);
    auto* d_narrow_band_flags = narrow_band_flags.data();

    amrex::Vector<amrex::Box> temp_boxes;
    int box_index = 0;

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.tilebox();
        auto const& phi_arr = ls_mf.Patch(lev, mfi);

        int* flag_ptr = d_narrow_band_flags + box_index;

        // Scan through the box to find narrow band region
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Set::Scalar phi = phi_arr(i, j, k, ls_id);

            #ifdef AMREX_USE_GPU
                if (std::abs(phi) <= tube_width) {
                    amrex::Gpu::Atomic::Or(flag_ptr, 1);  // GPU-safe atomic operation
                }
            #else
                if (std::abs(phi) <= tube_width) {
                    *flag_ptr = 1;  // Direct assignment for CPU
                }
            #endif
        });

        temp_boxes.push_back(valid_box);
        box_index++;
    }

    // Copy back results to check which boxes contain the narrow band
    amrex::Gpu::Device::synchronize();
    std::vector<int> h_narrow_band_flags(narrow_band_flags.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, narrow_band_flags.begin(), narrow_band_flags.end(), h_narrow_band_flags.begin());

    for (int i = 0; i < temp_boxes.size(); i++) {
        if (h_narrow_band_flags[i] != 0) {
            narrow_band_boxes.push_back(temp_boxes[i]);
        }
    }

    // Update the level set's narrow band boxes
    ls_data.narrowband_boxes = narrow_band_boxes;
}

/*void NarrowBandLevelset::ComputeNarrowBandBoxList(int lev, int ls_id){
    // Access structure
    auto & ls_data = level_sets[ls_id];

    // Define narrowband variables
    amrex::BoxList narrow_band_boxes;

    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    // Loop through boxes to find narrow band regions
    for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();
        auto const& phi_arr = ls_mf.Patch(lev, mfi);

        // Initialize min/max bounds to track narrow band region
        bool has_narrow_band = false;

        // Scan through the box to find narrow band region
        amrex::ParallelFor(valid_box, [=, &has_narrow_band] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Get phi
            Set::Scalar phi = phi_arr(i, j, k, ls_id);
            if (std::abs(phi) <= tube_width) {
                has_narrow_band = true;
            }
        });

        // If narrow band cells found, create a box for them
        if (has_narrow_band) {
            amrex::IntVect min_ind = valid_box.smallEnd();
            amrex::IntVect max_ind = valid_box.bigEnd();
            amrex::Box nb_box(min_ind, max_ind);
            narrow_band_boxes.push_back(nb_box);
        }
    }
    
    // Update the level set's narrow band boxes
    ls_data.narrowband_boxes = narrow_band_boxes; 
}*/

void NarrowBandLevelset::ComputeNarrowBandMapping(int lev, int ls_id){
    auto& ls_data = level_sets[ls_id];

    // Check if there are any narrowband boxes
    amrex::BoxList narrowband_boxes = ls_data.narrowband_boxes;
    if (narrowband_boxes.size()> 0) {
        amrex::BoxArray ba(narrowband_boxes);
        ls_data.narrowband_dm = amrex::DistributionMapping(ba);
        ls_data.has_narrowband = true;  // Flag indicating valid narrowband
    } else {
        ls_data.has_narrowband = false; // No narrowband region found
    }
}

void NarrowBandLevelset::UpdateNarrowBandFlags(int lev, int ls_id){
    // Get the size of the box array
    int ba_size = ls_mf[lev]->boxArray().size();

    // Get the narrowband boxes
    amrex::BoxList narrowband_boxes = level_sets[ls_id].narrowband_boxes;

    // Initialize the narrowband_flags
    amrex::Vector<int> narrowband_flags(ba_size, 0);

    for (int i = 0; i < ba_size; ++i) {
        const amrex::Box& vb = ls_mf[lev]->boxArray()[i];

        for (const amrex::Box& nb : narrowband_boxes) {
            if (vb.intersects(nb)) {	
                narrowband_flags[i] = 1;
                break;
            }
        }
    }

    // Save to level_sets
    level_sets[ls_id].narrowband_flags = narrowband_flags;
}

void NarrowBandLevelset::UpdateNarrowband(int lev, int ls_id){
    // Tag 0 levelset gridpoints - zero will be done outside loop
    UpdateZerols(lev, ls_id);

    // Define narrowband box list 
    ComputeNarrowBandBoxList(lev, ls_id);

    // Get the distribution mapping for processors
    ComputeNarrowBandMapping(lev, ls_id);

    // Get the narrowband flags
    UpdateNarrowBandFlags(lev, ls_id);
}

void NarrowBandLevelset::ComputeNormal(int lev, int ls_id){
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();

    // Update boundaries for ghost cells **before** iterating over tileboxes
    ls_mf[lev]->FillBoundary();

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Only perform computations in narrowband
        const int idx = mfi.index();  // Index into box array
        amrex::Vector<int> narrowband_flags = level_sets[ls_id].narrowband_flags; // flags

        if (narrowband_flags[idx] == 0) continue;  // Skip non-narrowband boxes

        // Add a ghost layer for central differencing scheme
        const amrex::Box& grown_bx = mfi.growntilebox(1);
        auto const& phi_arr = ls_mf.Patch(lev, mfi);
        auto const& normal = level_sets[ls_id].normal_mf.Patch(lev, mfi); 

        amrex::ParallelFor(grown_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Utilize Numeric::Gradient for simplified computations
            Set::Vector gradient = Numeric::Gradient(phi_arr, i, j, k, ls_id, DX);
            Set::Scalar grad_mag = gradient.norm();
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
        // Only perform computations in narrowband
        const int idx = mfi.index();  // Index into box array
        amrex::Vector<int> narrowband_flags = level_sets[ls_id].narrowband_flags; // flags
        
        if (narrowband_flags[idx] == 0) continue;  // Skip non-narrowband boxes
        
        // Add a ghost layer for central differencing scheme
        const amrex::Box& grown_bx = mfi.growntilebox(1);
        auto const& normal = level_sets[ls_id].normal_mf.Patch(lev, mfi); 
        auto const& curvature = level_sets[ls_id].curvature_mf.Patch(lev, mfi);

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

    // Clear the zero levelset
    ClearZerols();

    // Update the narrowband
    for (int ls_id = 0; ls_id < number_of_components; ls_id++){
        UpdateNarrowband(lev, ls_id);
        ComputeGeometryQuantities(lev, ls_id);
    }

    // Reinitialize
    if(current_timestep % 5 == 0) {
        Reinitialize(lev);

        // Update the narrowband
        for (int ls_id = 0; ls_id < number_of_components; ls_id++){
            UpdateNarrowband(lev, ls_id);
            ComputeGeometryQuantities(lev, ls_id);
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
                // Only perform computations in narrowband
                const int idx = mfi.index();  // Index into box array
                amrex::Vector<int> narrowband_flags = level_sets[ils].narrowband_flags; // flags
        
                if (narrowband_flags[idx] == 0) continue;  // Skip non-narrowband boxes

                const amrex::Box& valid_bx = mfi.tilebox();
                auto const& ls_arr = ls_mf.Patch(lev, mfi);
                auto const& zerols_arr = Zerols_imf->array(mfi);  // FIXED ACCESS

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (zerols_arr(i,j,k) == ils) return;
  
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
