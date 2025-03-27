// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "Integrator/NarrowBandLevelset_Backup.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/TimeStepper.H"

// BC
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "BC/Expression.H"

// IC
#include "IC/Constant.H"
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"

// Numeric
#include "Util/Util.H"
#include "Util/NarrowBandLevelset_Util.H"

namespace Integrator
{

// Define constructor functions
// Empty constructor
NarrowBandLevelset_Backup::NarrowBandLevelset_Backup(int a_nghost) 
    : Integrator(), 
    number_of_ghost_cells(a_nghost) 
{}
    
// Constructor that triggers Parse
NarrowBandLevelset_Backup::NarrowBandLevelset_Backup(IO::ParmParse& pp) : NarrowBandLevelset_Backup() // Call default constructor
{
    fluxHandler = std::make_shared<Numeric::FluxHandler<NarrowBandLevelset_Backup>>();
    timeStepper = std::make_shared<Numeric::TimeStepper<NarrowBandLevelset_Backup>>();

    pp.queryclass(*this); // Call the static Parse function
}

// Define Parse function
void NarrowBandLevelset_Backup::Parse(NarrowBandLevelset_Backup& value, IO::ParmParse& pp){
    {// Define initial and boundary conditions
    
    // Query the IC assuming either LS::Sphere or LS::Zalesak
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak>("ic.ls",value.ic_ls,value.geom);
    
    // Assume Neumann BC for levelset field
    value.bc_ls = new BC::Constant(value.number_of_components, pp, "bc.ls");
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    
    // Define narrowband multifab to track tube
    value.RegisterNewFab(value.iNarrowBandMask_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "NBand", true);
    }
    
    {// Following section will be deleted when velocity is taken from ScimitarX integrator class  
    // Initialize velocity IC and BC
    /*value.velocity_mf.Define(
        value.geom.size(),                // Number of AMR levels
        value.grids,                      // BoxArray for each level
        value.dmap,                       // DistributionMapping for each level
        1,                                 // nComp = 1 (each cell stores a full Set::Vector)
        value.number_of_ghost_cells        // Number of ghost cells
    );*/
    
    //std::cout << "Velocity MultiFab components: " << value.velocity_mf.NComp() << std::endl;
    //std::cout << "Velocity MultiFab ghost cells: " << value.number_of_ghost_cells << std::endl;

    
    value.ic_velocity = new IC::Constant(value.geom, pp, "ic.velocity");
    value.bc_velocity = new BC::Constant(AMREX_SPACEDIM, pp, "bc.velocity");
    
    // Define velocity multifab
    value.RegisterNewFab(value.velocity_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "velocity", false); 
    }
    
    // Register face-centered flux fields
    value.RegisterFaceFab<0>(value.XFlux_mf, &value.bc_nothing, 1, value.number_of_ghost_cells, "xflux", false);
#if AMREX_SPACEDIM >= 2
    value.RegisterFaceFab<1>(value.YFlux_mf, &value.bc_nothing, 1, value.number_of_ghost_cells, "yflux", false);
#endif
#if AMREX_SPACEDIM == 3
    value.RegisterFaceFab<2>(value.ZFlux_mf, &value.bc_nothing, 1, value.number_of_ghost_cells, "zflux", false);
#endif

    {
    // Access the maps from the FeatureMaps singleton
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
void NarrowBandLevelset_Backup::Initialize(int lev){
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    
    // COPY INSTEAD OF INITIALIZE!
    ic_ls->Initialize(lev, ls_mf);
    
    // Debug
    //Util::NarrowBandLevelset_Util::Debug::printNeighbors2D(124, 126, 251, 252, *ls_mf[lev], geom[lev], lev, "ls", "After Initialize");
    
    // Initialize narrowband
    UpdateNarrowBand(lev);
    //Util::NarrowBandLevelset_Util::Debug::printNeighbors2D(124, 126, 251, 252, *iNarrowBandMask_mf[lev], geom[lev], lev, "NB", "After NarrowBand");
    //Util::NarrowBandLevelset_Util::Debug::printNeighbors2D(124, 126, 251, 252, *ls_mf[lev], geom[lev], lev, "ls", "After NarrowBand");
    
    // Initialize velocity field -- will be deleted
    /*#if AMREX_SPACEDIM == 3
        Set::Vector dirichlet_value(1.0, 0.0, 0.0);  // 3D case
    #elif AMREX_SPACEDIM == 2
        Set::Vector dirichlet_value(1.0, 0.0);       // 2D case
    #endif

    std::cout << "Initializing velocity field on level " << lev << std::endl;
    std::cout.flush();
    
    if (!velocity_mf[lev]) {
        std::cout << "Error: velocity_mf[" << lev << "] is uninitialized!" << std::endl;
        return;
    }

    for (amrex::MFIter mfi(*velocity_mf[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        
            std::cout << "Rank " << amrex::ParallelDescriptor::MyProc() 
              << " - Processing box: " << bx << std::endl;

    if (bx.isEmpty()) {
        std::cout << "Tilebox is empty! ParallelFor will not execute." << std::endl;
        continue; // Skip empty boxes
    }
        amrex::Array4<Set::Vector> velocity = velocity_mf.Patch(lev, mfi);

        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            velocity(i, j, k) = dirichlet_value; // Example: Set velocity to (1,0,0)
            printf("vel(%d, %d, %d): <%f, %f\n>",i,j,k,velocity(i,j,k)[0],velocity(i,j,k)[1]); 
        });
    }*/
    ic_velocity->Initialize(lev, velocity_mf); 
}

void NarrowBandLevelset_Backup::UpdateNarrowBand(int lev) {
    // Define grid spacing min_DX to find 0 cells 
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);

    // Define tube widths similar to Fortran
    const Set::Scalar inner_tube = 4.0 * min_DX;
    const Set::Scalar outer_tube = 6.0 * min_DX;
    const Set::Scalar tube_buffer = 7.0 * min_DX;  

    // Define neighbor offsets based on dimension
    #if AMREX_SPACEDIM == 1
    constexpr int num_neighbors = 2;
    const int dx[] = {-1, 1};
    #endif

    #if AMREX_SPACEDIM == 2
    constexpr int num_neighbors = 8;
    const int dx[] = {-1, 1,  0, 0, -1, -1,  1,  1};
    const int dy[] = { 0, 0, -1, 1, -1,  1, -1,  1};
    #endif

    #if AMREX_SPACEDIM == 3
    constexpr int num_neighbors = 26;
    const int dx[] = {-1, 1,  0, 0,  0, 0, -1, -1,  1,  1, -1,  1, -1,  1, -1, -1,  1,  1,  0,  0,  0,  0, -1,  1, -1,  1};
    const int dy[] = { 0, 0, -1, 1,  0, 0, -1,  1, -1,  1,  0,  0,  0,  0, -1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1,  1};
    const int dz[] = { 0, 0,  0, 0, -1, 1,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0, -1,  1, -1,  1, -1, -1,  1,  1};
    #endif
    
    // Update ghost cells before iterating over tileboxes
    ls_mf[lev]->FillBoundary();
    
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& ls_arr = ls_mf.Patch(lev, mfi);
        auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Get current level set value
            Set::Scalar phi = ls_arr(i, j, k);
            phi = std::clamp(phi, -tube_buffer, tube_buffer);
            ls_arr(i, j, k) = phi;

            // Compute |phi|
            Set::Scalar abs_phi = std::abs(phi);

            // **If outside the tube, assign ±3 and return early**
            if (abs_phi > outer_tube) {
                nbmask_arr(i, j, k) = (phi > 0) ? 3 : -3;
                return;
            }

            // **Detect interface cells (zero-crossing)**
            bool is_interface = abs_phi <= min_DX; 

            // Get box limits for bounds checking
            const amrex::Box domain = bx;
            
            // **Loop through neighbors to check for opposite sign**
            for (int n = 0; n < num_neighbors; ++n) {
                #if AMREX_SPACEDIM == 1
                int ni = i + dx[n];
                int nj = j;  
                int nk = k;  
                #endif

                #if AMREX_SPACEDIM == 2
                int ni = i + dx[n];
                int nj = j + dy[n];
                int nk = k;  
                #endif

                #if AMREX_SPACEDIM == 3
                int ni = i + dx[n];
                int nj = j + dy[n];
                int nk = k + dz[n];
                #endif

                // **Check if neighbor is inside the computational domain**
                if (!domain.contains(amrex::IntVect(AMREX_D_DECL(ni, nj, nk)))) {
                    continue;  // Skip out-of-bounds neighbors
                }

                // Check if a neighbor has an opposite sign
                if ((phi) * (ls_arr(ni, nj, nk)) < 0) {
                    is_interface = true;
                    break;  // Stop checking once we find an interface
                }
            }

            // **Assign `nbmask = 0` only for true zero-crossings**
            if (is_interface) {
                nbmask_arr(i, j, k) = 0;
                return;
            }

            // **Otherwise, assign narrowband values**
            int sign_phi = (phi > 0) - (phi < 0);
            int inside_inner_tube = (abs_phi <= inner_tube);
            int inside_outer_tube = (abs_phi > inner_tube) & (abs_phi <= outer_tube);

            nbmask_arr(i, j, k) = inside_inner_tube * (sign_phi * 1) +
                                  inside_outer_tube * (sign_phi * 2);
        });
    }
}

void NarrowBandLevelset_Backup::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset_Backup::TimeStepComplete(Set::Scalar time, int lev) {
    ComputeAndSetNewTimeStep(); // Compute dt based on global `minDt   
}

void NarrowBandLevelset_Backup::ComputeAndSetNewTimeStep() {
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

Set::Scalar NarrowBandLevelset_Backup::GetTimeStep() {
    Set::Scalar minDt = std::numeric_limits<Set::Scalar>::max();  // Start with a large value 
      
    for (int lev = 0; lev <= maxLevel(); ++lev) {  // Iterate over AMR levels
        const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`
        const Set::Scalar min_DX = *std::min_element(dx, dx + AMREX_SPACEDIM);

        for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
            auto const& velocity_arr = velocity_mf.Patch(lev, mfi);  // Velocity field

            Set::Scalar minDt_local = std::numeric_limits<Set::Scalar>::max();  // Thread-local minDt

            amrex::ParallelFor(bx, [=, &minDt_local](int i, int j, int k) noexcept {
                // Extract velocity components
                //Set::Scalar u = velocity_arr(i, j, k)[0];
                //printf("u: %f\n",u);
                Set::Scalar u = velocity_arr(i, j, k, 0);
#if (AMREX_SPACEDIM >= 2)
                //Set::Scalar v = velocity_arr(i, j, k)[1];
                //printf("v: %f\n",v);
                Set::Scalar v = velocity_arr(i, j, k, 1);
#endif
#if (AMREX_SPACEDIM == 3)
                //Set::Scalar w = velocity_arr(i, j, k)[2];
                Set::Scalar w = velocity_arr(i, j, k, 2);
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

    // Reduce across MPI processes to find the global minimum timestep
    amrex::ParallelDescriptor::ReduceRealMin(minDt);
    
    return cflNumber * minDt;  // Return CFL-adjusted time step
}

void NarrowBandLevelset_Backup::Advance(int lev, Set::Scalar time, Set::Scalar dt) 
{
    // Update current timestep (for reinitialize)
    current_timestep ++;
    
    // Advect
    //Advect(lev, time, dt);
    
    // Update narrowband
    //UpdateNarrowBand(lev);
    //Util::NarrowBandLevelset_Util::Debug::printNeighbors2D(124, 126, 251, 252, *ls_mf[lev], geom[lev], lev, "ls", "After Advect");
    
    // Reinitialize the level set function
    if(current_timestep % 5 == 0) {
        Reinitialize(lev, time);

        // Update narrowband after reinitialization
        //UpdateNarrowBand(lev);
    }     
}

void NarrowBandLevelset_Backup::Advect(int lev, Set::Scalar time, Set::Scalar dt)
{
    // Swap the old ls fab and the new ls fab so we use
    // the new one.
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);

    // Update the velocity
    UpdateInterfaceVelocity(lev);

    // Compute the flux
    switch (temporal_scheme) {
        case TimeSteppingScheme::ForwardEuler: {
            
            fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<NarrowBandLevelset_Backup>>());
            fluxHandler->SetFluxMethod(std::make_shared<Numeric::LocalLaxFriedrichsMethod<NarrowBandLevelset_Backup>>());

            timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<NarrowBandLevelset_Backup>>());

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

void NarrowBandLevelset_Backup::UpdateInterfaceVelocity(int lev)
{
    /*// Loop through and set all velocity components
    for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        //Set::Patch<Set::Scalar>  velocity = velocity_mf.Patch(lev,mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Do computation here - right now placeholder to keep velocity constant
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                velocity(i,j,k,d) = velocity(i,j,k,d);
            }
        });
    }*/
}

/*void NarrowBandLevelset::Reinitialize(int lev, Set::Scalar time) {
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50; 

    // Grid spacing
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);

    // Time step constraint for reinitialization PDE
    const Set::Scalar tau = cflNumber * min_DX;

    // Create alias for in-place update
    amrex::MultiFab ls_reinit(*ls_mf[lev], amrex::make_alias, 0, 1);

    for (int iter = 0; iter < max_iterations; ++iter) {
        printf("iter: %d\n", iter);

        // Device scalar for GPU atomic convergence check
        amrex::Gpu::DeviceScalar<int> d_converged(1);  // Start as converged (1)
        int* d_converged_ptr = d_converged.dataPtr();

        for (amrex::MFIter mfi(ls_reinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& valid_bx = mfi.tilebox();
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

            amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                int nbmask_val = nbmask_arr(i, j, k, 0);

                // **Skip cells outside the narrowband (`nbmask == ±3`) and interface (`nbmask == 0`)**
                if (std::abs(nbmask_val) > 2 || nbmask_val == 0) return;

                // **Compute sign of phi (avoid zero cases)**
                Set::Scalar phi_val = ls_arr(i, j, k, 0);
                Set::Scalar sign_phi = (phi_val > 0) - (phi_val < 0); // Ensures ±1, avoids zero issues

                // **Compute first-order upwind gradients**
                Set::Scalar dx_p = (ls_arr(i + 1, j, k, 0) - phi_val) / DX[0];
                Set::Scalar dx_m = (phi_val - ls_arr(i - 1, j, k, 0)) / DX[0];

                Set::Scalar a_p = std::max(sign_phi * dx_m, 0.0);
                Set::Scalar a_m = std::min(sign_phi * dx_p, 0.0);

                Set::Scalar b_p = 0.0, b_m = 0.0, c_p = 0.0, c_m = 0.0;
                #if AMREX_SPACEDIM >= 2
                Set::Scalar dy_p = (ls_arr(i, j + 1, k, 0) - phi_val) / DX[1];
                Set::Scalar dy_m = (phi_val - ls_arr(i, j - 1, k, 0)) / DX[1];

                b_p = std::max(sign_phi * dy_m, 0.0);
                b_m = std::min(sign_phi * dy_p, 0.0);
                #endif

                #if AMREX_SPACEDIM == 3
                Set::Scalar dz_p = (ls_arr(i, j, k + 1, 0) - phi_val) / DX[2];
                Set::Scalar dz_m = (phi_val - ls_arr(i, j, k - 1, 0)) / DX[2];

                c_p = std::max(sign_phi * dz_m, 0.0);
                c_m = std::min(sign_phi * dz_p, 0.0);
                #endif

                // **Compute gradient magnitude**
                Set::Scalar nablaG = std::sqrt(
                    std::max(a_p * a_p, a_m * a_m) +
                    std::max(b_p * b_p, b_m * b_m) +
                    std::max(c_p * c_p, c_m * c_m)
                );

                // **Update level set function using PDE**
                Set::Scalar phi_new = phi_val - tau * sign_phi * (nablaG - 1.0);

                // **Compute local convergence check**
                int local_converged = (std::abs(phi_new - phi_val) <= reinit_tolerance) ? 1 : 0;

                // **Atomic update: Track if any thread detects non-convergence**
                amrex::Gpu::Atomic::Min(d_converged_ptr, local_converged);

                ls_arr(i, j, k, 0) = phi_new;
            });
        }

        // Check convergence across processors
        int converged = d_converged.dataValue();
        if (converged) break;  // Stop early if converged

        // Apply boundary conditions to ghost cells
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);
    }
}*/

void NarrowBandLevelset_Backup::Reinitialize(int lev, Set::Scalar time) {
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50; 

    // Grid spacing
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);

    // Time constraint for PDE update
    const Set::Scalar tau = cflNumber * min_DX;

    // Create alias for in-place update
    //amrex::MultiFab ls_reinit(*ls_mf[lev], amrex::make_alias, 0, 1);

    for (int iter = 0; iter < max_iterations; ++iter) {
        printf("iter: %d\n", iter);

        // Use integer flag for convergence: 0 (not converged), 1 (converged)
        amrex::Gpu::DeviceScalar<int> d_converged(1);  // Start as converged (1)
        int* d_converged_ptr = d_converged.dataPtr();
        
        // Update boundaries for ghost cells **before** iterating over tileboxes
        ls_mf[lev]->FillBoundary();

        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& valid_bx = mfi.tilebox();
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

            amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                int nbmask_val = nbmask_arr(i, j, k, 0);

                // **Skip cells outside the narrowband (`nbmask == ±3`)**
                //if (std::abs(nbmask_val) > 2 || nbmask_val == 0) return;

                // **Compute sign of phi (avoid 0 cases)**
                Set::Scalar phi_val = ls_arr(i, j, k, 0);
                Set::Scalar sign_phi = (phi_val > 0) ? 1.0 : -1.0;
                
                // Compute first-order upwind gradients
                Set::Scalar a_p = std::max(sign_phi * (phi_val - ls_arr(i - 1, j, k, 0)) / DX[0], 0.0);
                Set::Scalar a_m = std::min(sign_phi * (ls_arr(i + 1, j, k, 0) - phi_val) / DX[0], 0.0);

                Set::Scalar b_p = 0.0, b_m = 0.0, c_p = 0.0, c_m = 0.0;
                #if AMREX_SPACEDIM >= 2
                b_p = std::max(sign_phi * (phi_val - ls_arr(i, j - 1, k, 0)) / DX[1], 0.0);
                b_m = std::min(sign_phi * (ls_arr(i, j + 1, k, 0) - phi_val) / DX[1], 0.0);
                #endif
                #if AMREX_SPACEDIM == 3
                c_p = std::max(sign_phi * (phi_val - ls_arr(i, j, k - 1, 0)) / DX[2], 0.0);
                c_m = std::min(sign_phi * (ls_arr(i, j, k + 1, 0) - phi_val) / DX[2], 0.0);
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

                ls_arr(i, j, k, 0) = phi_new;
            });
        }

        // Check convergence across processors
        int converged = d_converged.dataValue();
        if (converged) break;  // Stop early if converged

        // Apply BCs to ghost cells
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);
    }
}

/*void NarrowBandLevelset::Reinitialize(int lev, Set::Scalar time) {
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50; 

    // Grid spacing
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);

    // Time constraint for PDE
    const Set::Scalar tau = cflNumber * min_DX;

    // Create alias for in-place update
    amrex::MultiFab ls_reinit(*ls_mf[lev], amrex::make_alias, 0, 1);

    for (int iter = 0; iter < max_iterations; ++iter) {
        printf("iter: %d\n", iter);
        bool converged = true;
        amrex::Gpu::DeviceScalar<bool> d_converged(true);
        bool* d_converged_ptr = d_converged.dataPtr();

        for (amrex::MFIter mfi(ls_reinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& valid_bx = mfi.tilebox();
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

            amrex::ParallelFor(valid_bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k) {
                // **Only update cells inside the narrowband (`-1, 0, 1`)**
                int nbmask_val = nbmask_arr(i, j, k, 0);
                if (std::abs(nbmask_val) > 2) return;  // Skip outside region

                // **Ensure interface (`nbmask == 0`) maintains correct sign**
                Set::Scalar phi_val = ls_arr(i, j, k, 0);
                if (nbmask_val == 0) {
                    ls_arr(i, j, k, 0) = 0.0;
                    return;  // No further processing for zero-crossing cells
                }

                Set::Scalar sign_phi = std::copysign(1.0, phi_val);

                // Compute first-order upwind gradients
                Set::Scalar a_p = std::max(sign_phi * (phi_val - ls_arr(i - 1, j, k, 0)) / DX[0], 0.0);
                Set::Scalar a_m = std::min(sign_phi * (ls_arr(i + 1, j, k, 0) - phi_val) / DX[0], 0.0);

                Set::Scalar b_p = 0.0, b_m = 0.0, c_p = 0.0, c_m = 0.0;
                #if AMREX_SPACEDIM >= 2
                b_p = std::max(sign_phi * (phi_val - ls_arr(i, j - 1, k, 0)) / DX[1], 0.0);
                b_m = std::min(sign_phi * (ls_arr(i, j + 1, k, 0) - phi_val) / DX[1], 0.0);
                #endif
                #if AMREX_SPACEDIM == 3
                c_p = std::max(sign_phi * (phi_val - ls_arr(i, j, k - 1, 0)) / DX[2], 0.0);
                c_m = std::min(sign_phi * (ls_arr(i, j, k + 1, 0) - phi_val) / DX[2], 0.0);
                #endif

                // Compute gradient magnitude
                Set::Scalar nablaG = std::sqrt(
                    std::max(a_p * a_p, a_m * a_m) +
                    std::max(b_p * b_p, b_m * b_m) +
                    std::max(c_p * c_p, c_m * c_m)
                );

                // **Reinitialize Level Set Function**
                Set::Scalar phi_new = phi_val - tau * sign_phi * (nablaG - 1.0);
                if (std::abs(phi_new - phi_val) > reinit_tolerance) converged = false;
                ls_arr(i, j, k, 0) = phi_new;
                
                // Debug
                printf("i: %d, j: %d\n", i, j);
                printf("phi(i,j): %f, phi(i-1,j): %f, phi(i+1,j): %f, phi(i,j-1): %f, phi(i,j+1): %f\n",
                    phi_val, ls_arr(i-1,j,k,0), ls_arr(i+1,j,k,0), ls_arr(i,j-1,k,0), ls_arr(i,j+1,k,0));
                printf("a_p: %f, a_m: %f, b_p: %f, b_m: %f c_p: %f, c_m: %f\n", a_p, a_m, b_p, b_m, c_p, c_m); 
                printf("nablaG: %f, sign_phi: %f, phi_new %f\n", nablaG, sign_phi, phi_new);
            });
        }

        // Check convergence across processors
        converged = d_converged.dataValue();
        if (converged) break;

        // Apply BCs to ghost cells
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);
    }
}*/
 
void NarrowBandLevelset_Backup::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], + DX[1] * DX[1], + DX[2] * DX[2]));
    Set::Scalar refinement_threshold = 10.0; // Set refinement threshold to 10

    //for (amrex::MFIter mfi(*Pressure_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (amrex::MFIter mfi(*ls_mf[lev], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<Set::Scalar> const& levelset = (*ls_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=](int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(levelset, i, j, k, 0, DX);
            Set::Scalar grad_magnitude = grad.lpNorm<2>();

            if (grad_magnitude * dr > refinement_threshold) {
                tags(i, j, k) = amrex::TagBox::SET;
            }
        });
    }
}

void NarrowBandLevelset_Backup::Regrid(int lev, Set::Scalar time) {

}

void NarrowBandLevelset_Backup::ApplyBoundaryConditions(int lev, Set::Scalar time) {

    Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);        
    //Integrator::ApplyPatch(lev, time, ls_old_mf, *ls_old_mf[lev], *bc_ls, 0); 
    Integrator::ApplyPatch(lev, time, iNarrowBandMask_mf, *iNarrowBandMask_mf[lev], *bc_ls, 0);
        
    //Integrator::ApplyPatch(lev, time, velocity_mf, *velocity_mf[lev], *bc_velocity, 0);        
 
    //Integrator::ApplyPatch(lev, time, XFlux_mf, *XFlux_mf[lev], bc_nothing, 0);        
    //Integrator::ApplyPatch(lev, time, YFlux_mf, *YFlux_mf[lev], bc_nothing, 0);
#if AMREX_SPACEDIM == 3        
    //Integrator::ApplyPatch(lev, time, ZFlux_mf, *ZFlux_mf[lev], bc_nothing, 0);        
#endif

}

} // namespace Integrator
