// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/TimeStepper.H"
#include "Integrator/NarrowBandLevelset.H"

// BC
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "BC/Expression.H"

// IC
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"
#include "IC/Constant.H"

// Numeric
#include "Util/Util.H"

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
    
    {// Following section will be deleted when velocity is taken from ScimitarX integrator class
    
    // Initialize velocity IC and BC 
    value.ic_velocity = new IC::Constant(value.geom, pp, "ic.velocity");
    value.bc_velocity = new BC::Constant(AMREX_SPACEDIM, pp, "bc.velocity");
    
    // Define velocity multifab
    value.RegisterNewFab(value.velocity_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "Velocity", true); // "true" for debug
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
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    ic_ls->Initialize(lev, ls_mf);
    
    // Initialize velocity field -- will be deleted
    ic_velocity->Initialize(lev, velocity_mf); 
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
    ComputeAndSetNewTimeStep(); // Compute dt based on global `minDt   
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

    for (int lev = 0; lev <= maxLevel(); ++lev) {  // Iterate over AMR levels
        const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`

        for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
            auto const& velocity_arr = velocity_mf[lev]->array(mfi);  // Velocity field

            Set::Scalar minDt_local = std::numeric_limits<Set::Scalar>::max();  // Thread-local minDt

            amrex::ParallelFor(bx, [=, &minDt_local](int i, int j, int k) noexcept {
                // Extract velocity components
                Set::Scalar u = velocity_arr(i, j, k, 0);
#if (AMREX_SPACEDIM >= 2)
                Set::Scalar v = velocity_arr(i, j, k, 1);
#endif
#if (AMREX_SPACEDIM == 3)
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
            });

            // Update the global minDt
            minDt = std::min(minDt, minDt_local);
        }
    }

    // Reduce across MPI processes to find the global minimum timestep
    amrex::ParallelDescriptor::ReduceRealMin(minDt);
    
    return cflNumber * minDt;  // Return CFL-adjusted time step
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) 
{
    // Advance time step
    current_timestep++;
    
    // Apply boundary conditions
    ApplyBoundaryConditions(lev, time);
    
    // Advect
    Advect(lev, time, dt);
    
    /*Set::Scalar max_phi = 0;
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box valid_bx = mfi.tilebox();
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (ls(i,j,k,0) > 0)
            max_phi = ls(i,j,k,0);
            //printf("End timestep ls(%d, %d, %d) = %f\n", i, j, k, ls(i,j,k,0));
        });
    }
    printf("max_phi: %f\n", max_phi);*/
    
    // Reinitialize the level set function
    if(current_timestep % 1 == 0) {
        Reinitialize(lev, time);
    }
    
    // Redefine narrowband
    ApplyNarrowBanding(lev);
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar time, Set::Scalar dt)
{
    // Swap the old ls fab and the new ls fab so we use
    // the new one.
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);

    // Update the velocity
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

                ApplyBoundaryConditions(lev, time);

                // 3. Compute sub-step using the chosen time-stepping scheme
                timeStepper->ComputeSubStep(lev, dt, stage, this);
            }
            break;
        }
    
    } 
}

void NarrowBandLevelset::UpdateInterfaceVelocity(int lev)
{
    // Loop through and set all velocity components
    for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar>  velocity = velocity_mf.Patch(lev,mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Do computation here - right now placeholder to keep velocity constant
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                velocity(i,j,k,d) = velocity(i,j,k,d);
            }
        });
    }
}

void NarrowBandLevelset::Reinitialize(int lev, Set::Scalar time)
{
    const Set::Scalar reinit_tolerance = 1e-3;     
    const int max_iterations = 50;                 
    const Set::Scalar epsilon = 1e-6;              
    const Set::Scalar* DX = geom[lev].CellSize();
    printf("DX: (%f, %f)\n", DX[0], DX[1]);
    Set::Scalar dt = 0.5 * DX[0];  

    // **Tag narrow-band cells (1 = update, 0 = frozen)**
    /*amrex::iMultiFab tag_mask(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, 0);

    for (amrex::MFIter mfi(tag_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& valid_bx = mfi.tilebox();
        auto const& ls_arr = ls_mf[lev]->array(mfi);
        auto const& tag_arr = tag_mask.array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar phi_val = ls_arr(i, j, k, 0);
            tag_arr(i, j, k) = (std::abs(phi_val) < 6.0 * DX[0]) ? 1 : 0;
        });
    }*/
    
    // Print before copy
    /*printf("before copy: \n");
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box valid_bx = mfi.tilebox();
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> ls_old = ls_old_mf.Patch(lev, mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            printf("i:%d, j: %d, ls: %f, ls_old: %f\n", i, j, ls(i,j,k,0), ls_old(i,j,k,0)); 
        });
    }*/
    
    // Copy the old/new levelset values before reinitializing
    amrex::MultiFab::Copy(*ls_old_mf[lev], *ls_mf[lev], 0, 0, ls_mf[lev]->nComp(), ls_mf[lev]->nGrow());
    
    // Apply boundary conditions
    Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);        
    Integrator::ApplyPatch(lev, time, ls_old_mf, *ls_old_mf[lev], *bc_ls, 0);
    
    // **Reinitialize only tagged narrow-band cells**
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        bool converged = true;
        //printf("iter: %d\n", iter);

        // Main reinitialization loop
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& valid_bx = mfi.tilebox();
            //const amrex::Box& ghost_bx = mfi.growntilebox(ls_mf[lev]->nGrow());
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& ls_old_arr = ls_old_mf.Patch(lev, mfi); 
            //auto const& tag_arr = tag_mask.array(mfi);

            amrex::ParallelFor(valid_bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                //if (tag_arr(i, j, k) == 0)
                    //return;

                Set::Scalar phi_val = ls_old_arr(i, j, k, 0);  // Use old value for gradient calculation
                auto stencil = Numeric::GetStencil(i, j, k, valid_bx);
                Set::Vector grad = Numeric::Gradient(ls_old_arr, i, j, k, 0, DX, stencil);  // Use ls_old for gradient
                Set::Scalar grad_mag = std::max(grad.lpNorm<2>(), epsilon);
                
                Set::Scalar sign_phi = phi_val / std::sqrt(phi_val * phi_val + epsilon);

                Set::Scalar phi_new = phi_val - dt * sign_phi * (grad_mag - 1.0);
                
                // **Fix: Add debug print**
                /*if (std::abs(phi_val) < 1e-2)
                {
                    printf("DEBUG: Iter=%d, i=%d, j=%d, phi=%.6f, sign_phi=%.6f, grad_mag=%.6f, phi_new=%.6f\n",
                           iter, i, j, phi_val, sign_phi, grad_mag, phi_new);
                }*/
                
                //printf("i: %d, j: %d, grad: (%f %f), grad_mag: %f, phi_old: %f, phi_new: %f\n", i, j, grad[0], grad[1], grad_mag, phi_val, phi_new);
                //printf("Iteration %d, i: %d, j: %d, phi_old: %.6f, phi_new: %.6f, diff: %.6f\n", 
                    //iter, i, j, phi_val, phi_new, std::abs(phi_new - phi_val));


                if (std::abs(phi_new - phi_val) > reinit_tolerance)
                    converged = false;

                ls_arr(i, j, k, 0) = phi_new;  // Update the level set with the new value
                
                //printf("i: %d, j: %d, ls_arr: %f, ls_old_arr: %f\n", i, j, ls_arr(i,j,k,0), ls_old_arr(i,j,k,0));
            });
        }

        if (converged)
            break;
            
        // Swap the old/new levelset each iteration
        std::swap(*ls_mf[lev], *ls_old_mf[lev]);
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);        
        Integrator::ApplyPatch(lev, time, ls_old_mf, *ls_old_mf[lev], *bc_ls, 0);
    }
}

/*void NarrowBandLevelset::Reinitialize(int lev)
{
    const Set::Scalar reinit_tolerance = 1e-3;     
    const int max_iterations = 1;                 
    const Set::Scalar epsilon = 1e-6;              
    const Set::Scalar* DX = geom[lev].CellSize();
    printf("DX: (%f, %f)\n", DX[0], DX[1]);
    Set::Scalar dt = 0.5 * DX[0];  

    // **Tag narrow-band cells (1 = update, 0 = frozen)**
    amrex::iMultiFab tag_mask(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, 0);

    for (amrex::MFIter mfi(tag_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& valid_bx = mfi.tilebox();
        auto const& ls_arr = ls_mf[lev]->array(mfi);
        auto const& tag_arr = tag_mask.array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar phi_val = ls_arr(i, j, k, 0);
            tag_arr(i, j, k) = (std::abs(phi_val) < 6.0 * DX[0]) ? 1 : 0;
        });
    }
    
    printf("before re: \n");
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box valid_bx = mfi.tilebox();
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            printf("i:%d, j: %d, ls: %f\n", i, j, ls(i,j,k,0)); 
        });
    }

    // **Reinitialize only tagged narrow-band cells**
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        bool converged = true;
        //printf("iter: %d\n", iter);

        amrex::MultiFab ls_temp(*ls_mf[lev], amrex::make_alias, 0, 1);  // Read-only buffer for old values

        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& valid_bx = mfi.tilebox();
            auto const& ls_old_arr = ls_temp.array(mfi);
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& tag_arr = tag_mask.array(mfi);

            amrex::ParallelFor(valid_bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (tag_arr(i, j, k) == 0)
                    return;

                Set::Scalar phi_val = ls_old_arr(i, j, k, 0);
                auto stencil = Numeric::GetStencil(i, j, k, valid_bx);
                Set::Vector grad = Numeric::Gradient(ls_old_arr, i, j, k, 0, DX, stencil);
                Set::Scalar grad_mag = std::max(grad.lpNorm<2>(), epsilon);
                
                Set::Scalar sign_phi = (phi_val > 0) ? 1.0 : (phi_val < 0) ? -1.0 : 0.0;
                Set::Scalar sign_phi = (std::abs(phi_val) > epsilon) ? 
                                       (phi_val / std::sqrt(phi_val * phi_val + epsilon)) : 
                                       0.0;
                                       
                // Add Debugging Prints
                if (std::abs(phi_val) < 1e-2)  // Near zero level set
                {
                    printf("DEBUG: Iter=%d, i=%d, j=%d, phi=%.6f, sign_phi=%.6f, grad_mag=%.6f\n", 
                           iter, i, j, phi_val, sign_phi, grad_mag);
                }                  

                
                Set::Scalar phi_new = phi_val - dt * sign_phi * (grad_mag - 1.0);
                
                // DEBUG
                //printf("i: %d, j: %d\n", i, j);
                //printf("ls(i,j): %f, grad: (%f, %f), grad_mag: %f, sign_phi: %f\n", phi_val, grad[0], grad[1], grad_mag, sign_phi);
                //printf("phi_new: %f, phi_old: %f\n", phi_new, phi_val);

                if (std::abs(phi_new - phi_val) > reinit_tolerance)
                    converged = false;

                ls_arr(i, j, k, 0) = phi_new;
                printf("i: %d, j: %d, ls_arr: %f, temp_array: %f\n", i, j, ls_arr(i,j,k,0), ls_old_arr(i,j,k,0));
            });
        }

        if (converged)
            break;
            
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box valid_bx = mfi.tilebox();
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar> temp = ls_temp.array(mfi);

            amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                printf("i:%d, j: %d, ls: %f, ls_old: %f\n", i, j, ls(i,j,k,0), temp(i,j,k,0)); 
            });
        }
    }
}*/

void NarrowBandLevelset::ApplyNarrowBanding(int lev)
{
    const Set::Scalar* DX = geom[lev].CellSize();  // Cell size array for finite difference operations
    const Set::Scalar Narrow_Band_Width = 6.0 * DX[0];  
    const Set::Scalar InnerTube = -Narrow_Band_Width;
    const Set::Scalar OuterTube = Narrow_Band_Width;

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& valid_bx = mfi.tilebox();
        auto const& ls_arr = ls_mf[lev]->array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar phi_val = ls_arr(i, j, k, 0);
            ls_arr(i, j, k, 0) = (std::abs(phi_val) <= Narrow_Band_Width) ? phi_val : (phi_val < 0 ? InnerTube : OuterTube);
        });
    }
}

/*void NarrowBandLevelset::Reinitialize(int lev)
{
    const Set::Scalar reinit_tolerance = 1e-3;     // Tolerance for stopping criteria
    const int max_iterations = 50;                 // Maximum number of reinitialization iterations
    const Set::Scalar epsilon = 1e-6;              // Small value to prevent division by zero
    const Set::Scalar* DX = geom[lev].CellSize();  // Access the geometry at level `lev`
    const Set::Scalar Narrow_Band_Width = 6.0 * DX[0];
    const Set::Scalar InnerTube = -Narrow_Band_Width;
    const Set::Scalar OuterTube = Narrow_Band_Width;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        bool converged = true;

        // Iterate over all the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& valid_bx = mfi.tilebox();
            Set::Patch<Set::Scalar> ls_old = ls_old_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);

            amrex::ParallelFor(valid_bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Get lower/upper bounds to define proper stencil for gradient
                auto stencil = Numeric::GetStencil(i, j, k, valid_bx);
                //if (i==0 && j==0){
                //printf("At i= %d, j= %d, Stencil: %d %d\n", i, j, stencil[0], stencil[1]);
                //}
            
                // Default stencil: central difference
                std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType;
                
                // Set stencils dynamically based on valid_bx boundaries
                if (i == valid_bx.smallEnd(0)) {
                    stencil[0] = Numeric::StencilType::Hi;  // Forward difference at left boundary
                } 
                else if (i == valid_bx.bigEnd(0)) {
                    stencil[0] = Numeric::StencilType::Lo;  // Backward difference at right boundary
                }

                #if AMREX_SPACEDIM > 1
                if (j == valid_bx.smallEnd(1)) {
                    stencil[1] = Numeric::StencilType::Hi;  // Forward difference at bottom boundary
                } 
                else if (j == valid_bx.bigEnd(1)) {
                    stencil[1] = Numeric::StencilType::Lo;  // Backward difference at top boundary
                }
                #endif

                #if AMREX_SPACEDIM > 2
                if (k == valid_bx.smallEnd(2)) {
                    stencil[2] = Numeric::StencilType::Hi;  // Forward difference at front boundary
                } 
                else if (k == valid_bx.bigEnd(2)) {
                    stencil[2] = Numeric::StencilType::Lo;  // Backward difference at back boundary
                }
                #endif
                
                
                // Calculate the gradient magnitude |∇phi|
                if (i==1 && j==0){
                    printf("phi(i,j): %f\n", ls(i,j,k,0));
                    printf("phi(i+1,j): %f\n", ls(i+1,j,k,0));
                    printf("phi(i-1,j): %f\n", ls(i-1,j,k,0));
                    printf("phi(i,j+1): %f\n", ls(i,j+1,k,0));
                    //printf("phi(i,j-1): %f\n", ls(i,j-1,k,0));
                    printf("DX[0]: %f\n", DX[0]);
                    printf("DX[1]: %f\n", DX[1]);
                }
                Set::Vector grad = Numeric::Gradient(ls_old, i, j, k, 0, DX, stencil);
                /*if (i==0 && j==0){
                    printf("grad[0]: %f\n", grad[0]);
                    printf("grad[1]: %f\n", grad[1]);
                }
                if (std::isinf(grad[0]) || std::isinf(grad[0])){
                    //printf("Grad (%f, %f) is inf at i=%d j=%d\n", grad[0], grad[1], i, j);
                    //printf("Stencil: %d %d\n", stencil[0], stencil[1]);
                }
                Set::Scalar grad_mag = std::max(grad.lpNorm<2>(), epsilon);  // Prevent division by zero

                // Update the level set function phi using a smoothed sign function
                Set::Scalar sign_phi = ls_old(i, j, k, 0) / std::sqrt(ls_old(i, j, k, 0) * ls_old(i, j, k, 0) + epsilon);
                Set::Scalar phi_new = ls_old(i, j, k, 0) - sign_phi * (grad_mag - 1.0);
                if (std::isinf(phi_new)){
                    printf("Phi_new is inf at i=%d j=%d\n", i, j);
                }*
                
                if (std::abs(phi_new - ls_old(i, j, k, 0)) > reinit_tolerance)
                    converged = false;

                //ls(i, j, k, 0) = phi_new;
                ls(i, j, k, 0) = (std::abs(phi_new) <= Narrow_Band_Width) ? phi_new : (phi_new < 0 ? InnerTube : OuterTube);

                
                
                // Debug
                //printf("i: %d, j: %d, LS: %f, LS_old: %f\n", i, j, ls(i,j,k,0), ls_old(i,j,k,0));
            });
        }

        // Stop reinitialization if converged
        if (converged)
            break;
            
        // Swap old ls value with new one
        std::swap(*ls_old_mf[lev], *ls_mf[lev]);
    }
}*/
    
void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));
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

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {

}

void NarrowBandLevelset::ApplyBoundaryConditions(int lev, Set::Scalar time) {

        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);        
        Integrator::ApplyPatch(lev, time, ls_old_mf, *ls_old_mf[lev], *bc_ls, 0); 
        
        Integrator::ApplyPatch(lev, time, velocity_mf, *velocity_mf[lev], *bc_velocity, 0);        
 
        Integrator::ApplyPatch(lev, time, XFlux_mf, *XFlux_mf[lev], bc_nothing, 0);        
        Integrator::ApplyPatch(lev, time, YFlux_mf, *YFlux_mf[lev], bc_nothing, 0);
#if AMREX_SPACEDIM == 3        
        Integrator::ApplyPatch(lev, time, ZFlux_mf, *ZFlux_mf[lev], bc_nothing, 0);        
#endif

}

} // namespace Integrator
