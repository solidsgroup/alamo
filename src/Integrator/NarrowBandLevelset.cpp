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
    
    // Define narrowband multifab to track tube
    value.RegisterNewFab(value.iNarrowBandMask_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "NBand", true);
    }
    
    {// Following section will be deleted when velocity is taken from ScimitarX integrator class
    
    // Initialize velocity IC and BC 
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
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    ic_ls->Initialize(lev, ls_mf);
    
    // Initialize narrowband
    UpdateNarrowBand(lev);
    
    // Initialize velocity field -- will be deleted
    ic_velocity->Initialize(lev, velocity_mf); 
}

void NarrowBandLevelset::UpdateNarrowBand(int lev) {
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar narrow_band_width = 6.0 * min_DX;
    const Set::Scalar inner_tube = -narrow_band_width;
    const Set::Scalar outer_tube = narrow_band_width;

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& ls_arr = ls_mf.Patch(lev, mfi);
        auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar phi = ls_arr(i, j, k);
            phi = std::clamp(phi, inner_tube, outer_tube); // Set innertube <= phi <= outertube
            ls_arr(i, j, k) = phi;

            nbmask_arr(i, j, k) = (std::abs(phi) <= min_DX) ? 0 : (phi > 0 ? (phi < outer_tube ? 1 : 2) : (phi > inner_tube ? -1 : -2));
        });
    }
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
            auto const& velocity_arr = velocity_mf.Patch(lev, mfi);  // Velocity field

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
    
    /*printf("After Advect: \n");
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& ghost_bx = mfi.growntilebox(1);//ls_mf[lev]->nGrow());
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> nb = iNarrowBandMask_mf.Patch(lev, mfi);

        amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            printf("ls(%d, %d, %d) = %f\n", i, j, k, ls(i,j,k,0));
            printf("nb(%d, %d, %d) = %f\n", i, j, k, nb(i,j,k,0));
        });
    }*/
    
    // Reinitialize the level set function
    if(current_timestep % 100 == 0) {
        Reinitialize(lev, time);
    }     
    
    // Update narrowband after reinitialization
    UpdateNarrowBand(lev);  
    
    

    /*printf("After Reinitialize: \n");
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& ghost_bx = mfi.growntilebox(1); //ls_mf[lev]->nGrow());
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> nb = iNarrowBandMask_mf.Patch(lev, mfi);

        amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            printf("ls(%d, %d, %d) = %f\n", i, j, k, ls(i,j,k,0));
            printf("nb(%d, %d, %d) = %f\n", i, j, k, nb(i,j,k,0));
        });
    }*/
    
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

void NarrowBandLevelset::Reinitialize(int lev, Set::Scalar time) {
    // Loop constants
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50; 
    
    // Grid spacing
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    
    // Time constraint for PDE
    const Set::Scalar tau = cflNumber * min_DX;
    
    // Narrowband constants
    const Set::Scalar Narrow_Band_Width = 6.0 * min_DX;
    const Set::Scalar inner_tube = -Narrow_Band_Width;
    const Set::Scalar outer_tube = Narrow_Band_Width;
    
    // Update the levelset value after the flux computation
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& ls_arr = ls_mf.Patch(lev, mfi);
        auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar mask_val = nbmask_arr(i, j, k);
            if (mask_val == 0) {
                //ls_arr(i, j, k) = 0.0;
            } else if (mask_val == 2 || mask_val == -2) {
                ls_arr(i, j, k) = (mask_val < 0) ? inner_tube : outer_tube;
            }
        });
    }
    
    /*//  Create duplicate to update levelset in place
    amrex::MultiFab ls_reinit(*ls_mf[lev], amrex::make_alias, 0, 1); // Alias, so we update in place

    // Main loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        printf("iter: %d\n", iter);
        bool converged = true;
        amrex::Gpu::DeviceScalar<bool> d_converged(true);
        bool* d_converged_ptr = d_converged.dataPtr(); // Use dataPtr() to get device-accessible pointer
        
        for (amrex::MFIter mfi(ls_reinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& valid_bx = mfi.tilebox();
            auto const& ls_arr = ls_mf.Patch(lev, mfi);
            auto const& nbmask_arr = iNarrowBandMask_mf.Patch(lev, mfi);
            amrex::ParallelFor(valid_bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                // Skip non-narrow band cells
                if (nbmask_arr(i, j, k) == 2 || nbmask_arr(i, j, k) == -2) return;
                
                // Compute first order differences
                Set::Scalar phi_val = ls_arr(i, j, k, 0);
                Set::Scalar sign_phi = std::copysign(1.0, phi_val);
                
                // Compute first-order upwind gradients with branch-free min/max logic
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

                // Compute gradient magnitude without branches
                Set::Scalar nablaG = std::sqrt(
                    std::max(a_p * a_p, a_m * a_m) +
                    std::max(b_p * b_p, b_m * b_m) +
                    std::max(c_p * c_p, c_m * c_m)
                );

                // Update phi_new
                Set::Scalar phi_new = phi_val - tau * sign_phi * (nablaG - 1.0);
                if (std::abs(phi_new - phi_val) > reinit_tolerance) converged = false;
                ls_arr(i, j, k, 0) = phi_new;
                
                // Check converged variable across processors
                if (std::abs(phi_new - phi_val) > reinit_tolerance) *d_converged_ptr = false;
                
                // Debug
                /*printf("i: %d, j: %d\n", i, j);
                printf("phi(i,j): %f, phi(i-1,j): %f, phi(i+1,j): %f, phi(i,j-1): %f, phi(i,j+1): %f\n",
                    phi_val, ls_arr(i-1,j,k,0), ls_arr(i+1,j,k,0), ls_arr(i,j-1,k,0), ls_arr(i,j+1,k,0));
                printf("a_p: %f, a_m: %f, b_p: %f, b_m: %f c_p: %f, c_m: %f\n", a_p, a_m, b_p, b_m, c_p, c_m); 
                printf("nablaG: %f, sign_phi: %f, phi_new %f\n", nablaG, sign_phi, phi_new);
            });
        }
        
        // Check convergence across processors
        converged = d_converged.dataValue(); // Fix: Use dataValue()
        if (converged) break;
        
        // Apply BCs to ghost cells
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0); 
    }*/
}  
void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
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

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {

}

void NarrowBandLevelset::ApplyBoundaryConditions(int lev, Set::Scalar time) {

    Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);        
    Integrator::ApplyPatch(lev, time, ls_old_mf, *ls_old_mf[lev], *bc_ls, 0); 
    Integrator::ApplyPatch(lev, time, iNarrowBandMask_mf, *iNarrowBandMask_mf[lev], *bc_ls, 0);
        
    Integrator::ApplyPatch(lev, time, velocity_mf, *velocity_mf[lev], *bc_velocity, 0);        
 
    Integrator::ApplyPatch(lev, time, XFlux_mf, *XFlux_mf[lev], bc_nothing, 0);        
    Integrator::ApplyPatch(lev, time, YFlux_mf, *YFlux_mf[lev], bc_nothing, 0);
#if AMREX_SPACEDIM == 3        
    Integrator::ApplyPatch(lev, time, ZFlux_mf, *ZFlux_mf[lev], bc_nothing, 0);        
#endif

}

} // namespace Integrator
