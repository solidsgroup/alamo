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
    
    // Following section will be deleted when velocity is taken from ScimitarX integrator class
    {// Manually parse ic.velocity values
    std::vector<amrex::Real> velocity_values;
    pp_queryarr_required("ic.velocity.value", velocity_values);
    
    // Trim extra dimensions beyond AMREX_SPACEDIM
    velocity_values.resize(AMREX_SPACEDIM);
    
    // Initialize velocity IC manually
    value.ic_velocity = new IC::Constant(value.geom, velocity_values);
    
    // Assume Neumann BC for levelset field
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
    // Initialize levelset field
    ic_ls->Initialize(lev, ls_old_mf);
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);
    
    // Initialize velocity field -- will be deleted
    ic_velocity->Initialize(lev, velocity_mf); 
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
    //ComputeAndSetNewTimeStep(); // Compute dt based on global `minDt   
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
        //std::cout << "ls value is: " << ls_mf[lev] << std::endl;
        //std::cout << " velocity value is: " << velocity_mf[lev] << std::endl;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velocity_mf[lev] != nullptr, "velocity_mf[lev] is uninitialized!");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velocity_mf[lev]->ok(), "velocity_mf[lev] MultiFab is not allocated!");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(grids[lev].size() > 0, "Empty BoxArray detected!");
    
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
    Advect(lev, dt);
    
    //Reinitialize the level set function
    if(current_timestep % 100 == 0) {
        Reinitialize(lev);
    
    }
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar dt)
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

void NarrowBandLevelset::Reinitialize(int lev)
{
    const Set::Scalar reinit_tolerance = 1e-3; // Tolerance for stopping criteria
    const int max_iterations = 50;             // Maximum number of reinitialization iterations
    const Set::Scalar epsilon = 1e-6;          // Small value to prevent division by zero

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        bool converged = true;

        // Iterate over all the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);

            amrex::ParallelFor(bx, [=, &converged] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the gradient magnitude |∇phi|
                Set::Vector grad = Numeric::Gradient(ls, i, j, k, 0, geom[lev].CellSize());
                Set::Scalar grad_mag = std::max(grad.lpNorm<2>(), epsilon);  // Prevent division by zero

                // Update the level set function phi using a smoothed sign function
                Set::Scalar sign_phi = ls(i, j, k, 0) / std::sqrt(ls(i, j, k, 0) * ls(i, j, k, 0) + epsilon);
                Set::Scalar phi_new = ls(i, j, k, 0) - sign_phi * (grad_mag - 1.0);
                
                if (std::abs(phi_new - ls(i, j, k, 0)) > reinit_tolerance)
                    converged = false;

                ls(i, j, k, 0) = phi_new;
            });
        }

        // Stop reinitialization if converged
        if (converged)
            break;
    }
}
    
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
