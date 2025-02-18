// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
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
    pp.queryclass(*this); // Call the static Parse function
}

// Define Parse function
void NarrowBandLevelset::Parse(NarrowBandLevelset& value, IO::ParmParse& pp){
    {// Define initial and boundary conditions
    
    // Query the IC assuming either LS::Sphere or LS::Zalesak
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak>("ic.ls",value.ls_ic,value.geom);
    
    // Assume Neumann BC for levelset field
    value.ls_bc = new BC::Constant(value.number_of_components, pp, "bc.ls");
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    }
    
    // Following section will be deleted when velocity is taken from ScimitarX integrator class
    {// Manually parse ic.velocity values
    std::vector<amrex::Real> velocity_values;
    pp_queryarr_required("ic.velocity.value", velocity_values);
    
    // Trim extra dimensions beyond AMREX_SPACEDIM
    velocity_values.resize(AMREX_SPACEDIM);
    
    // Initialize velocity IC manually
    value.velocity_ic = new IC::Constant(value.geom, velocity_values);
    
    // Assume Neumann BC for levelset field
    value.velocity_bc = new BC::Constant(AMREX_SPACEDIM, pp, "bc.velocity");
    
    // Define velocity multifab
    value.RegisterNewFab(value.velocity_mf, value.velocity_bc, AMREX_SPACEDIM, value.number_of_ghost_cells, "Velocity", true); // "true" for debug
    }
}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset field
    ls_ic->Initialize(lev, ls_old_mf);
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);
    
    // Initialize velocity field -- will be deleted
    velocity_ic->Initialize(lev, velocity_mf);
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
     
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) 
{
    // Advance time step
    current_timestep++;
    
    // Advect
    Advect(lev, dt);
    
    //Reinitialize the level set function
    if(current_timestep % 5 == 0) {
        Reinitialize(lev);
    }
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar dt)
{
    // Swap the old ls fab and the new ls fab so we use
    // the new one.
    std::swap(*ls_mf[lev], *ls_old_mf[lev]);
    
    /*// Iterate over all of the patches on this level
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Get the box (index dimensions) for this patch
        const amrex::Box& bx = mfi.tilebox();

        // Get an array-accessible handle to the data on this patch.
        Set::Patch<const Set::Scalar>  ls_old   = ls_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>        ls       = ls_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>        velocity = ls_mf.Patch(lev,mfi);
    }*/
    
    // Update the velocity
    UpdateInterfaceVelocity(lev);
    
    // Compute the flux
    
   
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
                Set::Scalar sign_phi = ls(i, j, k) / std::sqrt(ls(i, j, k) * ls(i, j, k) + epsilon);
                Set::Scalar phi_new = ls(i, j, k) - sign_phi * (grad_mag - 1.0);
                
                if (std::abs(phi_new - ls(i, j, k)) > reinit_tolerance)
                    converged = false;

                ls(i, j, k) = phi_new;
            });
        }

        // Stop reinitialization if converged
        if (converged)
            break;
    }
}
    
void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{

}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {

}
} // namespace Integrator
