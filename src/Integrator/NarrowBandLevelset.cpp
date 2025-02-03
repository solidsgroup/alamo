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
    pp.select<IC::LS::Sphere,IC::LS::Zalesak>("ls.ic",value.ls_ic,value.geom);
    
    // Assume Neumann BC
    value.ls_bc = new BC::Constant(1, pp, "ls.bc");
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    }

}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    ls_ic->Initialize(lev, ls_old_mf);
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
     
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) {

}

void Reinitialize(int lev)
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
