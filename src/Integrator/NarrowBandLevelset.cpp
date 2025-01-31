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

    
// Constructor that triggers Parse
NarrowBandLevelset::NarrowBandLevelset(IO::ParmParse& pp):NarrowBandLevelset() // Call default constructor
{
    pp.queryclass(*this); // Call the static Parse function
}

// Define Parse function
void NarrowBandLevelset::Parse(NarrowBandLevelset& value, IO::ParmParse& pp){
    { // Query the boundary conditions
    std::string bc_type;
    
    // Select BC object for levelset
    pp.select<BC::Constant,BC::Expression>("ls.bc",value.ls_bc,1);
    // pp_query_validate("ls.bc.type",bc_type,{"constant","expression"});
    // if (bc_type == "expression")      value.ls_bc = new BC::Expression(1, pp, "ls.bc.expression");
    // else if (bc_type == "constant")   value.ls_bc = new BC::Constant(1, pp, "ls.bc");
    }

    {    
    // Register the levelset and old levelset fields
    value.RegisterNewFab(value.ls_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);

    }

    { // Query the levelset initial conditions
    std::string ic_type;
    
    // Validate the initial condition type is a correct string
    pp_query_validate("ls.ic.type", ic_type, {"constant", "spherels","zalesakls","expression"});
    if (ic_type == "spherels")        value.ls_ic = new IC::LS::Sphere(value.geom, pp, "ic.spherels");
    else if (ic_type == "zalesakls")  value.ls_ic = new IC::LS::Zalesak(value.geom, pp, "ic.zalesakls");
    }


}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    ls_ic->Initialize(lev, ls_mf);

/*    // Iterate over all the patches on this level
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> nbmask_patch = iNarrowBandMask_mf.Patch(lev, mfi);


        const Set::Scalar narrow_band_width = 6.0 * geom[lev].CellSize()[0];
        const Set::Scalar inner_tube = -narrow_band_width;
        const Set::Scalar outer_tube = narrow_band_width;

        // Update interface_mf to reflect the level set initialization
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (std::abs(ls(i, j, k)) <= geom[lev].CellSize()[0])
            {
                nbmask_patch(i, j, k) = 0; // Cells closest to the interface (considered zero level set)
            }
            else if (ls(i, j, k) > inner_tube && ls(i, j, k) < outer_tube)
            {
                nbmask_patch(i, j, k) = (ls(i, j, k) > 0) ? 1 : -1; // Narrow band cells
            }
            else
            {
                nbmask_patch(i, j, k) = (ls(i, j, k) > 0) ? 2 : -2; // Cells outside the narrow band
            }
        });
    }*/
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
     
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) {

}

void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{

}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {

}
} // namespace Integrator
