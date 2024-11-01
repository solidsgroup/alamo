
#ifndef INTEGRATOR_SCIMITARXLEVELSET_H // Include guards
#define INTEGRATOR_SCIMITARXLEVELSET_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// AMReX Includes
#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

// Alamo Includes
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "IC/IC.H"
#include "IC/SphereLS.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
class ScimitarxLevelSet : virtual public Integrator
{
public:
    // Empty constructor
    ScimitarxLevelSet(int a_nghost = 2) : 
        Integrator(),
        number_of_ghost_cells(a_nghost)
    {}

    // Constructor that triggers parse
    ScimitarxLevelSet(IO::ParmParse& pp) : ScimitarxLevelSet()
    {
        Parse(*this, pp);
    }

    virtual ~ScimitarxLevelSet()
    {
        delete ls_ic;
        delete ls_bc;
    }

    // The Parse function initializes the HeatConduction object using
    // a parser, pp. 
    // Note that this is a static function, which means it does not have
    // direct access to member variables. Instead, it initializes the variables
    // inside the argument, "value", and so all references to member items are
    // prefixed by "value."
    static void Parse(ScimitarxLevelSet& value, IO::ParmParse& pp)
    {
        // Diffusion coefficient :math:`\alpha`.
        //   *This is an example of a required input variable -
        //    - program will terminate unless it is provided.*
        pp_query_required("heat.alpha", value.alpha);

        // Criterion for mesh refinement.
        //   *This is an example of a default input variable.
        //    The default value is provided here, not in the 
        //    definition of the variable.*
        pp_query_default("heat.refinement_threshold", value.refinement_threshold, 0.01);

        std::string type;

        // Initial condition type.
        //   *This is an example of type validation:
        //    the input variable must be one of the three provided values, 
        //    and will error if not. The default selection is the first
        //    argument.*
        pp_query_validate("ls.ic.type", type, {"constant", "spherels","expression"});
        if (type == "spherels")          value.ls_ic = new IC::SphereLS(value.geom, pp, "ic.spherels");
        else  Util::Abort(INFO, "Invalid ic.type ", type);

        std::string bc_type;
        // Select BC object for temperature
        pp_query_validate("ls.bc.type",bc_type,{"constant","expression"});
        if (bc_type == "expression")      value.ls_bc = new BC::Expression(1, pp, "ls.bc.expression");
        else if (bc_type == "constant")   value.ls_bc = new BC::Constant(1, pp, "ls.bc");

        // Register the temperature and old temperature fields.
        // temp_mf and temp_old_mf are defined near the bottom of this Header file.
        value.RegisterNewFab(value.ls_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS", true);
        value.RegisterNewFab(value.ls_old_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS_old", true);
        value.RegisterNewFab(value.lsflux_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LSFlux", false);
    }

protected:

    // Use the ic object to initialize the temperature field
    void Initialize(int lev)
    {
        ls_ic->Initialize(lev, ls_mf);
    }

    // Integrate the heat equation
    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
           std::swap(*ls_mf[lev], *ls_old_mf[lev]);

        // Get the cell size corresponding to this level
        const Set::Scalar* DX = geom[lev].CellSize();
	/*const Set::Scalar Narrow_Band_Width = 6.0 * DX[0];
        const Set::Scalar InnerTube = -Narrow_Band_Width;
        const Set::Scalar OuterTube = Narrow_Band_Width;
        const Set::Scalar radius    = 0.25;*/	

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box& bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            // amrex::Array4<const Set::Scalar> const& temp_old = (*temp_old_mf[lev]).array(mfi);
            // amrex::Array4<Set::Scalar>       const& temp = temp_mf.Patch(lev,mfi);//(*temp_mf[lev]).array(mfi);
            Set::Patch<const Set::Scalar>        ls_old = ls_old_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar>        ls     = ls_mf.Patch(lev,mfi);
	    Set::Patch<Set::Scalar>        flux   = lsflux_mf.Patch(lev,mfi);

            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
	     /*Set::Scalar  x = (i + 0.5) * DX[0];
	     Set::Scalar  y = (j + 0.5) * DX[1];
	     Set::Scalar  z = (AMREX_SPACEDIM == 3) ? (k + 0.5) * DX[2] : 0.0;
	     Set::Scalar  dist = std::sqrt((x - 0.5) * (x - 0.5) +
                                      (y - 0.5) * (y - 0.5) +
                                      (z - 0.5) * (z - 0.5));*/
                Set::Scalar ls_xm = (i > 0) ? ls_old(i-1, j, k) : ls_old(i, j, k);
                Set::Scalar ls_xp = ls_old(i, j, k);
                Set::Scalar vel_x = 0.1;
                Set::Scalar flux_x = (vel_x > 0.0) ? vel_x * ls_xm : vel_x * ls_xp;

                Set::Scalar ls_ym = (j > 0) ? ls_old(i, j-1, k) : ls_old(i, j, k);
                Set::Scalar ls_yp = ls_old(i, j, k);
                Set::Scalar vel_y = 0.0;
                Set::Scalar flux_y = (vel_y > 0.0) ? vel_y * ls_ym : vel_y * ls_yp;

                Set::Scalar flux_z = 0.0;
                if (AMREX_SPACEDIM == 3)
                {
                    Set::Scalar ls_zm = (k > 0) ? ls_old(i, j, k-1) : ls_old(i, j, k);
                    Set::Scalar ls_zp = ls_old(i, j, k);
                    Set::Scalar vel_z = 0.0;
                    flux_z = (vel_z > 0.0) ? vel_z * ls_zm : vel_z * ls_zp;
                }

                flux(i, j, k) = -(flux_x / DX[0] + flux_y / DX[1] + flux_z / DX[2]);

	        ls(i, j, k) = ls_old(i, j, k) + dt*flux(i, j, k);     
            });
        }
    }

    // Tag cells for mesh refinement based on temperature gradient
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
    {
        // Get cell dimensions as done above.
        const Set::Scalar* DX = geom[lev].CellSize();
        // Calculate the diagonal.
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Iterate over the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char>         const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const& ls = (*ls_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the temperature gradient.
                Set::Vector grad = Numeric::Gradient(ls, i, j, k, 0, DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad.lpNorm<2>() * dr > 0.0)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<Set::Scalar> ls_mf;         // Levelset field variable (current timestep)
    Set::Field<Set::Scalar> ls_old_mf;     // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar> lsflux_mf;     // Levelset flux field variable 

private:

    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    const int number_of_components = 1;      // Number of components
    const int number_of_ghost_cells = 2;     // Number of ghost cells


    Set::Scalar alpha = NAN;                 // Thermal diffusivity
    Set::Scalar refinement_threshold = NAN ; // Criterion for cell refinement

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC* ls_ic = nullptr;                    // Object used to initialize temperature field
    BC::BC<Set::Scalar>* ls_bc = nullptr;       // Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
