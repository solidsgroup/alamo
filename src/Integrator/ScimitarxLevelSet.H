
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
#include "BC/IntBC.H"
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
	delete nbmask_bc;
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
        // Select BC object for levelset
        pp_query_validate("ls.bc.type",bc_type,{"constant","expression"});
        if (bc_type == "expression")      value.ls_bc = new BC::Expression(1, pp, "ls.bc.expression");
        else if (bc_type == "constant")   value.ls_bc = new BC::Constant(1, pp, "ls.bc");

        // Select BC object for narrowband mask
        pp_query_validate("nbmask.bc.type",bc_type,{"constant"});
        if (bc_type == "constant")   value.nbmask_bc = new BC::NothingInt(1, pp, "nbmask.bc");

        // Register the levelset and old levelset fields.
        value.RegisterNewFab(value.ls_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS", true);
        value.RegisterNewFab(value.ls_old_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
        value.RegisterNewFab(value.lsflux_mf, value.ls_bc, value.number_of_components, value.number_of_ghost_cells, "LSFlux", false);
        value.RegisterNewFab(value.iNarrowBandMask_mf, value.nbmask_bc, value.number_of_components, value.number_of_ghost_cells, "NBand", false);

    }

protected:

// Use the ic object to initialize the level set field
void Initialize(int lev)
{
    ls_ic->Initialize(lev, ls_mf);

    // Iterate over all the patches on this level
    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
        Set::Patch<Set::IntScalar> nbmask_patch = iNarrowBandMask_mf.Patch(lev, mfi);


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
    }
}

    // Integrate the heat equation
    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
    {
            // Increment the timestep
               current_timestep++;

    	    // Advect the level set field
               Advect(lev, dt);

	     //Reinitialize the level set function
	     if(current_timestep % 5 == 0) {
               Reinitialize(lev);
	      }

            // Update the narrow band
               UpdateNarrowBand(lev);
       
    }

    // Advect the level set field
    void Advect(int lev, Set::Scalar dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
           std::swap(*ls_mf[lev], *ls_old_mf[lev]);

   	// Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box& bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch
            Set::Patch<const Set::Scalar> ls_old = ls_old_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar> flux = lsflux_mf.Patch(lev, mfi);

            // Calculate the flux
            CalculateFlux(ls_old, flux, bx, lev);

            // Update the level set function
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                ls(i, j, k) = ls_old(i, j, k) + dt * flux(i, j, k);
            });
        }
    }

    // Calculate flux using an upwind scheme
    void CalculateFlux(const Set::Patch<const Set::Scalar>& ls_old,
                       Set::Patch<Set::Scalar>& flux,
                       const amrex::Box& bx,
                       int lev)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar ls_xm = (i > 0) ? ls_old(i - 1, j, k) : ls_old(i, j, k);
            Set::Scalar ls_xp = ls_old(i, j, k);
            Set::Scalar vel_x = 0.0;
            Set::Scalar flux_x = (vel_x > 0.0) ? vel_x * ls_xm : vel_x * ls_xp;

            Set::Scalar ls_ym = (j > 0) ? ls_old(i, j - 1, k) : ls_old(i, j, k);
            Set::Scalar ls_yp = ls_old(i, j, k);
            Set::Scalar vel_y = 0.0;
            Set::Scalar flux_y = (vel_y > 0.0) ? vel_y * ls_ym : vel_y * ls_yp;

            Set::Scalar flux_z = 0.0;
            if (AMREX_SPACEDIM == 3)
            {
                Set::Scalar ls_zm = (k > 0) ? ls_old(i, j, k - 1) : ls_old(i, j, k);
                Set::Scalar ls_zp = ls_old(i, j, k);
                Set::Scalar vel_z = 0.0;
                flux_z = (vel_z > 0.0) ? vel_z * ls_zm : vel_z * ls_zp;
            }

            flux(i, j, k) = -(flux_x / DX[0] + flux_y / DX[1] + flux_z / DX[2]);
        });
    }

    void Reinitialize(int lev)
    {
        const Set::Scalar narrow_band_width = 6.0 * geom[lev].CellSize()[0];
        const Set::Scalar inner_tube = -narrow_band_width;
        const Set::Scalar outer_tube = narrow_band_width;

        // Set the initial condition for reinitialization from Interface_mf
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
            Set::Patch<Set::IntScalar> nbmask_patch = iNarrowBandMask_mf.Patch(lev, mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                     // Set initial condition based on the interface_patch values
            if (nbmask_patch(i, j, k) == 0)
            {
                ls(i, j, k) = 0.0; // Interface cells set to 0
            }
            else if (nbmask_patch(i, j, k) == -1 || nbmask_patch(i, j, k) == 1)
            {
                ls(i, j, k) = (nbmask_patch(i, j, k) == -1) ? inner_tube : outer_tube; // Narrow band cells set to inner_tube or outer_tube
            }
            else if (nbmask_patch(i, j, k) == -2 || nbmask_patch(i, j, k) == 2)
            {
                // Cells outside the narrow band maintain their signed distance value
                ls(i, j, k) = (nbmask_patch(i, j, k) == -2) ? inner_tube : outer_tube;
            }

            });
        }
       
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

    // Update the narrow band of the level set function
    void UpdateNarrowBand(int lev)
    {
        const Set::Scalar narrow_band_width = 6.0 * geom[lev].CellSize()[0];
        const Set::Scalar inner_tube = -narrow_band_width;
        const Set::Scalar outer_tube = narrow_band_width;

        // Iterate over all the patches on this level
        for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> ls = ls_mf.Patch(lev, mfi);
            Set::Patch<Set::IntScalar> nbmask_patch = iNarrowBandMask_mf.Patch(lev, mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (ls(i, j, k) < inner_tube || ls(i, j, k) > outer_tube)
                {
                    ls(i, j, k) = (ls(i, j, k) > 0) ? outer_tube : inner_tube;
                }
               
	     // Update interface_mf to reflect the current level set values
            
	    if (std::abs(ls(i, j, k)) <= geom[lev].CellSize()[0])
            {
                nbmask_patch(i, j, k) = 0; // Interface cells (closest to zero level set)
            }
            else if (ls(i, j, k) > inner_tube && ls(i, j, k) < outer_tube)
            {
                nbmask_patch(i, j, k) = (ls(i, j, k) > 0) ? 1 : -1; // Narrow band cells
            }
            else
            {
                nbmask_patch(i, j, k) = (ls(i, j, k) > 0) ? 2 : -2; // Outside narrow band
            }

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
    Set::Field<Set::Scalar>    ls_mf;                 // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    ls_old_mf;             // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar>    lsflux_mf;             // Levelset flux field variable 
    Set::Field<Set::Scalar>    interface_mf;          // Variable for Interface Labeling
    Set::Field<Set::IntScalar> iNarrowBandMask_mf;     // Narrow Band Mask

private:

    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    int current_timestep           = 0;
    const int number_of_components = 1;      // Number of components
    const int number_of_ghost_cells = 2;     // Number of ghost cells

    static constexpr int NarrowBand_Innertube   = -3; // InnerTube
    static constexpr int NarrowBand_Outertube   =  3; // OuterTube
    static constexpr int NarrowBand_InsideLS    = -2; // Inside the Levelset
    static constexpr int NarrowBand_OutsideLS   =  2; // Outside the Levelset
    static constexpr int NarrowBand_ZeroLS      =  1; // Closest to the Levelset


    Set::Scalar alpha = NAN;                 // Thermal diffusivity
    Set::Scalar refinement_threshold = NAN ; // Criterion for cell refinement

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC* ls_ic = nullptr;                    // Object used to initialize temperature field
    BC::BC<Set::Scalar>* ls_bc = nullptr;       // Object used to update temp field boundary ghost cells
    BC::BC<Set::IntScalar>* nbmask_bc = nullptr; // Object used to update narrowband mask ghost cells 
};
} // namespace Integrator
#endif
