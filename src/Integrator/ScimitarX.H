
#ifndef INTEGRATOR_SCIMITARX_H // Include guards
#define INTEGRATOR_SCIMITARX_H // 

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
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
class ScimitarX : virtual public Integrator
{
public:
    // Empty constructor
    ScimitarX(int a_nghost = 2) : 
        Integrator(),
        number_of_ghost_cells(a_nghost)
    {}

    // Constructor that triggers parse
    ScimitarX(IO::ParmParse& pp) : ScimitarX()
    {
        Parse(*this, pp);
    }

    // The Parse function initializes the HeatConduction object using
    // a parser, pp. 
    // Note that this is a static function, which means it does not have
    // direct access to member variables. Instead, it initializes the variables
    // inside the argument, "value", and so all references to member items are
    // prefixed by "value."
    static void Parse(ScimitarX& value, IO::ParmParse& pp)
    {
        // Diffusion coefficient [1.0]
        pp.query("heat.alpha", value.alpha);
        // Criterion for mesh refinement [0.01]
        pp.query("heat.refinement_threshold", value.refinement_threshold);

        std::string type = "sphere";
        // Initial condition type ([sphere], constant, expression)
        pp.query("ic.type", type);
        if (type == "sphere")          value.ic = new IC::Sphere(value.geom, pp, "ic.sphere");
        else if (type == "constant")   value.ic = new IC::Constant(value.geom, pp, "ic.constant");
        else if (type == "expression") value.ic = new IC::Expression(value.geom, pp, "ic.expression");
        else  Util::Abort(INFO, "Invalid ic.type ", type);

        // Use a constant BC object for temperature
        value.bc = new BC::Constant(1);
        // :ref:`BC::Constant` parameters with prefix bc.temp
        pp.queryclass("bc.temp", *static_cast<BC::Constant*>(value.bc));

        // Register the temperature and old temperature fields.
        // temp_mf and temp_old_mf are defined near the bottom of this Header file.
        value.RegisterNewFab(value.temp_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp", true);
        value.RegisterNewFab(value.temp_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp_old", false);
    }

protected:

    // Use the ic object to initialize the temperature field
    void Initialize(int lev)
    {
        ic->Initialize(lev, temp_old_mf);
    }

    // Integrate the heat equation
    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);

        // Get the cell size corresponding to this level
        const Set::Scalar* DX = geom[lev].CellSize();

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box& bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            amrex::Array4<const Set::Scalar> const& temp_old = (*temp_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& temp = (*temp_mf[lev]).array(mfi);

            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Do the physics!
                // Note that Numeric::Laplacian is an inlined function so there is no overhead.
                // You can calculate the derivatives yourself if you want.
                temp(i, j, k) = temp_old(i, j, k) + dt * alpha * Numeric::Laplacian(temp_old, i, j, k, 0, DX);
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
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char>         const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const& temp = (*temp_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the temperature gradient.
                Set::Vector grad = Numeric::Gradient(temp, i, j, k, 0, DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<Set::Scalar> temp_mf;         // Temperature field variable (current timestep)
    Set::Field<Set::Scalar> temp_old_mf;     // Temperature field variable (previous timestep)

private:
    int number_of_components = 1;            // Number of components
    int number_of_ghost_cells = 2;           // Number of ghost cells


    Set::Scalar alpha = 1.0;                 // Thermal diffusivity
    Set::Scalar refinement_threshold = 0.01; // Criterion for cell refinement

    IC::IC* ic;                              // Object used to initialize temperature field
    BC::BC<Set::Scalar>* bc;                 // Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
