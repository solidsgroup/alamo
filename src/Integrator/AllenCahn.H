//
// This is a simple implementation of an Allen-Cahn equation governed by
// 
// .. math::
//
//    \frac{\partial\alpha}{\partial t} = 
//               - L \Big(\frac{\lambda}{\epsilon}(2.0\alpha 
//                        - 6\alpha^2 + 4\alpha^3) + \epsilon\,\kappa\,\Delta \alpha \Big)
// 
// where :math:`\alpha(x,t)` is the order parameter.
//

#ifndef INTEGRATOR_ALLENCAHN_H // Include guards
#define INTEGRATOR_ALLENCAHN_H // 

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
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "IC/PSRead.H"
#include "Numeric/Stencil.H"
#include "IC/Random.H"

namespace Integrator
{
class AllenCahn : virtual public Integrator
{
public:
    static constexpr const char* name = "allencahn";

    // Empty constructor
    AllenCahn() : Integrator()
    {}

    // Constructor that triggers parse
    AllenCahn(IO::ParmParse& pp) : AllenCahn()
    {
        Parse(*this, pp);
    }

    // The Parse function initializes the AllenCahn object using
    // a parser, pp. 
    // Note that this is a static function, which means it does not have
    // direct access to member variables. Instead, it initializes the variables
    // inside the argument, "value", and so all references to member items are
    // prefixed by "value."
    static void Parse(AllenCahn& value, IO::ParmParse& pp)
    {
        // Criterion for mesh refinement [0.01]
        pp_query_default("refinement_threshold", value.refinement_threshold, 0.01);

        // Value for :math:`L` (mobility)
        pp_query_default("ch.L",value.ch.L, 1.0);
        // Value for :math:`\epsilon` (diffuse boundary width)
        pp_query_default("ch.eps",value.ch.eps, 0.1);
        // Value for :math:`\kappa` (Interface energy parameter)
        pp_query_default("ch.grad",value.ch.grad, 1.0);
        // Value for :math:`\lambda` (Chemical potential coefficient)
        pp_query_default("ch.chempot",value.ch.chempot, 1.0);
        // Force directional growth: 0=no growth, 1=only positive, -1=only negative
        pp_query_validate("ch.direction",value.ch.direction, {0,1,-1});
        // Time to start forcing directional growth
        pp_query_default("ch.direction_tstart",value.ch.direction_tstart, 0.0);

        std::string type = "sphere";
        // Initial condition type ([sphere], constant, expression, bmp, random)
        pp_query_validate("alpha.ic.type", type, {"sphere","constant","expression","bmp","png","random","psread"});
        if (type == "sphere")          value.ic = new IC::Sphere(value.geom, pp, "alpha.ic.sphere");
        else if (type == "constant")   value.ic = new IC::Constant(value.geom, pp, "alpha.ic.constant");
        else if (type == "expression") value.ic = new IC::Expression(value.geom, pp, "alpha.ic.expression");
        else if (type == "bmp")        value.ic = new IC::BMP(value.geom, pp, "alpha.ic.bmp");
        else if (type == "png")        value.ic = new IC::PNG(value.geom, pp, "alpha.ic.png");
        else if (type == "random")     value.ic = new IC::Random(value.geom, pp, "alpha.ic.random");
        else if (type == "psread")     value.ic = new IC::PSRead(value.geom, pp, "alpha.ic.psread");
        else  Util::Abort(INFO, "Invalid ic.type ", type);

        // Use a constant BC object for temperature
        value.bc = new BC::Constant(1);
        // :ref:`BC::Constant` parameters 
        pp_queryclass("alpha.bc", *static_cast<BC::Constant*>(value.bc));

        // Register the temperature and old temperature fields.
        // alpha_mf and alpha_old_mf are defined near the bottom of this Header file.
        value.RegisterNewFab(value.alpha_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "alpha", true);
        value.RegisterNewFab(value.alpha_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "alpha_old", false);
    }

protected:

    // Use the ic object to initialize the order parameter
    void Initialize(int lev)
    {
        ic->Initialize(lev, alpha_mf);
        ic->Initialize(lev, alpha_old_mf);
    }

    // Integrate the Allen Cahn equation
    void Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {
        std::swap(*alpha_mf[lev], *alpha_old_mf[lev]);

        const Set::Scalar *DX = this->geom[lev].CellSize();
        // Evolve alpha
        for (amrex::MFIter mfi(*alpha_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar>       alpha_new = alpha_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> alpha     = alpha_old_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar driving_force = 0.0;
                // Chemical potential
                Set::Scalar alpha_2 = alpha(i, j, k) * alpha(i, j, k);
                Set::Scalar alpha_3 = alpha_2 * alpha(i, j, k);
                driving_force += (ch.chempot / ch.eps) * (2.0 * alpha(i, j, k) - 6.0 * alpha_2 + 4.0 * alpha_3);
                // Gradient
                Set::Scalar alpha_lap = Numeric::Laplacian(alpha, i, j, k, 0, DX);
                driving_force -= ch.eps * ch.grad * alpha_lap;

                if (time >= ch.direction_tstart)
                {
                    if (ch.direction == 1)
                        driving_force = std::min(0.0, driving_force);
                    else if (ch.direction == -1)
                        driving_force = std::max(0.0, driving_force);
                }

                // Update
                alpha_new(i, j, k) =
                    alpha(i, j, k) - ch.L * dt * driving_force;
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
        for (amrex::MFIter mfi(*alpha_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char>         const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const& alpha = (*alpha_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the temperature gradient.
                Set::Vector grad = Numeric::Gradient(alpha, i, j, k, 0, DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<Set::Scalar> alpha_mf;         // Temperature field variable (current timestep)
    Set::Field<Set::Scalar> alpha_old_mf;     // Temperature field variable (previous timestep)
    IC::IC* ic;                              // Object used to initialize temperature field

    Set::Field<Set::Scalar> &eta_mf = alpha_mf;
    Set::Field<Set::Scalar> &eta_old_mf = alpha_old_mf;

    struct {
        Set::Scalar L = NAN;
        Set::Scalar eps = NAN;
        Set::Scalar grad = NAN;
        Set::Scalar chempot = NAN;
        Set::Scalar direction_tstart = NAN;
        int direction = 0;
    } ch;

private:
    int number_of_components = 1;            // Number of components
    int number_of_ghost_cells = 2;           // Number of ghost cells

    Set::Scalar refinement_threshold = NAN; // Criterion for cell refinement

    BC::BC<Set::Scalar>* bc;                 // Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
