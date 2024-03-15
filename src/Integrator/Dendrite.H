//
// This is a BASIC phase field method for dendritic growth as described in
// "A Numerical Approach to Three-Dimensional Dendritic Solidification" by Ryo Kobayashi (1994)
// Reference: https://doi.org/10.1080/10586458.1994.10504577
//

#ifndef INTEGRATOR_DENDRITE_H // Include guards
#define INTEGRATOR_DENDRITE_H // 

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

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
class Dendrite : virtual public Integrator
{
public:
    // Empty constructor
    Dendrite() : Integrator()
    {}

    // Constructor that triggers parse
    Dendrite(IO::ParmParse& pp) : Dendrite()
    {
        Parse(*this, pp);
    }

    static void Parse(Dendrite& value, IO::ParmParse& pp)
    {
        pp.query("alpha", value.alpha); // Pre-multiplier of "m" barrier height
        pp.query("delta", value.delta); // Anisotropy factor
        pp.query("gamma", value.gamma); // Anisotropic temperature coupling factor

        pp.query("eps", value.eps); // Diffuse boundary width
        pp.query("tau", value.tau); // Diffusive timescale

        // Refinement criteria
        pp.query("heat.refinement_threshold", value.refinement_threshold_temp);
        pp.query("phi.refinement_threshold", value.refinement_threshold_phi);

        // Expressions for ICs
        value.ic_temp = new IC::Expression(value.geom, pp, "ic.temp");
        value.ic_phi = new IC::Expression(value.geom, pp, "ic.phi");

        // Use a constant BC object for temperature
        value.bc_temp = new BC::Constant(1);
        value.bc_phi = new BC::Constant(1);
        pp.queryclass("bc.temp", *static_cast<BC::Constant*>(value.bc_temp));
        pp.queryclass("bc.phi", *static_cast<BC::Constant*>(value.bc_phi));

        value.RegisterNewFab(value.temp_mf, value.bc_temp, value.number_of_components, value.number_of_ghost_cells, "Temp", true);
        value.RegisterNewFab(value.temp_old_mf, value.bc_temp, value.number_of_components, value.number_of_ghost_cells, "Temp_old", false);
        value.RegisterNewFab(value.phi_mf, value.bc_phi, value.number_of_components, value.number_of_ghost_cells, "phi", true);
        value.RegisterNewFab(value.phi_old_mf, value.bc_phi, value.number_of_components, value.number_of_ghost_cells, "phi_old", false);
    }
protected:

    // Use the ic object to initialize the temperature field
    void Initialize(int lev)
    {
        ic_temp->Initialize(lev, temp_old_mf);
        ic_phi->Initialize(lev, phi_old_mf);
        ic_temp->Initialize(lev, temp_mf);
        ic_phi->Initialize(lev, phi_mf);
    }


    // Integrate the heat equation
    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        std::swap(*phi_mf[lev], *phi_old_mf[lev]);

        // Get the cell size corresponding to this level
        const Set::Scalar* DX = geom[lev].CellSize();

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box& bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            amrex::Array4<const Set::Scalar> const& temp = (*temp_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& temp_new = (*temp_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = (*phi_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& phi_new = (*phi_mf[lev]).array(mfi);

            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {

                Set::Vector gradphi = Numeric::Gradient(phi, i, j, k, 0, DX);

                // 
                // Calculate anisotropy function
                // [ equation 2.12 in reference]
                //
                Set::Scalar v44 = AMREX_D_TERM(
                    gradphi(0) * gradphi(0) * gradphi(0) * gradphi(0),
                    +gradphi(1) * gradphi(1) * gradphi(1) * gradphi(1),
                    +gradphi(2) * gradphi(2) * gradphi(2) * gradphi(2));
                Set::Scalar v22 = gradphi.squaredNorm(); v22 *= v22;
                Set::Scalar sigma = 1.0 - delta * (1.0 - v44 / (v22 + 1E-12));

                // 
                // Calculate anisotropic barrier height
                // [ unnumbered equation on page 67 in the reference ]
                //
                Set::Scalar m = -(alpha / Set::Constant::Pi) * std::atan(gamma * sigma * temp(i, j, k));


                //
                // Evolve order paramter
                // [ equation 2.10 in the reference ]
                //                
                Set::Scalar lapphi = Numeric::Laplacian(phi, i, j, k, 0, DX);
                Set::Scalar Dphi = (eps * eps * lapphi + phi(i, j, k) * (1.0 - phi(i, j, k)) * (phi(i, j, k) - 0.5 + m)) / tau;
                phi_new(i, j, k) = phi(i, j, k) + dt * Dphi;

                //
                // Evolve temperature
                // [ equation 2.11 in the reference ]
                //
                Set::Scalar laptemp = Numeric::Laplacian(temp, i, j, k, 0, DX);
                temp_new(i, j, k) = temp(i, j, k) + dt * (laptemp + Dphi);
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
            amrex::Array4<Set::Scalar>  const& phi = (*phi_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the temperature gradient.
                Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX);
                Set::Vector grad_phi = Numeric::Gradient(phi, i, j, k, 0, DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad_temp.lpNorm<2>() * dr > refinement_threshold_temp)
                    tags(i, j, k) = amrex::TagBox::SET;
                if (grad_phi.lpNorm<2>() * dr > refinement_threshold_phi)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> phi_old_mf;

private:
    int number_of_components = 1;            // Number of components
    int number_of_ghost_cells = 1;           // Number of ghost cells


    // calculation of m
    Set::Scalar alpha, delta, gamma;
    // evolution of p
    Set::Scalar eps, tau;

    Set::Scalar refinement_threshold_temp = 0.01; // Criterion for cell refinement
    Set::Scalar refinement_threshold_phi = 0.01; // Criterion for cell refinement

    IC::IC* ic_temp, * ic_phi;
    BC::BC<Set::Scalar>* bc_temp, * bc_phi;
};
} // namespace Integrator
#endif
