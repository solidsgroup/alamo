#ifndef INTEGRATOR_HEATCONDUCTION_H
#define INTEGRATOR_HEATCONDUCTION_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <functional>
#include <cmath>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/Constant.H"

#include "IC/Cylinder.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"

#include "Numeric/Stencil.H"

#define TEMP_OLD(i, j, k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i, j, k)))
#define TEMP(i, j, k) Temp_box(amrex::IntVect(AMREX_D_DECL(i, j, k)))

/// \class HeatConduction::Integrator
/// \brief Solve the heat diffusion equation in 2D or 3D.
///
/// This is a nontrivial example class that demonstrates how to write an explicit AMR
/// code using the #Integrator virtual class that abstracts the AmrBase class
/// from Amrex.
///
/// For more details:
///    - See documentation on #Initialize for input parameters
///    - See documentation on #Advance for equations and discretization
///    - See documentation on #TagCellsForRefinement for cell refinement criteria
/// For boundary conditions:
///    - See #BC
/// For initial conditions:
///    - See #Cylinder or #Constant
namespace Integrator
{
class HeatConduction : public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    HeatConduction() : Integrator()
    {
        IO::ParmParse pp("heat");
        pp.query("alpha", alpha);
        pp.query("refinement_threshold", refinement_threshold);
        IO::ParmParse pph("physics");
            pph.query("M",M);
            pph.query("kappa",kappa);
            pph.query("w1",w1);
            pph.query("w12",w12);
            pph.query("w0",w0);
            pph.query("rho1",rho1);
            pph.query("rho0",rho0);
            pph.query("k1",k1);
            pph.query("k0",k0);
            pph.query("P",P);
            pph.query("r",r);
            pph.query("n",n);
        // Read in initial condition
        {



 
            IO::ParmParse pp("ic");
            std::string type = "sphere";
            pp.query("type",type);
            if (type == "sphere")
            {
                ic = new IC::Sphere(geom);
                pp.queryclass("sphere",*static_cast<IC::Sphere*>(ic));
            }
            else
            {
                ic = new IC::Constant(geom);
                pp.queryclass("constant",*static_cast<IC::Constant*>(ic));
            }

        }
        {
            IO::ParmParse pp("bc");
            bc = new BC::Constant(1);
            pp.queryclass("temp",*static_cast<BC::Constant *>(bc));


            EtaBC = new BC::Constant(1);////////
                pp.queryclass("eta", *static_cast<BC::Constant *>(EtaBC));///////
        }



        RegisterNewFab(temp_mf,     bc, number_of_components, number_of_ghost_cells, "Temp",true);
        RegisterNewFab(temp_old_mf, bc, number_of_components, number_of_ghost_cells, "Temp_old",false);


        RegisterNewFab(Eta_mf, EtaBC, number_of_components, number_of_ghost_cells, "Eta", true);/////
        RegisterNewFab(Eta_old_mf, EtaBC, number_of_components, number_of_ghost_cells, "Eta_old", false);///////
    }

protected:

    /// \brief Use the #ic object to initialize #Temp
    void Initialize(int lev)
    {
        ic->Initialize(lev,temp_old_mf);
        Eta_mf[lev]->setVal(1.0);  //initializing afte removing wedge BC
        Eta_old_mf[lev]->setVal(1.0); // initializing after removing wedge BC
    }


    /// \brief Integrate the heat equation
    void Advance(int lev, amrex::Real time, amrex::Real dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        std::swap(Eta_old_mf[lev], Eta_mf[lev]);////
        ///std::swap(Temp_old_mf[lev], Temp_mf[lev]);////
        Set::Scalar a0 = w0, a1 = 0.0, a2 = -5 * w1 + 16 * w12 - 11 * a0, a3 = 14 * w1 - 32 * w12 + 18 * a0, a4 = -8 * w1 + 16 * w12 - 8 * a0;

        // Get the cell size corresponding to this level
        const amrex::Real *DX = geom[lev].CellSize();

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box &bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            amrex::Array4<const Set::Scalar> const &temp_old = (*temp_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const &temp     = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const &Eta              = (*Eta_mf[lev]).array(mfi);////
            amrex::Array4<const Set::Scalar> const &Eta_old    = (*Eta_old_mf[lev]).array(mfi);////

            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {


                // Do the physics!
                // Note that Numeric::Laplacian is an inlined function so there is no overhead.
                // You can calculate the derivatives yourself if you want.
                //Set::Vector eta_grad = Numeric::Gradient(Eta_old,i,j,k,0,DX);
                //Set::Vector temp_grad = Numeric::Gradient(temp_old,i,j,k, 0, DX);
                //Set::Vector normvec = eta_grad/Eta_old(i,j,k);
                //double in = inner_product(begin(eta_grad),end(eta_grad),begin(eta_grad),0.0);
                //double abso = sqrt(in);
                //double test = inner_product(begin(normvec),end(normvec),begin(temp_grad),0.0);
                ///for(i = 0; i < 3; i++) Eta_old[i] *= abso/Eta_old(i,j,k);
                Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);

                Eta(i, j, k) = Eta_old(i, j, k) -
                                (2.4*((r*pow(P,n)))) * dt * (a1 + 2 * a2 * Eta_old(i, j, k) + 3 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k) - kappa * eta_lap);
                if(time>0)
                {
                temp(i,j,k) = temp_old(i,j,k); //+ dt*alpha*(test+Numeric::Laplacian(temp_old,i,j,k,0,DX)+abso/Eta_old(i,j,k)*0);
                }
            });
        }
    }

    /// \brief Tag cells for mesh refinement based on temperature gradient
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {
        // Get cell dimensions as done above.
        const Set::Scalar *DX = geom[lev].CellSize();
        // Calculate the diagonal.
        Set::Scalar dr  = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Iterate over the patches on this level
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char>         const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const &temp = (*temp_mf[lev]).array(mfi);
            //amrex::Array4<const Set::Scalar> const &Eta  = (*Eta_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                // Calculate the temperature gradient.
                Set::Vector grad = Numeric::Gradient(temp,i,j,k,0,DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i,j,k) = amrex::TagBox::SET;

                //Set::Vector gradeta = Numeric::Gradient(Eta,i,j,k,0,DX);
                //if (gradeta.lpNorm<2>() * dr > refinement_threshold)
                //tags(i,j,k) = amrex::TagBox::SET;

            });
        }
    }
    void Regrid(int lev, Set::Scalar )
    {
        if (lev < finest_level) return;
        temp_mf[lev]->setVal(0.0);
               
        Util::Message(INFO, "Regridding on level ", lev);
    }

private:
    int number_of_components = 1;            ///< Number of components
    int number_of_ghost_cells = 1;           ///< Number of ghost cells
    Set::Scalar M;                  // Burn speed parameter
    Set::Scalar kappa;              // Burn width parameter
    Set::Scalar w1, w12, w0;        // Energy (before, during, after)
    Set::Scalar rho1, rho0;         // Density (before and after)
    Set::Scalar k1, k0;             // Thermal conductivity (before and after)
    Set::Scalar cp1, cp0;           // Specific heat (before and after)
    Set::Scalar qdotburn;

    Set::Scalar P;
    Set::Scalar r;
    Set::Scalar n;

    Set::Field<Set::Scalar> temp_mf;     ///< Temperature field variable (current timestep)
    Set::Field<Set::Scalar> temp_old_mf;     ///< Temperature field variable (previous timestep)
    Set::Field<Set::Scalar> Eta_mf;/////temp - adding eta into file
    Set::Field<Set::Scalar> Eta_old_mf;


    amrex::Real alpha = 1.0; ///< Thermal diffusivity
    amrex::Real refinement_threshold = 0.01; ///< Criterion for cell refinement

    IC::IC *ic;                 ///< Object used to initialize temperature field
    IC::IC *EtaIC; ///

    BC::BC<Set::Scalar> *bc;                 ///< Object used to update temp field boundary ghost cells
    BC::BC<Set::Scalar> *EtaBC;///


};
} // namespace Integrator
#endif
