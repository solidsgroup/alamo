#ifndef INTEGRATOR_HEATCONDUCTION_H
#define INTEGRATOR_HEATCONDUCTION_H
#include <iostream>
#include <fstream>
#include <iomanip>

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
#include "IC/Expression.H"

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
class HeatConduction : virtual public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    HeatConduction() : Integrator()
    {}
    HeatConduction(IO::ParmParse &pp) : HeatConduction()
    {Parse(*this,pp);}

    static void Parse(HeatConduction &value, IO::ParmParse &pp)
    {
        pp.query("heat.alpha", value.alpha);
        pp.query("heat.refinement_threshold", value.refinement_threshold);

        std::string type = "sphere";
        pp.query("ic.type",type);
        if (type == "sphere")          value.ic = new IC::Sphere(value.geom,pp,"ic.sphere");
        else if (type == "constant")   value.ic = new IC::Constant(value.geom, pp, "ic.constant");
        else if (type == "expression") value.ic = new IC::Expression(value.geom,pp,"ic.expression");
        else  Util::Abort(INFO,"Invalid ic.type ", type);

            
        value.bc = new BC::Constant(1);
        pp.queryclass("bc.temp",*static_cast<BC::Constant *>(value.bc));
                
        value.RegisterNewFab(value.temp_mf,     value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp",true);
        value.RegisterNewFab(value.temp_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp_old",false);
    }

protected:

    /// \brief Use the #ic object to initialize #Temp
    void Initialize(int lev)
    {
        ic->Initialize(lev,temp_old_mf);
    }

    /// \brief Integrate the heat equation
    void Advance(int lev, amrex::Real /*time*/, amrex::Real dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        
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
        
            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                    {
                                        // Do the physics!
                                        // Note that Numeric::Laplacian is an inlined function so there is no overhead.
                                        // You can calculate the derivatives yourself if you want.
                                        temp(i,j,k) = temp_old(i,j,k) + dt * alpha * Numeric::Laplacian(temp_old,i,j,k,0,DX);
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
            
            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                    {
                                        // Calculate the temperature gradient.
                                        Set::Vector grad = Numeric::Gradient(temp,i,j,k,0,DX);

                                        // Is the gradient * cell_size too big? If so, then
                                        // mark this cell as needing refinement.
                                        if (grad.lpNorm<2>() * dr > refinement_threshold)
                                            tags(i,j,k) = amrex::TagBox::SET;
                                    });
        }
    }

protected:
    Set::Field<Set::Scalar> temp_mf;         ///< Temperature field variable (current timestep)
    Set::Field<Set::Scalar> temp_old_mf;     ///< Temperature field variable (previous timestep)

private:
    int number_of_components = 1;            ///< Number of components
    int number_of_ghost_cells = 2;           ///< Number of ghost cells


    Set::Scalar alpha = 1.0;                 ///< Thermal diffusivity
    Set::Scalar refinement_threshold = 0.01; ///< Criterion for cell refinement

    IC::IC *ic;                              ///< Object used to initialize temperature field
    BC::BC<Set::Scalar> *bc;                 ///< Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
