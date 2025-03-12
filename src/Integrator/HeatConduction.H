//
// This implements a basic heat conduction method in Alamo.
// The partial differential equation to be solved is
//
// .. math::
//
//    \frac{\partial T}{\partial t} = \alpha\,\Delta T
//
// where :math:`T` is temperature, :math:`t` is time, and :math:`alpha` 
// is the thermal diffusivity.
// Integration is performed explicitly in time using forward Euler, and
// differentiation is performed using the finite difference method.
//

#ifndef INTEGRATOR_HEATCONDUCTION_H // Include guards
#define INTEGRATOR_HEATCONDUCTION_H // 

// AMReX Includes
#include "AMReX_Array4.H"
#include "AMReX_GpuComplex.H"
#include "AMReX_MFIter.H"
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#ifdef ALAMO_FFT
#include <AMReX_FFT.H>
#endif

// Alamo Includes
#include "Set/Base.H"
#include "Integrator.H"
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
class HeatConduction : virtual public Integrator
{
public:
    // Empty constructor
    HeatConduction(int a_nghost = 2) : 
        Integrator(),
        number_of_ghost_cells(a_nghost)
    {}

    // Constructor that triggers parse
    HeatConduction(IO::ParmParse& pp) : HeatConduction()
    {
        Parse(*this, pp);
    }

    virtual ~HeatConduction()
    {
        delete ic;
        delete bc;
    }

    // The Parse function initializes the HeatConduction object using
    // a parser, pp. 
    // Note that this is a static function, which means it does not have
    // direct access to member variables. Instead, it initializes the variables
    // inside the argument, "value", and so all references to member items are
    // prefixed by "value."
    static void Parse(HeatConduction& value, IO::ParmParse& pp)
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

        // Initial condition type.
        pp.select_default<IC::Constant,IC::Sphere,IC::Expression>("ic",value.ic,value.geom);

        // Select BC object for temperature
        pp.select_default<BC::Constant,BC::Expression>("bc.temp",value.bc,1);

        // Register the temperature and old temperature fields.
        // temp_mf and temp_old_mf are defined near the bottom of this Header file.
        value.RegisterNewFab(value.temp_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp", true);
        value.RegisterNewFab(value.temp_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp_old", false);

        #ifdef ALAMO_FFT
        value.RegisterNewFab(value.temp_fft_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "TempFFT", true);
        #endif

    }

protected:

    // Use the ic object to initialize the temperature field
    void Initialize(int lev)
    {
        ic->Initialize(lev, temp_old_mf);
        ic->Initialize(lev, temp_mf);
        #ifdef ALAMO_FFT
        ic->Initialize(lev, temp_fft_mf);
        #endif
    }

#ifdef ALAMO_FFT
    virtual void
    TimeStepBegin(Set::Scalar /*time*/, int iter)
    {
        if (iter % 10)
            return;
        Set::Scalar timestep = 10.0*this->timestep;

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::FFT::R2C my_fft(this->geom[lev].Domain());
            auto const& [cba, cdm] = my_fft.getSpectralDataLayout();
            amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > Temp_hat(cba, cdm, 1, 0);
            my_fft.forward(*temp_fft_mf[lev], Temp_hat);

            const Set::Scalar* DX = geom[lev].CellSize();

            // Set::Scalar
            //     pi_Lx = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(0) / DX[0],
            //     pi_Ly = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(1) / DX[1];

            Set::Scalar
                AMREX_D_DECL(
                    pi_Lx = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(0) / DX[0],
                    pi_Ly = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(1) / DX[1],
                    pi_Lz = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(2) / DX[2]);


            Set::Scalar scaling = 1.0 / geom[lev].Domain().d_numPts();
            //Util::Message(INFO,geom[lev].Domain().d_numPts());

            //for (amrex::MFIter mfi(Temp_hat, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            for (amrex::MFIter mfi(Temp_hat); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & T_hat =  Temp_hat.array(mfi);
                Util::Message(INFO, bx);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {
                    Set::Scalar dt = timestep;

                    // Set::Scalar k1 = m * pi_Lx;

                    // Set::Scalar k2;
                    // if (n < bx.length(1)/2)
                    //     k2 = n * pi_Ly;
                    // else
                    //     k2 = (n - bx.length(1)) * pi_Ly;

                    // Get spectral coordinates
                    AMREX_D_TERM(
                        Set::Scalar k1 = m * pi_Lx;,
                        Set::Scalar k2 = (n < bx.length(1)/2 ? n * pi_Ly : (n - bx.length(1)) * pi_Ly);,
                        Set::Scalar k3 = (p < bx.length(2)/2 ? p * pi_Lz : (p - bx.length(2)) * pi_Lz););

                    Set::Scalar lap = AMREX_D_TERM(k1 * k1, + k2 * k2, + k3*k3);


                    Set::Scalar factor = exp( - alpha * dt * lap); 
                    T_hat(m,n,p) *= factor*scaling;
                });
            }

            my_fft.backward(Temp_hat, *temp_fft_mf[lev]);
        }
    };
#endif

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
            const amrex::Box &bx = mfi.tilebox();
            //Util::Message(INFO,bx);

            // Get an array-accessible handle to the data on this patch.
            // amrex::Array4<const Set::Scalar> const& temp_old = (*temp_old_mf[lev]).array(mfi);
            // amrex::Array4<Set::Scalar>       const& temp = temp_mf.Patch(lev,mfi);//(*temp_mf[lev]).array(mfi);
            Set::Patch<const Set::Scalar>  temp_old = temp_old_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar>        temp     = temp_mf.Patch(lev,mfi);

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

    #ifdef ALAMO_FFT
    Set::Field<Set::Scalar> temp_fft_mf;
    #endif

private:

    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    const int number_of_components = 1;      // Number of components
    const int number_of_ghost_cells = 2;     // Number of ghost cells

    //
    // Definition of user-determined variables.
    //
    // Instantiate all variables to NAN if possible.
    // Default values may be set in Parse using query_default.
    //

    Set::Scalar alpha = NAN;                 // Thermal diffusivity
    Set::Scalar refinement_threshold = NAN ; // Criterion for cell refinement

    //
    // Definition of user-determined pointer variables.
    //
    // These should be set to nullptr. Make sure that they are deleted
    // in the ~HeatConduction destructor.
    //

    IC::IC<Set::Scalar>* ic = nullptr;                    // Object used to initialize temperature field
    BC::BC<Set::Scalar>* bc = nullptr;       // Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
