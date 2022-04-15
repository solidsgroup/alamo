#ifndef INTEGRATOR_DYNAMICS_H
#define INTEGRATOR_DYNAMICS_H
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

#include "Model/Solid/Linear/Isotropic.H"

#include "Numeric/Stencil.H"

#define TEMP_OLD(i, j, k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i, j, k)))
#define TEMP(i, j, k) Temp_box(amrex::IntVect(AMREX_D_DECL(i, j, k)))

namespace Integrator
{
class Dynamics : virtual public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    Dynamics() : Integrator() {}

    static void Parse(Dynamics &value, IO::ParmParse &pp)
    {
        pp.queryclass("model",value);
    }

protected:

    /// \brief Use the #ic object to initialize #Temp
    void Initialize(int lev)
    {
        Set::Vector b0(AMREX_D_DECL(1.0,0.0,0.0));

        u_mf[lev]->setVal(Set::Vector::Zero());
        eps_mf[lev]->setVal(Set::Matrix::Zero());
        sig_mf[lev]->setVal(Set::Matrix::Zero());
        b_mf[lev]->setVal(b0);

        RegisterGeneralFab(unew_mf,1,2);

        RegisterGeneralFab(u_mf,1,2);
        RegisterGeneralFab(eps_mf,1,2);
        RegisterGeneralFab(sig_mf,1,2);
        RegisterGeneralFab(b_mf,1,2);

        //ic->Initialize(lev,temp_old_mf);
    }

    /// \brief Integrate the heat equation
    void Advance(int lev, amrex::Real /*time*/, amrex::Real dt)
    {
        std::swap(*unew_mf[lev], *u_mf[lev]);
        
        const amrex::Real *DX = geom[lev].CellSize();

        for (amrex::MFIter mfi(*unew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<Set::Vector>         const &unew = (*unew_mf[lev]).array(mfi);
            amrex::Array4<const Set::Vector>   const &u    = (*u_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const &eps  = (*eps_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const &sig  = (*sig_mf[lev]).array(mfi);
            amrex::Array4<const Set::Vector>   const &b    = (*b_mf[lev]).array(mfi);
        
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                eps(i,j,k) = Numeric::Gradient(u,i,j,k,DX);
                sig(i,j,k) = model.DW(eps(i,j,k));
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                Set::Vector rhoudotdot = Numeric::Divergence(sig,i,j,k,DX) + b(i,j,k);
            });

        }

    }

    /// \brief Tag cells for mesh refinement based on temperature gradient
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {}

protected:

    Set::Field<Set::Vector> unew_mf;

    Set::Field<Set::Vector> u_mf;
    Set::Field<Set::Matrix> eps_mf;
    Model::Solid::Linear::Isotropic model;    
    Set::Field<Set::Matrix> sig_mf;
    Set::Field<Set::Vector> b_mf;

private:
    int number_of_components = 1;            ///< Number of components
    int number_of_ghost_cells = 2;           ///< Number of ghost cells

    //IC::IC *ic;                              ///< Object used to initialize temperature field
    //BC::BC<Set::Vector> *bc;                 ///< Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
