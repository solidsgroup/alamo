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

namespace Integrator
{
class Dynamics : virtual public Integrator
{
public:
    Dynamics() : Integrator() 
    {
        //RegisterNewFab(empty_mf,bc,1,2,"empty",true);
        //RegisterNewFab(eta_new_mf, mybc, number_of_grains, number_of_ghost_cells, "Eta",true);
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 0, "disp",true);
        
        RegisterGeneralFab(unew_mf,1,2);
        RegisterGeneralFab(vnew_mf,1,2);
        RegisterGeneralFab(u_mf,1,2,"u");
        RegisterGeneralFab(v_mf,1,2,"v");
        RegisterGeneralFab(eps_mf,1,2);
        RegisterGeneralFab(sig_mf,1,2);
        RegisterGeneralFab(b_mf,1,2,"b");
    }

    static void Parse(Dynamics &value, IO::ParmParse &pp)
    {
        //value.bc = new BC::Constant(1);
        //pp.queryclass(*static_cast<BC::Constant *>(value.bc)); // See :ref:`BC::Constant`
        
        pp.queryclass("model",value.model);
    }

protected:

    void Initialize(int lev)
    {
        Set::Vector b0(AMREX_D_DECL(1.0,0.0,0.0));

        unew_mf[lev]->setVal(Set::Vector::Zero());
        vnew_mf[lev]->setVal(Set::Vector::Zero());

        u_mf[lev]->setVal(Set::Vector::Zero());
        v_mf[lev]->setVal(Set::Vector::Zero());
        eps_mf[lev]->setVal(Set::Matrix::Zero());
        sig_mf[lev]->setVal(Set::Matrix::Zero());
        
        b_mf[lev]->setVal(b0);

        disp_mf[lev]->setVal(0.0);
        //ic->Initialize(lev,temp_old_mf);
    }

    void Advance(int lev, amrex::Real /*time*/, amrex::Real dt)
    {
        std::swap(*unew_mf[lev], *u_mf[lev]);
        std::swap(*vnew_mf[lev], *v_mf[lev]);

        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        domain.grow(-1);
        
        const amrex::Real *DX = geom[lev].CellSize();

        for (amrex::MFIter mfi(*unew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox() & domain;

            amrex::Array4<Set::Vector>         const &unew = (*unew_mf[lev]).array(mfi);
            amrex::Array4<Set::Vector>         const &vnew = (*vnew_mf[lev]).array(mfi);

            amrex::Array4<const Set::Vector>   const &u    = (*u_mf[lev]).array(mfi);
            amrex::Array4<Set::Vector>         const &v    = (*v_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const &eps  = (*eps_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const &sig  = (*sig_mf[lev]).array(mfi);
            amrex::Array4<const Set::Vector>   const &b    = (*b_mf[lev]).array(mfi);
        
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                //std::array<Numeric::StencilType,AMREX_SPACEDIM> sten = Numeric::GetStencil(i,j,k,bx);
                
                eps(i,j,k) = Numeric::Gradient(u,i,j,k,DX);
                sig(i,j,k) = model.DW(eps(i,j,k));
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                //std::array<Numeric::StencilType,AMREX_SPACEDIM> sten = Numeric::GetStencil(i,j,k,bx);

                //Set::Vector udotdot = (Numeric::Divergence(sig,i,j,k,DX) + b(i,j,k)) / rho;
                Set::Matrix3 gradgradu = Numeric::Hessian(u,i,j,k,DX);
                Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> C = model.DDW(eps(i,j,k));

                Set::Vector udotdot = (C * gradgradu + b(i,j,k)) / rho;

                vnew(i,j,k) = v(i,j,k) + dt * udotdot;
                unew(i,j,k) = u(i,j,k) + dt*v(i,j,k);
            });

            amrex::Array4<Set::Scalar> const disp = (*disp_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                disp(i,j,k,0) = u(i,j,k)[0];
                disp(i,j,k,1) = u(i,j,k)[1];
            });
            
        }
    }

    void TagCellsForRefinement(int /*lev*/, amrex::TagBoxArray &/*a_tags*/, amrex::Real /*time*/, int /*ngrow*/)
    {}

protected:

    Set::Field<Set::Scalar> disp_mf;

    Set::Field<Set::Vector> unew_mf;
    Set::Field<Set::Vector> vnew_mf;

    Set::Field<Set::Vector> u_mf;
    Set::Field<Set::Vector> v_mf;
    Set::Field<Set::Matrix> eps_mf;
    Set::Field<Set::Matrix> sig_mf;
    Set::Field<Set::Vector> b_mf;

    Model::Solid::Linear::Isotropic model;
    Set::Scalar rho = 1.0;

private:
    int number_of_components = 1;            ///< Number of components
    int number_of_ghost_cells = 2;           ///< Number of ghost cells

    //IC::IC *ic;                              ///< Object used to initialize temperature field
//    BC::BC<Set::Vector> *bc;                 ///< Object used to update temp field boundary ghost cells
    BC::BC<Set::Scalar> *bc;                 ///< Object used to update temp field boundary ghost cells

};
} // namespace Integrator
#endif
