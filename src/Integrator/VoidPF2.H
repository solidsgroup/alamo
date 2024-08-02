#ifndef INTEGRATOR_VOIDPF2_H
#define INTEGRATOR_VOIDPF2_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Base/Mechanics.H"


#include "IC/IC.H"
#include "BC/BC.H"
#include "IC/PSRead.H"
#include "IC/PS.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "IC/BMP.H"
#include "IC/EulerAngles.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Linear/Isotropic.H"
#include "Operator/Operator.H"



namespace Integrator
{
template<class MODEL>
class VoidPF2: virtual public Base::Mechanics<MODEL>
{
public:

    VoidPF2(): Base::Mechanics<MODEL>() {}
    VoidPF2(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(VoidPF2& value, IO::ParmParse& pp)
    {
        Base::Mechanics<MODEL>::Parse(value, pp);

        pp_queryclass("mechanics",static_cast<Base::Mechanics<MODEL>&>(value));
        pp_queryclass("model", value.model);
        // Read in IC for psi
        if (pp.contains("psi.ic.type"))
        {
            std::string type;
            pp_query("psi.ic.type", type); // Read IC type for the eta field
            //if (type == "ellipse") value.ic_psi = new IC::Ellipse(value.geom, pp, "psi.ic.ellipse");
            if (type == "constant") value.ic_psi = new IC::Constant(value.geom, pp, "psi.ic.constant");
            else if (type == "psread") value.ic_psi = new IC::PSRead(value.geom,pp,"psi.ic.psread");
            
            else Util::Abort(INFO, "Invalid value for psi.ic.type: ", type);

            value.bc_psi = new BC::Constant(1, pp, "psi.bc.constant");
            
               
                   
            value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", true);
            value.RegisterNewFab(value.psi_old_mf, value.bc_psi, 1, 2, "psiold", false);
            value.RegisterNewFab(value.p_mf, value.bc_psi, 1, 3, "p", true);
            value.psi_on = true;
            value.bc_p = new BC::Nothing();   
        }
        pp_query("eta_ref_threshold", value.m_eta_ref_threshold); // Refinement threshold based on eta

        pp_query_required("spall.M",value.spall.M);
        pp_query_required("spall.kappa",value.spall.kappa);
        pp_query_required("spall.cp",value.spall.cp);
        pp_query_required("spall.po",value.spall.po);
        pp_query_required("spall.elastic_mult",value.spall.elastic_mult);
        pp_query_required("spall.small",value.spall.small); // regularization to set minimum modulus multiplier

    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        if (psi_on) ic_psi->Initialize(lev, psi_mf);
        if (psi_on) ic_psi->Initialize(lev, psi_old_mf);
        p_mf[lev]->setVal(0.0);
    }

    virtual void UpdateModel(int a_step, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (a_step > 0) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            model_mf[lev]->setVal(model);
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_old_mf[lev], geom[lev]);
        }

    }


    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        BL_PROFILE("VoidPF2::Advance");
        Base::Mechanics<MODEL>::Advance(lev, time, dt);
        std::swap(psi_old_mf[lev], psi_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        //Set::Scalar Lnow = L(time);
        //Set::Scalar lambdaT = lambda(time);

        for (amrex::MFIter mfi(*psi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(1);
            bx = bx & domain;
            Set::Patch<const Set::Matrix> sig  = stress_mf.Patch(lev,mfi);
            Set::Patch<const Set::Matrix> eps  = strain_mf.Patch(lev,mfi);
            amrex::Array4<const Set::Scalar> const& psi = (*psi_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& psinew = (*psi_mf[lev]).array(mfi);
            Set::Patch<Set::Scalar>       p    = p_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar driving_force = 0.0;

                driving_force += spall.cp * 2.0 * psi(i,j,k)*(1 - 2.0*psi(i,j,k))*(1.0 - psi(i,j,k));
                //driving_force += spall.cp * 2.0 * psi(i,j,k)*(2.0 * psi(i, j, k) * psi(i, j, k) - 3.0 * psi(i, j, k) + 1.0);
                driving_force += spall.kappa * Numeric::Laplacian(psi,i,j,k,0,DX);

                Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k, 0);
                Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);

                p(i,j,k) = 0.5 * sig_avg.trace();
      
                if (p(i,j,k) > spall.po)
                    driving_force += spall.elastic_mult*(p(i,j,k) - spall.po);


                psinew(i, j, k) = psi(i,j,k) - spall.M*dt*driving_force;
                if (psinew(i, j, k) < 0.0) psinew(i, j, k) = 0.0;
                if (psinew(i, j, k) > 1.0) psinew(i, j, k) = 1.0;
            });
        }
    }


    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            if (psi_on)
            {
                amrex::Array4<Set::Scalar> const& psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    {
                        Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(), sten);
                        if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }
        }
    }


protected:
    MODEL model;
    IC::IC* ic_psi = nullptr;
    BC::BC<Set::Scalar>* bc_psi = nullptr;
    BC::BC<Set::Scalar>*bc_p;
    Set::Scalar m_eta_ref_threshold = 1E-2;
    Set::Field<Set::Scalar> psi_old_mf;
    Set::Field<Set::Scalar> p_mf;


    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;

     struct {
        Set::Scalar M = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar po = NAN;
        Set::Scalar cp = NAN;
        Set::Scalar elastic_mult = NAN;
        Set::Scalar small = NAN;
    } spall;


};










} // namespace Integrator
#endif
