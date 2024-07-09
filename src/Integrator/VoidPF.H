#ifndef INTEGRATOR_VOIDPF_H
#define INTEGRATOR_VOIDPF_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Mechanics.H"


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

#include "Operator/Operator.H"


namespace Integrator
{
template<class MODEL>
class VoidPF: virtual public Mechanics<MODEL>
{
public:

    VoidPF(): Base::Mechanics<MODEL>() {}
    VoidPF(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(VoidPF& value, IO::ParmParse& pp)
    {
        pp_queryclass("mechanics",static_cast<Base::Mechanics<MODEL>&>(value));

        pp_queryclass("model", value.model);
            
        if (pp.contains("euler.ic.type"))
        {
            std::string type;
            pp_query("euler.ic.type", type); // Read IC type for the eta field
            if (type == "eulerangles") value.ic_euler = new IC::EulerAngles(value.geom, pp, "euler.ic.eulerangles");
            else Util::Abort(INFO, "Invalid value for euler.ic.type: ", type);

            
        }

        value.bc_d = new BC::Constant(1,pp,"d.bc.constant");
        value.bc_p = new BC::Nothing();

        value.RegisterNewFab(value.d_mf, value.bc_d, 1, 3, "d",true);
        value.RegisterNewFab(value.d_old_mf, value.bc_d,  1, 3,"dold",false);
        value.RegisterNewFab(value.p_mf, value.bc_d, 1, 3, "p",true);
        //value.spall.ic = new IC::PSRead(value.geom, pp, "spall.ic.psread");
        
        // Initialize PSRead with parameters from the input file
        //std::string filename;
        //pp.query_file("spall.ic.psread", filename);
        std::string ictype;
        pp_query_validate("spall.ic.type",ictype,{"psread"});
        if (ictype == "psread") value.spall.ic = new IC::PSRead(value.geom, pp,"spall.ic.psread");

        pp_query_required("spall.M",value.spall.M);
        pp_query_required("spall.kappa",value.spall.kappa);
        pp_query_required("spall.cp",value.spall.cp);
        pp_query_required("spall.po",value.spall.po);
        pp_query_required("spall.elastic_mult",value.spall.elastic_mult);
        pp_query_required("spall.small",value.spall.small); // regularization to set minimum modulus multiplier

        pp_query_required("eta_ref_threshold", value.m_eta_ref_threshold); // Refinement threshold based on eta
        //pp.query("alpha", value.alpha); // :math:`\alpha` parameter
        //pp.query("beta", value.beta); // :math:`\beta` parameter
        //pp.query("gamma", value.gamma); // :math:`\gamma` parameter
        //pp.queryclass("L", value.L); // Mobility

        //pp.queryclass("lambda", value.lambda); // Lagrange multiplier (can be interplated)


        //value.RegisterIntegratedVariable(&value.volume, "volume");
        //value.RegisterIntegratedVariable(&value.w_chem_potential, "chem_potential");
        //value.RegisterIntegratedVariable(&value.w_bndry, "bndry");
        //value.RegisterIntegratedVariable(&value.w_elastic, "elastic");
        
        // std::string ictype;
        // pp_query_validate("ic.type",{"constant","bmp"});
        // if      (ictype == "constant") value.ic = new IC::Constant(value.geom, pp, "ic.constant");
        // else if (ictype == "bmp")      value.ic = new IC::BMP(value.geom, pp, "ic.bmp");
    }

    void Initialize(int lev) override
    {
        spall.ic->Initialize(lev,d_mf);
        spall.ic->Initialize(lev,d_old_mf);

        p_mf[lev]->setVal(0.0);
    }




    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    { 
        Base::Mechanics<MODEL>::TimeStepBegin(a_time, a_step);

    }

    virtual void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        for (int lev = 0; lev <= this->finest_level; ++lev)
            {
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());

                //this->rhs_mf[lev]->setVal(0.0);
                d_mf[lev]->FillBoundary();
                this->model_mf[lev]->setVal(this->model);
                Set::Vector DX(this->geom[lev].CellSize());

                for (MFIter mfi(*this->model_mf[lev], false); mfi.isValid(); ++mfi)
                    {
                        amrex::Box bx = mfi.grownnodaltilebox() & domain;
                        Set::Patch<MODEL> model_sum = model_mf.Patch(lev,mfi);
                        Set::Patch<const Set::Scalar> d = d_mf.Patch(lev,mfi);
                
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            Set::Scalar d_inv_avg = 1.0 - Numeric::Interpolate::CellToNodeAverage(d,i,j,k,0);
                            if (d_inv_avg > 1.0)   d_inv_avg = 1.0;
                            if (d_inv_avg < spall.small) d_inv_avg = spall.small;

                            model_sum(i, j, k).ddw = model.ddw * d_inv_avg;
                        });
                    }
                Util::RealFillBoundary(*this->model_mf[lev], this->geom[lev]);
            }

    }

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        BL_PROFILE("VoidPF::Advance");
        Base::Mechanics<MODEL>::Advance(lev, time, dt);
        std::swap(d_old_mf[lev], d_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        // Set::Scalar Lnow = L(time);
        //Set::Scalar lambdaT = lambda(time);

        for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            //bx.grow(1);
            //bx = bx & domain;
            
            Set::Patch<const Set::Matrix> sig  = stress_mf.Patch(lev,mfi);
            Set::Patch<const Set::Matrix> eps  = strain_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> dold = d_old_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar>       d    = d_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar>       p    = p_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar df = 0.0;

                df += spall.cp * 2.0 * dold(i,j,k)*(1 - 2.0*dold(i,j,k))*(1.0 - dold(i,j,k));
                df += spall.kappa * Numeric::Laplacian(dold,i,j,k,0,DX);

                Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k,0);
                Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k,0);
  
                p(i,j,k) = 0.5 * sig_avg.trace();
      
                if (p(i,j,k) > spall.po)
                    df += spall.elastic_mult*(p(i,j,k) - spall.po);

                d(i,j,k) = dold(i,j,k) - spall.M * dt * df;

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

            amrex::Array4<Set::Scalar> const& d = d_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                {
                    Set::Vector gradpsi = Numeric::Gradient(d, i, j, k, 0, DX.data(), sten);
                    if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
            
        }
    }


protected:
    MODEL model;
    IC::IC* ic;
    IC::IC* ic_euler = nullptr;
    BC::BC<Set::Scalar>*bc_d, *bc_p;
    Set::Scalar m_eta_ref_threshold = 1E-2;
    Set::Field<Set::Scalar> d_old_mf;
    Set::Field<Set::Scalar> d_mf;
    Set::Field<Set::Scalar> p_mf;


    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    //using Base::Mechanics<MODEL>::psi_on;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;

    //Set::Scalar alpha = 1.0;
    //Set::Scalar beta = 1.0;
    //Set::Scalar gamma = 1.0;
    //Set::Scalar L = 1.0;
    //Numeric::Interpolator::Linear<Set::Scalar> L;
    //Set::Scalar volume0 = 0.5;
    //Set::Scalar lambda = 1.0;
    //Numeric::Interpolator::Linear<Set::Scalar> lambda;

    //Set::Scalar volume = 0.0;
    //Set::Scalar w_chem_potential = 0.0;
    //Set::Scalar w_bndry = 0.0;
    //Set::Scalar w_elastic = 0.0;
    
    std::string filename; 
    
    struct {
        IC::IC *ic;
        Set::Scalar M = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar po = NAN;
        Set::Scalar cp = NAN;
        Set::Scalar elastic_mult = NAN;
        Set::Scalar small = 0.0;
    } spall;
};










} // namespace Integrator
#endif
