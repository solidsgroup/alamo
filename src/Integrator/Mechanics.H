#ifndef INTEGRATOR_MECHANICS_H
#define INTEGRATOR_MECHANICS_H
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
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"


namespace Integrator
{
template<class MODEL>
class Mechanics : virtual public Base::Mechanics<MODEL>
{
    public:
    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(Mechanics &value, IO::ParmParse &pp)
    {
        Base::Mechanics<MODEL>::Parse(value,pp);
        int nmodels = 1;
        pp.query("nmodels",nmodels); // Number of elastic model varieties
        for (int i = 0; i<nmodels; i++)
        {
            std::string name = "model" + std::to_string(i+1);
            MODEL tmp_model;
            pp.queryclass(std::string(name.data()), tmp_model);
            value.models.push_back(tmp_model);
        }
        Util::Assert(INFO,TEST(value.models.size() > 0));
        value.RegisterNodalFab(value.eta_mf, value.models.size(), 2, "eta", true);    
        // Refinement threshold for eta field
        pp.query("eta_ref_threshold",value.m_eta_ref_threshold);
        // Refinement threshold for strain gradient
        pp.query("ref_threshold",value.m_elastic_ref_threshold);


        // Read in IC for eta
        {
            std::string type;
            pp.query("ic.type",type); // Read IC type for the eta field
            if (type=="ellipse") value.ic_eta = new IC::Ellipse(value.geom,pp,"ic.ellipse");
            if (type=="voronoi") value.ic_eta = new IC::Voronoi(value.geom,pp,"ic.voronoi");
        }

        // Read in IC for psi
        if (pp.contains("psi.ic.type"))
        {
            std::string type;
            pp.query("psi.ic.type",type); // Read IC type for the eta field
            if (type=="ellipse") value.ic_psi = new IC::Ellipse(value.geom,pp,"psi.ic.ellipse");
            else if (type=="constant") value.ic_psi = new IC::Constant(value.geom,pp,"psi.ic.constant");
            else Util::Abort(INFO,"Invalid value for psi.ic.type: ",type);
            
            value.bc_psi = new BC::Nothing();
            value.RegisterNewFab(value.psi_mf,value.bc_psi,1,2,"psi",true);
            value.psi_on = true;
        }
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        eta_mf[lev]->setVal(0.0);
        if (models.size() > 1 && ic_eta) ic_eta->Initialize(lev, eta_mf);
        else eta_mf[lev]->setVal(1.0);

        if (psi_on) ic_psi->Initialize(lev,psi_mf);
    }

    virtual void UpdateModel(int a_step) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (a_step > 0) return;
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            //rhs_mf[lev]->setVal(0.0);

            eta_mf[lev]->FillBoundary();

            //disp_mf[lev]->setVal(Set::Vector::Zero());

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();

                amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                            model(i,j,k) = MODEL::Zero();
                                            for (unsigned int n = 0; n < models.size(); n++) 
                                                model(i,j,k) += eta(i,j,k,n) * models[n];
                                        });
            }

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev,a_tags,a_time,a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                auto sten = Numeric::GetStencil(i,j,k,bx);
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data(),sten);
                    if (grad.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
            if (psi_on)
            {
                amrex::Array4<Set::Scalar> const &psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                {
                    auto sten = Numeric::GetStencil(i,j,k,bx);
                    {
                        Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(),sten);
                        if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }            
        }
    }


    protected:
    Set::Field<Set::Scalar> eta_mf;
    Set::Scalar m_eta_ref_threshold = 0.01;
    std::vector<MODEL> models;
    IC::IC *ic_eta = nullptr;
    IC::IC *ic_psi = nullptr;
    BC::BC<Set::Scalar> *bc_psi = nullptr;

    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;
};










} // namespace Integrator
#endif
