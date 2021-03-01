#ifndef INTEGRATOR_TENSIONTEST_H
#define INTEGRATOR_TENSIONTEST_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"

#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"

#include "IO/ParmParse.H"

namespace Integrator
{
template<class MODEL>
class TensionTest : public Integrator
{
    //using model_type = Model::Solid::Affine::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    TensionTest()
    {
        
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        

        {
            IO::ParmParse pp("elastic");
            int nmodels = 1;
            pp.query("nmodels",nmodels);
            for (int i = 0; i<nmodels; i++)
            {
                std::string name = "model" + std::to_string(i+1);
                MODEL tmp_model;
                pp.queryclass(name.data(), tmp_model);
                elastic.models.push_back(tmp_model);
            }
            Util::Assert(INFO,TEST(elastic.models.size() > 0));

            pp.queryclass("bc",elastic.bc);
            RegisterNodalFab(eta_mf, elastic.models.size(), 2, "eta",true);    
        }
        RegisterGeneralFab(model_mf, 1, 2);

        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                ic = new IC::Ellipse(geom);
                pp.queryclass("ellipse",static_cast<IC::Ellipse*>(ic));
            }
        }
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        eta_mf[lev]->setVal(0.0);
        if (elastic.models.size() > 1 && ic)
        {
            ic->Initialize(lev, eta_mf);
        }
        else 
        {
            eta_mf[lev]->setVal(1.0);
        }

        disp_mf[lev]->setVal(0.);
        rhs_mf[lev]->setVal(0.);
    }

    void TimeStepBegin(Set::Scalar, int) override
    {
        // Set linear elastic model
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            rhs_mf[lev]->setVal(0.0);

            eta_mf[lev]->FillBoundary();

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();

                amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    model(i,j,k) = MODEL::Zero();
                    for (int n = 0; n < elastic.models.size(); n++) 
                        model(i,j,k) += eta(i,j,k,n) * elastic.models[n];
                });
            }

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }

        elastic.bc.Init(rhs_mf,geom);

        amrex::LPInfo info;
        Operator::Elastic<MODEL::sym> elastic_op(Geom(0,finest_level), grids, DistributionMap(0,finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetBC(&elastic.bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        IO::ParmParse pp("elastic");
        elastic.solver = new Solver::Nonlocal::Newton<MODEL>(elastic_op);
        pp.queryclass("solver",*elastic.solver);

        elastic.solver->solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);

        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            const amrex::Real* DX = geom[lev].CellSize();
            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                amrex::Array4<MODEL>        const &model = model_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>       const &stress  = stress_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>       const &strain  = strain_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &disp  = disp_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix gradu = Numeric::Gradient(disp,i,j,k,DX);
                    Set::Matrix sigma = model(i,j,k).DW(gradu);
                    Set::Matrix eps = 0.5*(gradu + gradu.transpose());
                    Numeric::MatrixToField(stress,i,j,k,sigma);
                    Numeric::MatrixToField(strain,i,j,k,eps);
                });
            }
        }
    }

    void Advance(int /*lev*/, Set::Scalar /*time*/, Set::Scalar /*dt*/) override
    {
        // Nothing to do here.
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        a_tags.setVal(amrex::TagBox::CLEAR);
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &disp = eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
                {
                    Set::Vector grad = Numeric::Gradient(disp, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
        }
    }

private:
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> disp_mf;
    Set::Field<Set::Scalar> rhs_mf;
    Set::Field<Set::Scalar> res_mf;
    Set::Field<Set::Scalar> stress_mf;
    Set::Field<Set::Scalar> strain_mf;

    Set::Field<MODEL> model_mf;

    IC::IC *ic = nullptr;          ///< Pointer to abstract IC object
    BC::BC<Set::Scalar> *mybc;        ///< Stock generic boundary condition object
    
    struct {
        std::vector<MODEL> models;
        model_type model1, model2;
        Solver::Nonlocal::Newton<MODEL> *solver;//(elastic.op);
        BC::Operator::Elastic::Constant bc;
    } elastic;

};
} // namespace Integrator
#endif
