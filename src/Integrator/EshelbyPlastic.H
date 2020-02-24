#ifndef INTEGRATOR_ESHELBYPLASTIC_H
#define INTEGRATOR_ESHELBYPLASTIC_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"

#include "Operator/Operator.H"

#include "IO/ParmParse.H"

namespace Integrator
{
class EshelbyPlastic : public Integrator
{
    using model_type = Model::Solid::CrystalPlastic::CrystalPlastic;
public:
    /// \brief Read in parameters and register field variables
    EshelbyPlastic()
    {
        RegisterNodalFab(eta_mf, 1, 2, "eta",true);
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                ic = new IC::Ellipse(geom); 
                pp.queryclass("ellipse",static_cast<IC::Ellipse*>(ic));
            }
            if (pp.contains("Finclusion")) pp.queryarr("Finclusion",Finclusion);
            if (pp.contains("Fmatrix"))    pp.queryarr("Fmatrix",Fmatrix);
        }
        {
            IO::ParmParse pp("elastic");
            {
                // Initialize elastic model
                pp.queryclass("model",elastic.model);
                elastic.model.esp = Set::Matrix::Zero();

                // TODO: figure out how to initialize solvers in
                // here instead of TimeStepBegin...

                pp.queryclass("bc",elastic.bc);
            }
        }
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        eta_mf[lev]->setVal(0.0);
        ic->Initialize(lev, eta_mf);

        disp_mf[lev]->setVal(0.);
        rhs_mf[lev]->setVal(0.);
    }

    void TimeStepBegin(Set::Scalar, int) override
    {
        if (finest_level != rhs_mf.size() - 1)
        {
            Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
        }
        
        // Set linear elastic model
        Set::Field<model_type> model_mf;
        model_mf.resize(disp_mf.size());
        for (int lev = 0; lev < rhs_mf.size(); ++lev)
        {
            rhs_mf[lev]->setVal(0.0);

            amrex::Box domain(geom[lev].Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            model_mf.Define(lev,disp_mf[lev]->boxArray(), disp_mf[lev]->DistributionMap(), 1, 2);
            model_mf[lev]->setVal(elastic.model);

            eta_mf[lev]->FillBoundary();

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.growntilebox(2);

                amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix Fgb = (1.0 - eta(i, j, k)) * Fmatrix + eta(i, j, k) * Finclusion;
                    model(i, j, k) = elastic.model;
                    model(i, j, k).esp = Fgb;
                });
            }

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }

        elastic.bc.Init(rhs_mf,geom);

        amrex::LPInfo info;
        elastic.op.define(geom, grids, dmap, info);
        elastic.op.SetUniform(false);
        elastic.op.SetModel(model_mf);
        elastic.op.SetBC(&elastic.bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        IO::ParmParse pp("elastic");
        elastic.solver = new Solver::Nonlocal::Newton<Model::Solid::CrystalPlastic::CrystalPlastic>(elastic.op);
        pp.queryclass("solver",*elastic.solver);

        elastic.solver->solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);

        for (int lev = 0; lev < disp_mf.size(); lev++)
        {
            elastic.op.Stress(lev, *stress_mf[lev], *disp_mf[lev]);
            elastic.op.Strain(lev, *strain_mf[lev], *disp_mf[lev]);
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

    IC::IC *ic;          ///< Pointer to abstract IC object
    BC::BC *mybc;        ///< Stock generic boundary condition object
    
    Set::Matrix Fmatrix = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();

    struct {
        Model::Solid::CrystalPlastic::CrystalPlastic model;
        Operator::Elastic<Model::Solid::CrystalPlastic::CrystalPlastic> op;
        Solver::Nonlocal::Newton<Model::Solid::CrystalPlastic::CrystalPlastic> *solver;//(elastic.op);
        BC::Operator::Elastic<model_type> bc;
    } elastic;

};
} // namespace Integrator
#endif