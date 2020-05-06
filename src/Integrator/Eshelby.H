#ifndef INTEGRATOR_ESHELBY_H
#define INTEGRATOR_ESHELBY_H
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
#include "IC/BMP.H"
#include "IC/Notch.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Affine/Isotropic.H"

#include "Operator/Operator.H"

#include "IO/ParmParse.H"

namespace Integrator
{
class Eshelby : public Integrator
{
    using model_type = Model::Solid::Affine::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    Eshelby()
    {
        RegisterNodalFab(eta_mf, 1, 2, "eta",true);
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        {
            IO::ParmParse pp("elastic");
            {
                // Initialize elastic model
                //if (pp.contains("model"))
                //{
                //    pp.queryclass("model",elastic.model1);
                //    pp.queryclass("model",elastic.model2);
                //}
                //else if (pp.contains("model1") && pp.contains("model2"))
                //{
                pp.queryclass("model1",elastic.model1);
                pp.queryclass("model2",elastic.model2);
                //}
                //else
                //{
                //    Util::Abort(INFO,"Did not find (model) or (model1 & model2) when querying");
                //}
                

                // TODO: figure out how to initialize solvers in
                // here instead of TimeStepBegin...

                pp.queryclass("bc",elastic.bc);
            }
        }

        RegisterNodalFab(c_mf, 1, 2, "c",true);
        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                ic = new IC::Ellipse(geom);
                pp.queryclass("ellipse",static_cast<IC::Ellipse*>(ic));
            }
            else if (type=="bmp")
            {
                ic = new IC::BMP(geom);
                pp.queryclass("bmp",static_cast<IC::BMP*>(ic));
            }
            if (pp.contains("Fmatrix"))    {pp.queryarr("Fmatrix",   elastic.model1.F0); Util::Warning(INFO,"This is going to be depricated"); }
            if (pp.contains("Finclusion")) {pp.queryarr("Finclusion",elastic.model2.F0); Util::Warning(INFO,"This is going to be depricated"); }
        }

        {
            IO::ParmParse pp("fracture");
            fracture.ic = new IC::Notch(geom);
            pp.queryclass("ic",static_cast<IC::Notch*>(fracture.ic));
            pp.query("lower_bound",fracture.lower_bound);
            pp.query("ref_threshold",fracture.ref_threshold);
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

        fracture.ic->Initialize(lev,c_mf);
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
            model_mf[lev]->setVal(elastic.model1);

            eta_mf[lev]->FillBoundary();

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.growntilebox(2);

                amrex::Array4<model_type>        const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &c   = c_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    model(i, j, k).ddw = eta(i, j, k) * elastic.model1.ddw + (1.0 - eta(i, j, k)) * elastic.model2.ddw;
                    model(i, j, k).ddw *= c(i,j,k) + fracture.lower_bound;
                    model(i, j, k).F0  = eta(i, j, k) * elastic.model1.F0  + (1.0 - eta(i, j, k)) * elastic.model2.F0;
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
        elastic.solver = new Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic>(elastic.op);
        pp.queryclass("solver",*elastic.solver);

        elastic.solver->solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);

        elastic.solver->DW(stress_mf,disp_mf,model_mf);
        for (int lev = 0; lev < disp_mf.size(); lev++)
        {
            //elastic.op.Stress(lev, *stress_mf[lev], *disp_mf[lev]);
            //elastic.op.Strain(lev, *strain_mf[lev], *disp_mf[lev]);
            strain_mf[lev]->setVal(0.0);
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
            amrex::Array4<Set::Scalar> const &eta  = eta_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &disp = disp_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &c    = c_mf[lev]->array(mfi);
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
                {
                    Set::Vector grad = Numeric::Gradient(c, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > fracture.ref_threshold)
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

    Set::Field<Set::Scalar> c_mf;

    IC::IC *ic;          ///< Pointer to abstract IC object
    BC::BC *mybc;        ///< Stock generic boundary condition object
    
    Set::Matrix Fmatrix = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();

    struct {
        model_type model1, model2;
        Operator::Elastic<Model::Solid::Affine::Isotropic> op;
        Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic> *solver;//(elastic.op);
        BC::Operator::Elastic<model_type> bc;
    } elastic;

    struct {
        IC::IC *ic;        
        Set::Scalar lower_bound = 0.01;
        Set::Scalar ref_threshold = 0.01;
    } fracture;

};
} // namespace Integrator
#endif
