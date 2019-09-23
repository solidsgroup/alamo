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

#include "IC/HeatConduction/Cylinder.H"
#include "IC/HeatConduction/Constant.H"
#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/LinearElastic/MultiWell.H"
#include "Solver/Nonlocal/Linear.H"

#include "IO/ParmParse.H"

namespace Integrator
{
class Eshelby : public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    Eshelby()
    {
        RegisterNodalFab(eta_mf, 1, 2, "eta");
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp");
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs");
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress");
        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                Set::Matrix A = Set::Matrix::Identity();
                if (pp.contains("A")) pp.queryarr("A",A);
                Set::Vector x0 = Set::Vector::Zero();
                if (pp.contains("x0")) pp.queryarr("x0",x0);
                Set::Scalar eps = 0.01;
                pp.query("eps",eps);
                ic = new IC::Ellipse(geom, A, x0, IC::Ellipse::Mollifier::Gaussian, eps);

                if (pp.contains("Finclusion")) pp.queryarr("Finclusion",Finclusion);
                if (pp.contains("Fmatrix")) pp.queryarr("Fmatrix",Fmatrix);
            }
        }
        {
            IO::ParmParse pp("elastic");
            {
                if (pp.contains("lame") && pp.contains("shear"))
                {
                    pp.query("lame",lame);
                    pp.query("shear",shear);
                }
                else if (pp.contains("E") && pp.contains("nu"))
                {
                    Set::Scalar E, nu;
                    pp.query("E",E);
                    pp.query("nu",nu);
                    lame = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
                    shear = E / 2.0 / (1.0 + nu);
                    Util::Message(INFO,"lame = ", lame);
                    Util::Message(INFO,"shear = ", shear);
                }
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
        Util::Message(INFO, "finest_level = ", finest_level);
        Util::Message(INFO, "rhs_mf.size() = ", rhs_mf.size());

        if (finest_level != rhs_mf.size() - 1)
        {
            Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
        }
        for (int lev = 0; lev < rhs_mf.size(); lev++) rhs_mf[lev]->setVal(0.0);
        //for (int lev = 0; lev < rhs_mf.size(); lev++) disp_mf[lev]->setVal(0.0);

        using model_type = Model::Solid::LinearElastic::Multiwell;
        Operator::Elastic<model_type> elastic;
        elastic.SetUniform(false);
        amrex::LPInfo info;
        //info.setMaxCoarseningLevel(0);
        elastic.define(geom, grids, dmap, info);

        // Set linear elastic model
        model_type mymodel(lame, shear, Set::Matrix::Zero());
        amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type>>> model_mf;
        model_mf.resize(disp_mf.size());
        Finclusion(0, 0) = 0.1;
        for (int lev = 0; lev < rhs_mf.size(); ++lev)
        {
            amrex::Box domain(geom[lev].Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            model_mf[lev].define(disp_mf[lev]->boxArray(), disp_mf[lev]->DistributionMap(), 1, 2);
            model_mf[lev].setVal(mymodel);

            eta_mf[lev]->FillBoundary();

            Set::Vector DX(geom[lev].CellSize());

            //for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
            for (MFIter mfi(model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.growntilebox(2);

                amrex::Array4<model_type> const &model = model_mf[lev].array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix Fgb = (1.0 - eta(i, j, k)) * Fmatrix + eta(i, j, k) * Finclusion;
                    model(i, j, k) = model_type(lame, shear, Fgb);
                });

                bx = mfi.growntilebox(1);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    for (int p = 0; p < AMREX_SPACEDIM; p++)
                        for (int q = 0; q < AMREX_SPACEDIM; q++)
                        {
                            AMREX_D_TERM(model(i, j, k).gradFgb[p](q, 0) = ((model(i + 1, j, k).Fgb - model(i - 1, j, k).Fgb) / 2. / DX(0))(p, q);,
                                         model(i, j, k).gradFgb[p](q, 1) = ((model(i, j + 1, k).Fgb - model(i, j - 1, k).Fgb) / 2. / DX(1))(p, q);,
                                         model(i, j, k).gradFgb[p](q, 2) = ((model(i, j, k + 1).Fgb - model(i, j, k - 1).Fgb) / 2. / DX(2))(p, q);)
                        }
                });
            }

            amrex::Geometry geom = elastic.Geom(lev);
            for (int i = 0; i < 2; i++)
            {
                amrex::FabArray<amrex::BaseFab<model_type>> &mf = model_mf[lev];
                mf.FillBoundary(geom.periodicity());
                const int ncomp = mf.nComp();
                const int ng1 = 1;
                const int ng2 = 2;
                amrex::FabArray<amrex::BaseFab<model_type>> tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
                amrex::Copy(tmpmf, mf, 0, 0, ncomp, ng1);
                mf.ParallelCopy(tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
            }
        }
        elastic.SetModel(model_mf);

        BC::Operator::Elastic<model_type> bc;
        elastic.SetBC(&bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
        Solver::Nonlocal::Linear linearsolver(elastic);
        linearsolver.setVerbose(2);
        linearsolver.solveaffine(disp_mf, rhs_mf, tol_rel, tol_abs, true);

        for (int lev = 0; lev < disp_mf.size(); lev++)
        {
            elastic.Stress(lev, *stress_mf[lev], *disp_mf[lev]);
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
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> eta_mf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> disp_mf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs_mf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> res_mf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> stress_mf;

    amrex::Real alpha = 1.0;                 ///< Thermal diffusivity
    amrex::Real refinement_threshold = 0.01; ///< Criterion for cell refinement

    std::string ic_type; ///< String to determine what kind of initial condition to use
    IC::IC *ic;          ///< Pointer to abstract IC object
    BC::BC *mybc;        ///< Stock generic boundary condition object
    
    Set::Scalar lame = 2.6, shear = 6.0;

    Set::Matrix Fmatrix = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();

};
} // namespace Integrator
#endif
