#ifndef INTEGRATOR_FINITEKINEMATICS_H
#define INTEGRATOR_FINITEKINEMATICS_H
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

#include "Solver/Nonlocal/Linear.H"

#include "IO/ParmParse.H"

#include "Model/Solid/Elastic/Elastic.H"
#include "Model/Solid/Elastic/NeoHookean.H"
#include "Model/Solid/Linear/Isotropic.H"

#include "Solver/Nonlocal/Newton.H"

namespace Integrator
{
class FiniteKinematics : public Integrator
{
    using model_type = Model::Solid::Elastic::NeoHookean;
    //using model_type = Model::Solid::Linear::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    FiniteKinematics()
    {
        RegisterNodalFab(eta_mf, 1, 2, "eta",true);
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(res_mf, AMREX_SPACEDIM, 2, "res",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        {
            IO::ParmParse pp("elastic");
            {
                pp.queryarr("b",b);
                pp.query("mu",mu);
                pp.query("kappa",kappa);
                pp.queryclass("bc",elastic.bc);
            }
        }
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        eta_mf[lev]->setVal(0.0);
        //ic->Initialize(lev, eta_mf);

        disp_mf[lev]->setVal(0.);
        rhs_mf[lev]->setVal(0.);
    }

    void TimeStepBegin(Set::Scalar a_time, int iter) override
    {
        Util::Message(INFO, "finest_level = ", finest_level);
        Util::Message(INFO, "rhs_mf.size() = ", rhs_mf.size());

        if (finest_level != rhs_mf.size() - 1)
        {
            Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
        }
        for (int lev = 0; lev < rhs_mf.size(); lev++) rhs_mf[lev]->setVal(0.0);
        for (int lev = 0; lev < rhs_mf.size(); lev++) disp_mf[lev]->setVal(0.0);

        Operator::Elastic<model_type> op;
        op.SetUniform(false);
        amrex::LPInfo info;
        //info.setMaxCoarseningLevel(0);
        op.define(geom, grids, dmap, info);

        // Set linear elastic model
        model_type mymodel;//(shear, lame, Set::Matrix::Zero());
        mymodel.mu = mu; mymodel.kappa = kappa;
        
        
        //amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type>>> model_mf;
        Set::Field<model_type> model_mf;
        model_mf.resize(disp_mf.size());
        for (int lev = 0; lev < rhs_mf.size(); ++lev)
        {
            amrex::Box domain(geom[lev].Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            model_mf[lev].reset(new amrex::FabArray<amrex::BaseFab<model_type>>(disp_mf[lev]->boxArray(), disp_mf[lev]->DistributionMap(), 1, 2));
            model_mf[lev]->setVal(mymodel);

            if (iter == 0) disp_mf[lev]->setVal(0.0);

            eta_mf[lev]->FillBoundary();
            rhs_mf[lev]->setVal(b(0),0,1,2);
            rhs_mf[lev]->setVal(b(1),1,1,2);
            rhs_mf[lev]->setVal(b(2),2,1,2);
        }
        
        elastic.bc.SetTime(a_time);
        elastic.bc.Init(rhs_mf,geom);
        op.SetBC(&elastic.bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
        Solver::Nonlocal::Newton<model_type> nr(op);
        IO::ParmParse mypp("elastic");
        mypp.queryclass("newton",nr);
        
        nr.solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);

        nr.compResidual(res_mf,disp_mf,rhs_mf,model_mf);

        for (int lev = 0; lev < disp_mf.size(); lev++)
        {
            op.Stress(lev, *stress_mf[lev], *disp_mf[lev],false,true);
            op.Strain(lev, *strain_mf[lev], *disp_mf[lev],false);
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

    model_type mymodel;

    amrex::Real alpha = 1.0;                 ///< Thermal diffusivity
    amrex::Real refinement_threshold = 0.01; ///< Criterion for cell refinement

    std::string ic_type; ///< String to determine what kind of initial condition to use
    //IC::IC *ic;          ///< Pointer to abstract IC object
    BC::BC *mybc;        ///< Stock generic boundary condition object
    
    Set::Scalar lame = 121.154, shear = 80.7692;
    Set::Scalar kappa, mu;

    Set::Matrix Fmatrix = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();
    Set::Scalar disp = 0.0;
    Set::Vector b = Set::Vector::Zero();
    
    struct {
        BC::Operator::Elastic<model_type> bc;
    } elastic;
    

    int fixediter = 100;
    int m_nriters = 2;

};
} // namespace Integrator
#endif
