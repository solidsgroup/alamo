#ifndef INTEGRATOR_MECHANICS_H
#define INTEGRATOR_MECHANICS_H
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
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

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
class Mechanics : virtual public Integrator
{
    //using model_type = Model::Solid::Affine::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    Mechanics()
    {
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(Mechanics &value, IO::ParmParse &pp)
    {

        int nmodels = 1;
        pp.query("elastic.nmodels",nmodels); // Number of elastic model varieties
        for (int i = 0; i<nmodels; i++)
        {
            std::string name = "model" + std::to_string(i+1);
            MODEL tmp_model;
            pp.queryclass("elastic." + std::string(name.data()), tmp_model);
            value.elastic.models.push_back(tmp_model);
        }

        Util::Assert(INFO,TEST(value.elastic.models.size() > 0));
        value.RegisterNodalFab(value.eta_mf, value.elastic.models.size(), 2, "eta", true);    

        std::string bc_type = "constant";
        // Determine the boundary condition type (contant, tension_test, expression)
        pp.query("elastic.bc.type",bc_type);
        Util::Message(INFO,bc_type);
        if (bc_type == "constant")
        {
            value.elastic.bc = new BC::Operator::Elastic::Constant();
            // Query :ref:`BC::Operator::Elastic::Constant`
            pp.queryclass("elastic.bc",static_cast<BC::Operator::Elastic::Constant*>(value.elastic.bc)); 
        }
        else if (bc_type == "tension_test")
        {
            value.elastic.bc = new BC::Operator::Elastic::TensionTest();
            // Query :ref:`BC::Operator::Elastic::TensionTest`
            pp.queryclass("elastic.bc.tension_test",static_cast<BC::Operator::Elastic::TensionTest*>(value.elastic.bc)); 
        }
        else if (bc_type == "expression")
        {
            value.elastic.bc = new BC::Operator::Elastic::Expression();
            // Query :ref:`BC::Operator::Elastic::Expression`
            pp.queryclass("elastic.bc.expression",static_cast<BC::Operator::Elastic::Expression*>(value.elastic.bc));
        }
        // Refinement threshold for eta field
        pp.query("elastic.eta_ref_threshold",value.m_eta_ref_threshold);
        // Refinement threshold for strain gradient
        pp.query("elastic.ref_threshold",value.m_elastic_ref_threshold);

        value.RegisterGeneralFab(value.model_mf, 1, 2);

        // Read in IC for eta
        std::string type;
        pp.query("ic.type",type); // Read IC type for the eta field
        Util::Message(INFO,type);
        if (type=="ellipse")
        {
            value.ic = new IC::Ellipse(value.geom);
            // Initialize with :ref:`IC::Ellipse`
            pp.queryclass("ic.ellipse",static_cast<IC::Ellipse*>(value.ic));
        }

        // Read in IC for RHS
        std::string rhstype;
        pp.query("rhs.type",rhstype); // Initializer for RHS
        if (rhstype == "trig")
        {
            value.ic_rhs = new IC::Trig(value.geom);
            // See :ref:`IC::Trig`
            pp.queryclass("rhs.trig",static_cast<IC::Trig*>(value.ic_rhs));
        }

        // Timestep interval for elastic solves (default - solve every time)
        pp.query("interval",value.m_interval);

#if AMREX_SPACEDIM==3
        value.RegisterIntegratedVariable(&(value.trac_lo[0].data()[0]),"trac_lo_xx");
#endif

        // Read parameters for :ref:`Solver::Nonlocal::Newton` solver
        pp.queryclass("elastic.solver",value.elastic.solver);

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

        if (ic_rhs) ic_rhs->Initialize(lev,rhs_mf);
        else rhs_mf[lev]->setVal(0.);
    }

    virtual void UpdateModel()
    {
        // Set linear elastic model
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            //rhs_mf[lev]->setVal(0.0);

            eta_mf[lev]->FillBoundary();

            disp_mf[lev]->setVal(0.0);

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();

                amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                            model(i,j,k) = MODEL::Zero();
                                            for (unsigned int n = 0; n < elastic.models.size(); n++) 
                                                model(i,j,k) += eta(i,j,k,n) * elastic.models[n];
                                        });
            }

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {

        if (m_interval && a_step%m_interval) return;

        UpdateModel();        

        elastic.bc->SetTime(a_time);
        elastic.bc->Init(rhs_mf,geom);

        amrex::LPInfo info;
        Operator::Elastic<MODEL::sym> elastic_op(Geom(0,finest_level), grids, DistributionMap(0,finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetBC(elastic.bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        elastic.solver.Define(elastic_op);
        elastic.solver.solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);
        elastic.solver.Clear();

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
                                            std::array<Numeric::StencilType,AMREX_SPACEDIM>
                                                sten = Numeric::GetStencil(i,j,k,bx);

                                            if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                                            {
                                                Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp,i,j,k,DX,sten);
                                                Set::Matrix P = model(i,j,k).DW(F);
                                                Numeric::MatrixToField(stress,i,j,k,P);
                                                Numeric::MatrixToField(strain,i,j,k,F);
                                            }
                                            else
                                            {                    
                                                Set::Matrix gradu = Numeric::Gradient(disp,i,j,k,DX,sten);
                                                Set::Matrix sigma = model(i,j,k).DW(gradu);
                                                Set::Matrix eps = 0.5*(gradu + gradu.transpose());
                                                Numeric::MatrixToField(stress,i,j,k,sigma);
                                                Numeric::MatrixToField(strain,i,j,k,eps);
                                            }
                                        });
            }
            Util::RealFillBoundary(*stress_mf[lev],geom[lev]);
            Util::RealFillBoundary(*disp_mf[lev],geom[lev]);
        }
    }

    void Advance(int /*lev*/, Set::Scalar /*time*/, Set::Scalar /*dt*/) override
    {

    }

    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
                    const amrex::MFIter &mfi, const amrex::Box &box) override
    {
        const amrex::Real *DX = geom[amrlev].CellSize();
        amrex::Box domain = geom[amrlev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());


        //Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
#if AMREX_SPACEDIM == 2
        Set::Vector da(DX[1],DX[0]);
#elif AMREX_SPACEDIM == 3
        Set::Vector da(DX[1]*DX[2], DX[2]*DX[0], DX[0]*DX[1]);
#endif

        const Dim3 lo= amrex::lbound(domain)/*, hi = amrex::ubound(domain)*/;
        const Dim3 /*boxlo= amrex::lbound(box),*/ boxhi = amrex::ubound(box);

        amrex::Array4<const amrex::Real> const &stress = (*stress_mf[amrlev]).array(mfi);
        //amrex::Array4<const amrex::Real> const &disp   = (*disp_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

#if AMREX_SPACEDIM == 2
                                    if (i == lo.x && j < boxhi.y)
                                    {
                                        trac_lo[0](0) += 0.5 * (stress(i,j,k,0) + stress(i,j+1,k,0)) * da(0);
                                    } 
#elif AMREX_SPACEDIM == 3
                                    if (i == lo.x && (j < boxhi.y && k < boxhi.z))
                                    {
                                        trac_lo[0](0) += 0.25 * (stress(i,j,k,0) + stress(i,j+1,k,0) + stress(i,j,k+1,0) + stress(i,j+1,k+1,0)) * da(0);
                                    } 
#endif

                                });
            
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        //a_tags.setVal(amrex::TagBox::CLEAR);
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &strain = strain_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        {
                                            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                                            if (grad.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                                                tags(i, j, k) = amrex::TagBox::SET;
                                        }
                                        {
                                            Set::Matrix3 grad = Numeric::MatrixGradient(strain, i, j, k, DX.data());
                                            if (grad.norm() * DXnorm > m_elastic_ref_threshold)
                                                tags(i, j, k) = amrex::TagBox::SET;
                                        }
                                    });
        }
    }

protected:
    Set::Field<MODEL> model_mf;
    Set::Field<Set::Scalar> eta_mf;

private:
    int m_interval = 0;

    Set::Field<Set::Scalar> disp_mf;
    Set::Field<Set::Scalar> rhs_mf;
    Set::Field<Set::Scalar> res_mf;
    Set::Field<Set::Scalar> stress_mf;
    Set::Field<Set::Scalar> strain_mf;

    Set::Vector trac_lo[AMREX_SPACEDIM];
    Set::Vector trac_hi[AMREX_SPACEDIM];



    IC::IC *ic = nullptr;          ///< Pointer to abstract IC object
    IC::IC *ic_rhs = nullptr;
    BC::BC<Set::Scalar> *mybc;        ///< Stock generic boundary condition object
    
    struct {
        std::vector<MODEL> models;
        Solver::Nonlocal::Newton<MODEL> solver;//(elastic.op);
        BC::Operator::Elastic::Elastic *bc;
    } elastic;

    Set::Scalar m_eta_ref_threshold = 0.01;
    Set::Scalar m_elastic_ref_threshold = 0.01;

};
} // namespace Integrator
#endif
