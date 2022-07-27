#ifndef INTEGRATOR_BASE_MECHANICS_H
#define INTEGRATOR_BASE_MECHANICS_H

#include "AMReX.H"
#include "BC/Operator/Elastic/Elastic.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"
#include "Integrator/Integrator.H"
#include "Numeric/Stencil.H"
#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Operator/Operator.H"

namespace Integrator
{
namespace Base
{
template<class MODEL>
class Mechanics : virtual public Integrator
{
public:

    enum Type{Static, Dynamic, Disable};
    
    Mechanics() {}

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(Mechanics &value, IO::ParmParse &pp)
    {
        BL_PROFILE("Integrator::Base::Mechanics::Parse");
        if (pp.contains("type"))
        {
            std::string type_str = "static";
            pp.query("type",type_str);
            if (type_str == "static")        value.m_type = Type::Static;
            else if (type_str == "dynamic")  value.m_type = Type::Dynamic;
            else if (type_str == "disable")  value.m_type = Type::Disable;
            else Util::Abort(INFO,"Invalid type ", type_str, " specified");
        }
        if (value.m_type == Type::Disable) return;

        pp.query("time_evolving",value.m_time_evolving);
        
        value.RegisterGeneralFab(value.disp_mf, 1, 2, "disp",value.m_time_evolving);
        value.RegisterGeneralFab(value.rhs_mf, 1, 2, "rhs",value.m_time_evolving);
        value.RegisterGeneralFab(value.stress_mf, 1, 2, "stress",value.m_time_evolving);
        value.RegisterGeneralFab(value.strain_mf, 1, 2, "strain",value.m_time_evolving);

        if (value.m_type == Type::Static)
        {
            // Read parameters for :ref:`Solver::Nonlocal::Newton` solver
            pp.queryclass("solver",value.solver);
        }
        if (value.m_type == Type::Dynamic)
        {
                value.RegisterGeneralFab(value.vel_mf,1,2,"vel");
                value.RegisterGeneralFab(value.disp_old_mf,1,2,"dispold");
                value.RegisterGeneralFab(value.vel_old_mf,1,2,"velold");
                pp.query("viscous.mu",value.mu);            
        }

        std::string bc_type = "constant";
        // Determine the boundary condition type (contant, tension_test, expression)
        pp.query("bc.type",bc_type);
        Util::Message(INFO,bc_type);
        if (bc_type == "constant")          value.bc = new BC::Operator::Elastic::Constant(pp,"bc.constant");
        else if (bc_type == "tension_test") value.bc = new BC::Operator::Elastic::TensionTest(pp,"bc.tension_test");
        else if (bc_type == "expression")   value.bc = new BC::Operator::Elastic::Expression(pp,"bc.expression");
        else Util::Abort(INFO,"Invalid bc.type ",bc_type);
        
        pp.query("print_model",value.m_print_model);
        if (value.m_print_model) value.RegisterGeneralFab(value.model_mf, 1, 2,"model",value.m_time_evolving);
        else value.RegisterGeneralFab(value.model_mf, 1, 2,value.m_time_evolving);

        // Read in IC for RHS
        std::string rhstype;
        pp.query("rhs.type",rhstype); // Initializer for RHS
        if (rhstype == "trig") value.ic_rhs = new IC::Trig(value.geom,pp,"rhs.trig");

        // Timestep interval for elastic solves (default - solve every time)
        pp.query("interval",value.m_interval);

        value.RegisterIntegratedVariable(&(value.disp_hi[0].data()[0]),"disp_xhi_x");
        value.RegisterIntegratedVariable(&(value.trac_hi[0].data()[0]),"trac_xhi_x");

        // Maximum multigrid coarsening level (default - none, maximum coarsening)
        pp.query("max_coarsening_level",value.m_max_coarsening_level);

        pp.query("print_residual",value.m_print_residual);
        if (value.m_print_residual) value.RegisterGeneralFab(value.res_mf,1,2,"res",false);

        pp.query("elastic_ref_threshold",value.m_elastic_ref_threshold); // Whether to refine based on elastic solution
        
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Initialize");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;

        disp_mf[lev]->setVal(Set::Vector::Zero());

        if (ic_rhs) ic_rhs->Initialize(lev,rhs_mf);
        else rhs_mf[lev]->setVal(Set::Vector::Zero());
    }

    virtual void UpdateModel(int a_step) = 0;

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::TimeStepBegin");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;
        if (m_interval && a_step%m_interval) return;


        if (!m_type == Mechanics<MODEL>::Type::Static) return;

        for (int lev=0; lev <= finest_level; ++lev) 
        {
            rhs_mf[lev]->setVal(Set::Vector::Zero());
            if (ic_rhs) ic_rhs->Initialize(lev,rhs_mf);
        }

        UpdateModel(a_step);        

        bc->SetTime(a_time);
        bc->Init(rhs_mf,geom);

        amrex::LPInfo info;
        if (m_max_coarsening_level >= 0)
            info.setMaxCoarseningLevel(m_max_coarsening_level);
        Operator::Elastic<MODEL::sym> elastic_op(Geom(0,finest_level), grids, DistributionMap(0,finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetHomogeneous(false);
        elastic_op.SetBC(bc);
        IO::ParmParse pp("elasticop");
        pp.queryclass(elastic_op);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        solver.Define(elastic_op);
        if (psi_on) solver.setPsi(psi_mf);
        solver.solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);
        if (m_print_residual) solver.compLinearSolverResidual(res_mf,disp_mf,rhs_mf);
        solver.Clear();

        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());

            const amrex::Real* DX = geom[lev].CellSize();
            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                bx.grow(1);
                bx = bx & domain;
                amrex::Array4<MODEL>             const &model   = model_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const &stress  = stress_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const &strain  = strain_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector> const &disp  = disp_mf[lev]->array(mfi);


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                {
                    auto sten = Numeric::GetStencil(i,j,k,bx);
                    if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp,i,j,k,DX,sten);
                        stress(i,j,k) = model(i,j,k).DW(F);
                        strain(i,j,k) = F;
                    }
                    else
                    {                    
                        Set::Matrix gradu = Numeric::Gradient(disp,i,j,k,DX,sten);
                        stress(i,j,k) = model(i,j,k).DW(gradu);
                        strain(i,j,k) = 0.5*(gradu + gradu.transpose());
                    }
                });
            }
            Util::RealFillBoundary(*stress_mf[lev],geom[lev]);
            Util::RealFillBoundary(*strain_mf[lev],geom[lev]);
        }
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Advance");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;
        const amrex::Real *DX = geom[lev].CellSize();

        if (m_type == Type::Dynamic)
        {
            std::swap(*disp_mf[lev], *disp_old_mf[lev]);
            std::swap(*vel_mf[lev], *vel_old_mf[lev]);

            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            domain.grow(-1);

            for (amrex::MFIter mfi(*disp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Array4<Set::Vector> const &unew = (*disp_mf[lev]).array(mfi);
                amrex::Array4<Set::Vector> const &vnew = (*vel_mf[lev]).array(mfi);

                amrex::Array4<const Set::Vector> const &u = (*disp_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Vector> const &v = (*vel_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Matrix> const &eps = (*strain_mf[lev]).array(mfi);
                amrex::Array4<Set::Matrix> const &sig = (*stress_mf[lev]).array(mfi);
                amrex::Array4<const Set::Vector> const &b = (*rhs_mf[lev]).array(mfi);

                amrex::Array4<MODEL> const &model = (*model_mf[lev]).array(mfi);

                // Right now, we are forcing Dirichlet conditions on all boundaries
                amrex::Box bx = mfi.nodaltilebox().grow(1) & domain;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    //auto sten = Numeric::GetStencil(i,j,k,bx);
                    eps(i,j,k) = Numeric::Gradient(u,i,j,k,DX);
                    sig(i,j,k) = model(i,j,k).DW(eps(i,j,k)); 
                });

                bx = mfi.nodaltilebox().grow(0) & domain;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    //auto sten = Numeric::GetStencil(i,j,k,bx);
                    Set::Vector udotdot = Numeric::Divergence(sig,i,j,k,DX) + b(i,j,k) - mu * v(i,j,k);
                    udotdot /= rho;
                    vnew(i,j,k) = v(i,j,k) + dt*udotdot;
                    unew(i,j,k) = u(i,j,k) + dt*vnew(i,j,k); 
                });
            }
        }

        for (amrex::MFIter mfi(*disp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<Set::Matrix>         const &eps  = (*strain_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const &sig  = (*stress_mf[lev]).array(mfi);
            amrex::Array4<MODEL>               const &model= (*model_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                model(i,j,k).Advance(dt,eps(i,j,k),sig(i,j,k));
            });
        }
        
    }

    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
                    const amrex::MFIter &mfi, const amrex::Box &a_box) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Integrate");
        if (m_type==Type::Disable) return;

        const amrex::Real *DX = geom[amrlev].CellSize();
        amrex::Box domain = geom[amrlev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

        amrex::Box box = a_box;
        box.convert(amrex::IntVect::TheNodeVector());


        //Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        #if AMREX_SPACEDIM == 2
        Set::Vector da(DX[1],0);
        #elif AMREX_SPACEDIM == 3
        Set::Vector da(DX[1]*DX[2], 0, 0);
        #endif

        const Dim3 /*lo= amrex::lbound(domain),*/ hi = amrex::ubound(domain);
        const Dim3 /*boxlo= amrex::lbound(box),*/ boxhi = amrex::ubound(box);

        amrex::Array4<const Set::Matrix> const &stress = (*stress_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Vector> const &disp   = (*disp_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
        {
            #if AMREX_SPACEDIM == 2
            if (i == hi.x && j < boxhi.y)
            {
                trac_hi[0]+= (0.5 * (stress(i,j,k) + stress(i,j+1,k)) * da);
                disp_hi[0] = disp(i,j,k); 
            } 
            #elif AMREX_SPACEDIM == 3
            if (i == hi.x && (j < boxhi.y && k < boxhi.z))
            {
                trac_hi[0] += (0.25 * (stress(i,j,k) + stress(i,j+1,k)
                                       + stress(i,j,k+1) + stress(i,j+1,k+1)) * da);
                disp_hi[0] = disp(i,j,k);
            } 
            #endif
        });
            
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar /*time*/, int /*ngrow*/) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::TagCellsForRefinement");
        if (m_type == Type::Disable) return;

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*strain_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Matrix> const &eps = strain_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                auto sten = Numeric::GetStencil(i,j,k,bx);
                Set::Matrix3 grad = Numeric::Gradient(eps, i, j, k, DX.data(),sten);
                if (grad.norm() * DXnorm > m_elastic_ref_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<MODEL> model_mf;
    Set::Field<Set::Scalar> psi_mf;
    bool psi_on = false;

    int m_interval = 0;
    Type m_type = Type::Static;

    Set::Field<Set::Vector> disp_mf;
    Set::Field<Set::Vector> rhs_mf;
    Set::Field<Set::Vector> res_mf;
    Set::Field<Set::Matrix> stress_mf;
    Set::Field<Set::Matrix> strain_mf;
    
    // Only use these if using the "dynamics" option
    Set::Field<Set::Vector> disp_old_mf;
    Set::Field<Set::Vector> vel_mf;
    Set::Field<Set::Vector> vel_old_mf;
    //Set::Field<Set::Matrix4<AMREX_SPACEDIM,MODEL::sym>> ddw_mf;
    Set::Scalar rho = 1.0;
    Set::Scalar mu = 0.0;

    //Set::Vector trac_lo[AMREX_SPACEDIM];
    Set::Vector trac_hi[AMREX_SPACEDIM];
    Set::Vector disp_hi[AMREX_SPACEDIM];


    IC::IC *ic_rhs = nullptr;
    BC::BC<Set::Scalar> *mybc;
    
    Solver::Nonlocal::Newton<MODEL> solver;//(elastic.op);
    BC::Operator::Elastic::Elastic *bc;

    Set::Scalar m_elastic_ref_threshold = 0.01;
    bool m_print_model = false;
    bool m_print_residual = false;
    bool m_time_evolving = false;
    int m_max_coarsening_level = -1;
};
}
}
#endif