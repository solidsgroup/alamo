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
#include "IC/Constant.H"
#include "IC/Expression.H"

namespace Integrator
{
namespace Base
{
template<class MODEL>
class Mechanics: virtual public Integrator
{
public:

    enum Type { Static, Dynamic, Disable };

    Mechanics() : Integrator()
    {}

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(Mechanics& value, IO::ParmParse& pp)
    {
        BL_PROFILE("Integrator::Base::Mechanics::Parse");
        if (pp.contains("type"))
        {
            std::string type_str;
            // Type of mecahnics to use.
            // Static: do full implicit solve.
            // Dynamic: evolve dynamic equations with explicit dynamics
            // Disable: do nothing.
            pp_query_validate("type", type_str, {"disable","static","dynamic"}); 
            if (type_str == "static")        value.m_type = Type::Static;
            else if (type_str == "dynamic")  value.m_type = Type::Dynamic;
            else if (type_str == "disable")  value.m_type = Type::Disable;
            else Util::Abort(INFO, "Invalid type ", type_str, " specified");
        }
        if (value.m_type == Type::Disable) return;

        // Treat mechanics fields as changing in time. [false]
        // You should use this if you care about other physics driven by
        // the output of this integrator.
        pp_query_default("time_evolving", value.m_time_evolving,false);

        pp_query_default("plot_disp", value.plot_disp, true); // Include displacement field in output
        pp_query_default("plot_rhs", value.plot_rhs, true); // Include right-hand side in output
        pp_query_default("plot_psi", value.plot_psi, true); // Include :math:`\psi` field in output
        pp_query_default("plot_stress", value.plot_stress, true); // Include stress in output
        pp_query_default("plot_strain", value.plot_strain, true); // Include strain in output

        value.RegisterGeneralFab(value.disp_mf, 1, 2, value.plot_disp, "disp", value.m_time_evolving);
        value.RegisterGeneralFab(value.rhs_mf, 1, 2, value.plot_rhs, "rhs", value.m_time_evolving);
        value.RegisterGeneralFab(value.stress_mf, 1, 2, value.plot_stress, "stress", value.m_time_evolving);
        value.RegisterGeneralFab(value.strain_mf, 1, 2, value.plot_strain, "strain", value.m_time_evolving);

        if (value.m_type == Type::Static)
        {
            // Read parameters for :ref:`Solver::Nonlocal::Newton` solver
            pp_queryclass("solver", value.solver);
        }
        if (value.m_type == Type::Dynamic)
        {
            value.RegisterGeneralFab(value.vel_mf, 1, 2, "vel",true);
            value.RegisterGeneralFab(value.disp_old_mf, 1, 2, "dispold");
            value.RegisterGeneralFab(value.vel_old_mf, 1, 2, "velold");
            value.RegisterGeneralFab(value.ddw_mf, 1, 2);
            pp_forbid("viscous.mu","replaced with viscous.mu_dashpot")
            pp_query_default("viscous.mu_dashpot", value.mu_dashpot,0.0); // Dashpot damping (damps velocity)
            pp_forbid("viscous.mu2","replaced with viscous.mu_newton");
            pp.query_default("viscous.mu_newton", value.mu_newton,0.0); // Newtonian viscous damping (damps velocity gradient)

            std::string velocity_ic_str;
            pp_query_validate("velocity.ic.type",velocity_ic_str,{"none","expression"}); // Initializer for RHS
            if (velocity_ic_str == "expression") value.velocity_ic = new IC::Expression(value.geom, pp,"velocity.ic.expression");
        }

        std::string bc_type;
        // Determine the boundary condition type
        pp_query_validate("bc.type", bc_type,{"constant","tension_test","expression"});
        if (bc_type == "constant")          value.bc = new BC::Operator::Elastic::Constant(pp, "bc.constant");
        else if (bc_type == "tension_test") value.bc = new BC::Operator::Elastic::TensionTest(pp, "bc.tension_test");
        else if (bc_type == "expression")   value.bc = new BC::Operator::Elastic::Expression(pp, "bc.expression");
        else Util::Abort(INFO, "Invalid bc.type ", bc_type);

        pp_query_default("print_model", value.m_print_model, false); // Print out model variables (if enabled by model)
        if (value.m_print_model) value.RegisterGeneralFab(value.model_mf, 1, 2, "model", value.m_time_evolving);
        else value.RegisterGeneralFab(value.model_mf, 1, 2, value.m_time_evolving);

        // Read in IC for RHS
        std::string rhstype;
        pp_query_validate("rhs.type", rhstype,{"none","trig","constant"}); // Initializer for RHS
        if (rhstype == "trig") value.ic_rhs = new IC::Trig(value.geom, pp, "rhs.trig");
        else if (rhstype == "constant") value.ic_rhs = new IC::Constant(value.geom, pp, "rhs.constant");

        // Timestep interval for elastic solves (default - solve every time)
        pp_query_default("interval", value.m_interval, 0);

        value.RegisterIntegratedVariable(&(value.disp_hi[0].data()[0]), "disp_xhi_x");
        value.RegisterIntegratedVariable(&(value.disp_hi[0].data()[1]), "disp_xhi_y");
        value.RegisterIntegratedVariable(&(value.disp_hi[1].data()[0]), "disp_yhi_x");
        value.RegisterIntegratedVariable(&(value.disp_hi[1].data()[1]), "disp_yhi_y");
        value.RegisterIntegratedVariable(&(value.trac_hi[0].data()[0]), "trac_xhi_x");
        value.RegisterIntegratedVariable(&(value.trac_hi[0].data()[1]), "trac_xhi_y");
        value.RegisterIntegratedVariable(&(value.trac_hi[1].data()[0]), "trac_yhi_x");
        value.RegisterIntegratedVariable(&(value.trac_hi[1].data()[1]), "trac_yhi_y");

        // Maximum multigrid coarsening level (default - none, maximum coarsening)
        pp_query_default("max_coarsening_level", value.m_max_coarsening_level,-1);

        // Whether to include residual output field
        pp_query_default("print_residual", value.m_print_residual,false); 
        if (value.m_print_residual) value.RegisterGeneralFab(value.res_mf, 1, 2, "res", false);

        // Whether to refine based on elastic solution
        pp_query_default("elastic_ref_threshold", value.m_elastic_ref_threshold,0.01); 


        // Set this to true to zero out the displacement before each solve.
        // (This is a temporary fix - we need to figure out why this is needed.)
        pp_query_default("zero_out_displacement", value.m_zero_out_displacement,false);

        // Time to start doing the elastic solve (by default, start immediately)
        pp_query_default("tstart", value.tstart,-1.0);
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Initialize");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;

        disp_mf[lev]->setVal(Set::Vector::Zero());
        //disp_old_mf[lev]->setVal(Set::Vector::Zero());

        if (m_type == Type::Dynamic)
            if (velocity_ic)
            {
                velocity_ic->Initialize(lev,vel_mf);
                //velocity_ic->Initialize(lev,vel_old_mf);
                //ddw_mf[lev]->setVal(MATRIX4::Zero());
                //stress_mf[lev]->setVal(Set::Matrix::Zero());
                //strain_mf[lev]->setVal(Set::Matrix::Zero());
                //model_mf[lev]->setVal(MODEL::Zero());
            }

        if (ic_rhs) ic_rhs->Initialize(lev, rhs_mf);
        else rhs_mf[lev]->setVal(Set::Vector::Zero());
    }

    virtual void UpdateModel(int a_step, Set::Scalar a_time) = 0;

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::TimeStepBegin");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (m_zero_out_displacement) disp_mf[lev]->setVal(Set::Vector::Zero());
            rhs_mf[lev]->setVal(Set::Vector::Zero());
            if (ic_rhs) ic_rhs->Initialize(lev, rhs_mf);
        }

        UpdateModel(a_step, a_time);

        bc->SetTime(a_time);
        bc->Init(rhs_mf, geom);

        if (m_type != Mechanics<MODEL>::Type::Static) return;
        if (a_time < tstart) return;
        if (m_interval && a_step % m_interval) return;

        amrex::LPInfo info;
        if (m_max_coarsening_level >= 0)
            info.setMaxCoarseningLevel(m_max_coarsening_level);
        Operator::Elastic<MODEL::sym> elastic_op(Geom(0, finest_level), grids, DistributionMap(0, finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetHomogeneous(false);
        elastic_op.SetBC(bc);
        IO::ParmParse pp("elasticop");
        pp_queryclass(elastic_op);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        solver.Define(elastic_op);
        if (psi_on) solver.setPsi(psi_mf);
        solver.solve(disp_mf, rhs_mf, model_mf, tol_rel, tol_abs);
        if (m_print_residual) solver.compLinearSolverResidual(res_mf, disp_mf, rhs_mf);
        solver.Clear();

        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());

            const amrex::Real* DX = geom[lev].CellSize();
            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                bx.grow(2);
                bx = bx & domain;
                amrex::Array4<MODEL>             const& model = model_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& stress = stress_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& strain = strain_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector> const& disp = disp_mf[lev]->array(mfi);


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp, i, j, k, DX, sten);
                        stress(i, j, k) = model(i, j, k).DW(F);
                        strain(i, j, k) = F;
                    }
                    else
                    {
                        Set::Matrix gradu = Numeric::Gradient(disp, i, j, k, DX, sten);
                        stress(i, j, k) = model(i, j, k).DW(gradu);
                        strain(i, j, k) = 0.5 * (gradu + gradu.transpose());
                    }
                });
            }
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);
            Util::RealFillBoundary(*strain_mf[lev], geom[lev]);
        }
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Advance");
        if (m_type == Mechanics<MODEL>::Type::Disable) return;
        const amrex::Real* DX = geom[lev].CellSize();

        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

        if (m_type == Type::Dynamic)
        {

            std::swap(*disp_mf[lev], *disp_old_mf[lev]);
            std::swap(*vel_mf[lev], *vel_old_mf[lev]);


            for (amrex::MFIter mfi(*disp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Set::Patch<const Set::Vector> u     = disp_old_mf.Patch(lev,mfi);
                Set::Patch<const Set::Vector> v     = vel_old_mf.Patch(lev,mfi);
                Set::Patch<const Set::Vector> b     = rhs_mf.Patch(lev,mfi);

                Set::Patch<Set::Vector>       unew  = disp_mf.Patch(lev,mfi);
                Set::Patch<Set::Vector>       vnew  = vel_mf.Patch(lev,mfi);
                Set::Patch<Set::Matrix>       eps   = strain_mf.Patch(lev,mfi);
                Set::Patch<Set::Matrix>       sig   = stress_mf.Patch(lev,mfi);
                Set::Patch<MATRIX4>           ddw   = ddw_mf.Patch(lev,mfi);
                Set::Patch<MODEL>             model = model_mf.Patch(lev,mfi);

                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i,j,k,bx);
                    eps(i, j, k) = Numeric::Gradient(u, i, j, k, DX,sten);
                    sig(i, j, k) = model(i, j, k).DW(eps(i, j, k));
                    ddw(i,j,k) = model(i,j,k).DDW(eps(i,j,k));
                });

                bx = mfi.nodaltilebox() & domain;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {

                    bool AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z));
                    bool AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));

                    if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                    {
                        auto sten = Numeric::GetStencil(i,j,k,domain);                        

                        auto bctype = bc->getType(i,j,k,domain);

                        for (int d = 0; d < AMREX_SPACEDIM; d++)
                        {
                            if (bctype[d] == BC::Operator::Elastic::Elastic::Type::Displacement)
                            {
                                unew(i,j,k)(d) = b(i,j,k)(d);
                            }
                            else if (bctype[d] == BC::Operator::Elastic::Elastic::Type::Traction)
                            {
                                
                                Set::Vector N = Set::Normal(AMREX_D_DECL(xmin,ymin,zmin),
                                                            AMREX_D_DECL(xmax,ymax,zmax));
                                
                                auto [phi, offdiag] = Numeric::GradientSplit(u,i,j,k,DX,sten);
                                
                                Set::Matrix A = Set::Matrix::Zero();
                                Set::Vector rhs = b(i,j,k);
                                //Util::Message(INFO,b(i,j,k).transpose());

                                for (int p = 0; p < AMREX_SPACEDIM; p++)
                                    for (int q = 0; q < AMREX_SPACEDIM; q++)
                                        for (int r = 0; r < AMREX_SPACEDIM; r++)
                                            for (int s = 0; s < AMREX_SPACEDIM; s++)
                                            {
                                                A(p,r) += ddw(i,j,k)(p,q,r,s) * phi(s) * N(q);

                                                rhs(p) -= ddw(i,j,k)(p,q,r,s) * eps(i,j,k)(r,s) * N(q);
                                            }
                                Set::Vector delta_u = (A.inverse() * rhs);

                                unew(i,j,k)(d) = u(i,j,k)(d) + delta_u(d);
                                vnew(i,j,k)(d) = (unew(i,j,k)(d) - u(i,j,k)(d))/dt;
                            }
                            else
                            {
                                Util::Abort(INFO,"Elastic dynamics not supported for other BCs yet");
                            }
                        }
                    }
                    else
                    {
                        Set::Matrix      gradu = Numeric::Gradient(u,i,j,k,DX);
                        Set::Matrix3 gradgradu = Numeric::Hessian(u,i,j,k,DX);
                        
                        Set::Vector f = ddw(i,j,k)*gradgradu;

                        MATRIX4 AMREX_D_DECL(
                            Cgrad1 = (Numeric::Stencil<MATRIX4, 1, 0, 0>::D(ddw, i, j, k, 0, DX)),
                            Cgrad2 = (Numeric::Stencil<MATRIX4, 0, 1, 0>::D(ddw, i, j, k, 0, DX)),
                            Cgrad3 = (Numeric::Stencil<MATRIX4, 0, 0, 1>::D(ddw, i, j, k, 0, DX)));

                        f += AMREX_D_TERM(  ( Cgrad1*gradu).col(0),
                                            +(Cgrad2*gradu).col(1),
                                            +(Cgrad3*gradu).col(2));
                        
                        Set::Vector lapv = Numeric::Laplacian(v, i, j, k, DX);
                        f += mu_newton * lapv + b(i,j,k) - mu_dashpot*v(i,j,k);
                        
                        Set::Vector udotdot = f / rho;
                        vnew(i, j, k) = v(i, j, k) + dt * udotdot;
                        unew(i, j, k) = u(i, j, k) + dt * v(i, j, k);
                    }
                });
            }
        }

        for (amrex::MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            amrex::Array4<Set::Matrix>         const& eps = (*strain_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>         const& sig = (*stress_mf[lev]).array(mfi);
            amrex::Array4<MODEL>               const& model = (*model_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                model(i, j, k).Advance(dt, eps(i, j, k), sig(i, j, k));
            });
        }
        Util::RealFillBoundary(*model_mf[lev], geom[lev]);

    }

    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
        const amrex::MFIter& mfi, const amrex::Box& a_box) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::Integrate");
        if (m_type == Type::Disable) return;

        const amrex::Real* DX = geom[amrlev].CellSize();
        amrex::Box domain = geom[amrlev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

        amrex::Box box = a_box;
        box.convert(amrex::IntVect::TheNodeVector());


        //Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
#if AMREX_SPACEDIM == 2
        Set::Vector da0(DX[1], 0);
        Set::Vector da1(0, DX[0]);
#elif AMREX_SPACEDIM == 3
        Set::Vector da(DX[1] * DX[2], 0, 0);
#endif

        const Dim3 /*lo= amrex::lbound(domain),*/ hi = amrex::ubound(domain);
        const Dim3 /*boxlo= amrex::lbound(box),*/ boxhi = amrex::ubound(box);

        amrex::Array4<const Set::Matrix> const& stress = (*stress_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Vector> const& disp = (*disp_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
#if AMREX_SPACEDIM == 2
            if (i == hi.x && j < boxhi.y)
            {
                trac_hi[0] += (0.5 * (stress(i, j, k) + stress(i, j + 1, k)) * da0);
                disp_hi[0] = disp(i, j, k);
            }
            if (j == hi.y && i < boxhi.x)
            {
                trac_hi[1] += (0.5 * (stress(i, j, k) + stress(i + 1, j, k)) * da1);
                disp_hi[1] = disp(i, j, k);
            }
#elif AMREX_SPACEDIM == 3
            if (i == hi.x && (j < boxhi.y && k < boxhi.z))
            {
                trac_hi[0] += (0.25 * (stress(i, j, k) + stress(i, j + 1, k)
                    + stress(i, j, k + 1) + stress(i, j + 1, k + 1)) * da);
                disp_hi[0] = disp(i, j, k);
            }
#endif
        });

    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/) override
    {
        BL_PROFILE("Integrator::Base::Mechanics::TagCellsForRefinement");
        if (m_type == Type::Disable) return;

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*strain_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Matrix> const& eps = strain_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                Set::Matrix3 grad = Numeric::Gradient(eps, i, j, k, DX.data(), sten);
                if (grad.norm() * DXnorm > m_elastic_ref_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    typedef Set::Matrix4<AMREX_SPACEDIM,MODEL::sym> MATRIX4;
    Set::Field<MODEL> model_mf;
    Set::Field<MATRIX4> ddw_mf;
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
    Set::Scalar mu_dashpot = NAN;
    Set::Scalar mu_newton = NAN;

    //Set::Vector trac_lo[AMREX_SPACEDIM];
    Set::Vector trac_hi[AMREX_SPACEDIM];
    Set::Vector disp_hi[AMREX_SPACEDIM];


    IC::IC* ic_rhs = nullptr;
    BC::BC<Set::Scalar>* mybc;

    IC::IC* velocity_ic = nullptr;

    Solver::Nonlocal::Newton<MODEL> solver;//(elastic.op);
    BC::Operator::Elastic::Elastic* bc;

    Set::Scalar m_elastic_ref_threshold = 0.01;
    bool m_print_model = false;
    bool m_print_residual = false;
    bool m_time_evolving = false;
    int m_max_coarsening_level = -1;

    bool m_zero_out_displacement = false;

    bool plot_disp = true;
    bool plot_stress = true;
    bool plot_strain = true;
    bool plot_psi = true;
    bool plot_rhs = true;


    Set::Scalar tstart = -1;
};
}
}
#endif
