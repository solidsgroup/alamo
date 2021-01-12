#ifndef SOLVER_NONLOCAL_NEWTON
#define SOLVER_NONLOCAL_NEWTON

#include "Set/Set.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Elastic/NeoHookean.H"
#include "Numeric/Stencil.H"

namespace Solver
{
namespace Nonlocal
{
template <typename T>
class Newton : public Linear
{
public:
    Newton(Operator::Elastic<T::sym> &a_op) : 
        Linear(a_op), 
        m_elastic(dynamic_cast<Operator::Elastic<T::sym> &>(linop)),
        m_bc(m_elastic.GetBC())
    {};

    void setNRIters(int a_nriters) { m_nriters = a_nriters; }


private:

    void prepareForSolve(const Set::Field<Set::Scalar>& a_u_mf, 
                         const Set::Field<Set::Scalar>& a_b_mf,
                         Set::Field<Set::Scalar>& a_rhs_mf,
                         Set::Field<Set::Matrix> &a_dw_mf,
                         Set::Field<Set::Matrix4<AMREX_SPACEDIM,T::sym>> &a_ddw_mf,
                         Set::Field<T> &a_model_mf)
    {
            for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
            {
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar *dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());
                
                for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox();
                    bx = bx & domain;

                    amrex::Array4<const T>           const &model = a_model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &u     = a_u_mf[lev]->array(mfi);
                    amrex::Array4<Set::Matrix>       const &dw    = a_dw_mf[lev]->array(mfi);
                    amrex::Array4<Set::Matrix4<AMREX_SPACEDIM,T::sym>>  const &ddw = a_ddw_mf[lev]->array(mfi);

                    // Set model internal dw and ddw.
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> sten = Numeric::GetStencil(i, j, k, bx);

                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);

                        if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::gradu)
                        {
                            dw(i,j,k) = model(i, j, k).DW(gradu);
                            ddw(i,j,k) = model(i, j, k).DDW(gradu);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        {
                            Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                            dw(i,j,k) = model(i, j, k).DW(eps);
                            ddw(i,j,k) = model(i, j, k).DDW(eps);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                        {
                            Set::Matrix F = gradu + Set::Matrix::Identity();
                            dw(i,j,k) = model(i, j, k).DW(F);
                            ddw(i,j,k) = model(i, j, k).DDW(F);
                        }
                    });
                }

                Util::RealFillBoundary(*a_dw_mf[lev],m_elastic.Geom(lev));
                Util::RealFillBoundary(*a_ddw_mf[lev],m_elastic.Geom(lev));
            }

            m_elastic.SetModel(a_ddw_mf);

            for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
            {
                //amrex::MultiFab::Copy(a_rhs_mf[lev],a_b_mf[lev],0,0,a_rhs_mf[lev]->nComp(),a_rhs_mf[lev]->nGrow());
                m_elastic.Apply(lev,0,*a_rhs_mf[lev],*a_u_mf[lev]);
                a_rhs_mf[lev]->mult(-1.0,2);
                amrex::MultiFab::Add(*a_rhs_mf[lev],*a_b_mf[lev],0,0,a_rhs_mf[lev]->nComp(),2);

                
                Util::RealFillBoundary(*a_rhs_mf[lev],m_elastic.Geom(lev));
            }
            for (int lev = a_b_mf.finest_level; lev > 0; lev--)
            {
                linop.Reflux(lev-1, *a_rhs_mf[lev-1], *a_u_mf[lev-1], *a_b_mf[lev-1], *a_rhs_mf[lev], *a_u_mf[lev], *a_b_mf[lev]);
            }

    }

public:
    Set::Scalar solve (const Set::Field<Set::Scalar> & a_u_mf, 
                       const Set::Field<Set::Scalar> & a_b_mf,
                       Set::Field<T> &a_model_mf,
                       Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Set::Field<Set::Scalar> dsol_mf, rhs_mf;
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM,T::sym>> ddw_mf;

        dsol_mf.resize(a_u_mf.finest_level+1); dsol_mf.finest_level = a_u_mf.finest_level;
        dw_mf.  resize(a_u_mf.finest_level+1); dw_mf  .finest_level = a_u_mf.finest_level;
        ddw_mf. resize(a_u_mf.finest_level+1); ddw_mf .finest_level = a_u_mf.finest_level;
        rhs_mf. resize(a_u_mf.finest_level+1); rhs_mf .finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                                a_u_mf[lev]->DistributionMap(),
                                a_u_mf[lev]->nComp(), 
                                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev,   a_b_mf[lev]->boxArray(),
                                a_b_mf[lev]->DistributionMap(),
                                1, 
                                a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev,  a_b_mf[lev]->boxArray(),
                                a_b_mf[lev]->DistributionMap(),
                                1, 
                                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev,  a_b_mf[lev]->boxArray(),
                                a_b_mf[lev]->DistributionMap(),
                                a_b_mf[lev]->nComp(), 
                                a_b_mf[lev]->nGrow());
            
            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            ddw_mf[lev]->setVal(Set::Matrix4<AMREX_SPACEDIM,T::sym>::Zero());
            
            amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (m_verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter+1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);


            if (nriter == m_nriters) break;
            
            Solver::Nonlocal::Linear::solve(GetVecOfPtrs(dsol_mf), GetVecOfConstPtrs(rhs_mf), a_tol_rel, a_tol_abs,checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp,0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    Set::Scalar tmpsolnorm = a_u_mf[lev]->norm0(comp,0);
                    if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (m_verbose > 1) Util::Message(INFO,"NR relative norm(ddisp) = ",relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        return 0.0;
    }
    Set::Scalar solve (const Set::Field<Set::Scalar> & a_u_mf, 
                       const Set::Field<Set::Scalar> & a_b_mf,
                       Set::Field<T> &a_model_mf)
    {
        return solve(a_u_mf,a_b_mf,a_model_mf,m_tol_rel,m_tol_abs);
    }

    void compResidual(Set::Field<Set::Scalar> & a_res_mf,
                      Set::Field<Set::Scalar> & a_u_mf,
                      Set::Field<Set::Scalar> & a_b_mf,
                      Set::Field<T> &a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM,T::sym>> ddw_mf;
        dw_mf.resize(a_u_mf.size());
        ddw_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                              a_b_mf[lev]->DistributionMap(),
                              1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                               a_b_mf[lev]->DistributionMap(),
                               1, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }
        
        //for (int lev = 0; lev < a_b_mf.size(); ++lev)
        //m_elastic.GetBC().Init(a_b_mf[lev].get(),m_elastic.Geom(lev),true);
            
        prepareForSolve(a_u_mf, a_b_mf, a_res_mf, dw_mf, ddw_mf, a_model_mf);
        
    }

    void W(Set::Field<Set::Scalar> & a_w_mf,
           Set::Field<Set::Scalar> & a_u_mf,
           Set::Field<T> &a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
        	BL_PROFILE("Solver::Nonlocal::Newton::DW()");

        	const amrex::Real* DX = linop.Geom(lev).CellSize();
        	amrex::Box domain(linop.Geom(lev).Domain());
        	domain.convert(amrex::IntVect::TheNodeVector());

        	for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        	{
        		const Box& bx = mfi.tilebox();
        		amrex::Array4<T> const& C                 = a_model_mf[lev]->array(mfi);
        		amrex::Array4<amrex::Real> const& w       = a_w_mf[lev]->array(mfi);
        		amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
        		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        				    {
        					    Set::Matrix gradu;

        					    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
        						    = Numeric::GetStencil(i,j,k,domain);

        					    // Fill gradu
        					    for (int p = 0; p < AMREX_SPACEDIM; p++)
        					    {
        						    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
        						     		     gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
        						      		     gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
        					    }

                                if      (C(i,j,k).kinvar == Model::Solid::KinematicVariable::gradu) 
                                    w(i,j,k) = C(i,j,k).W(gradu);
                                else if (C(i,j,k).kinvar == Model::Solid::KinematicVariable::epsilon) 
                                    w(i,j,k) = C(i,j,k).W(0.5*(gradu + gradu.transpose()));
                                else if (C(i,j,k).kinvar == Model::Solid::KinematicVariable::F) 
                                    w(i,j,k) = C(i,j,k).W(gradu + Set::Matrix::Identity());
        				    });
        	}
        }
    }

    void DW(Set::Field<Set::Scalar> & a_dw_mf,
            Set::Field<Set::Scalar> & a_u_mf,
            Set::Field<T> &a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
        	BL_PROFILE("Solver::Nonlocal::Newton::DW()");

        	const amrex::Real* DX = linop.Geom(lev).CellSize();
        	amrex::Box domain(linop.Geom(lev).Domain());
        	domain.convert(amrex::IntVect::TheNodeVector());

        	for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        	{
        		const Box& bx = mfi.tilebox();
        		amrex::Array4<T> const& C                 = a_model_mf[lev]->array(mfi);
        		amrex::Array4<amrex::Real> const& dw      = a_dw_mf[lev]->array(mfi);
        		amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
        		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        				    {
        					    Set::Matrix gradu;

        					    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
        						    = Numeric::GetStencil(i,j,k,domain);

        					    // Fill gradu
        					    for (int p = 0; p < AMREX_SPACEDIM; p++)
        					    {
        						    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
        						     		     gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
        						      		     gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
        					    }

                                Set::Matrix sig = Set::Matrix::Zero();

                                if      (C(i,j,k).kinvar == Model::Solid::KinematicVariable::gradu) 
                                    sig = C(i,j,k).DW(gradu);
                                else if (C(i,j,k).kinvar == Model::Solid::KinematicVariable::epsilon) 
                                    sig = C(i,j,k).DW(0.5*(gradu + gradu.transpose()));
                                else if (C(i,j,k).kinvar == Model::Solid::KinematicVariable::F) 
                                    sig = C(i,j,k).DW(gradu + Set::Matrix::Identity());

        					    // = C(i,j,k)(gradu,m_homogeneous);

        						AMREX_D_PICK(dw(i,j,k,0) = sig(0,0);
        						             ,
        						             dw(i,j,k,0) = sig(0,0); dw(i,j,k,1) = sig(0,1); 
        						             dw(i,j,k,2) = sig(1,0); dw(i,j,k,3) = sig(1,1);
        						             ,
        						             dw(i,j,k,0) = sig(0,0); dw(i,j,k,1) = sig(0,1); dw(i,j,k,2) = sig(0,2); 
        						             dw(i,j,k,3) = sig(1,0); dw(i,j,k,4) = sig(1,1); dw(i,j,k,5) = sig(1,2); 
        						             dw(i,j,k,6) = sig(2,0); dw(i,j,k,7) = sig(2,1); dw(i,j,k,8) = sig(2,2););

        				    });
        	}
        }
    }


private:
    int m_nriters = 1;
    Operator::Elastic<T::sym> &m_elastic;
    BC::Operator::Elastic::Elastic &m_bc;
public:
    static void Parse(Newton<T> & value, amrex::ParmParse & pp)
    {
        Linear::Parse(value,pp);
        
        pp.query("nriters",value.m_nriters);
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif
