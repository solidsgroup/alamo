#ifndef SOLVER_NONLOCAL_NEWTONRAPHSON
#define SOLVER_NONLOCAL_NEWTONRAPHSON

#include "Set/Set.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"

namespace Solver
{
namespace Nonlocal
{
template <class T>
class NewtonRaphson : public Linear
{
public:
    NewtonRaphson(Operator::Elastic<T> & a_op) : Linear(a_op) {} ;
    
    void setNRIters(int a_nriters) {m_nriters = a_nriters;}
    Set::Scalar solve( amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_sol,                      
                       amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_rhs,
                       amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > & a_model,
                       Set::Scalar a_tol_rel, Set::Scalar a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Operator::Elastic<T> &elastic = dynamic_cast<Operator::Elastic<T> &>(linop);
        
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > a_dsol;
        a_dsol.resize(a_sol.size());
        for (int i = 0; i < a_sol.size(); i++) 
        {
            a_dsol[i].reset(new amrex::MultiFab(a_sol[i]->boxArray(), 
                                             a_sol[i]->DistributionMap(),
                                             a_sol[i]->nComp(), a_sol[i]->nGrow()));
        }

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            Util::Message(INFO,"Newton Raphson Iteration ", nriter);
            for (int lev = 0; lev < a_rhs.size(); ++lev)
            {
                a_dsol[lev]->setVal(0.0);
                
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar * dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());

                for (MFIter mfi(a_model[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.growntilebox(2);

                    amrex::Array4<Set::Scalar> const &u     = a_sol[lev]->array(mfi);
                    amrex::Array4<T>           const &model = a_model[lev].array(mfi);

                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                    {
    	                std::array<Numeric::StencilType,AMREX_SPACEDIM> sten = Numeric::GetStencil(i,j,k,bx);

                        Set::Matrix gradu = Numeric::Gradient(u,i,j,k,dx,sten);
                        Set::Matrix F = gradu + Set::Matrix::Identity();

                        model(i,j,k).dw  = model(i,j,k).DW(F);
                        model(i,j,k).ddw = model(i,j,k).DDW(F);
                    });

                    bx = mfi.growntilebox(1);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int p = 0; p < AMREX_SPACEDIM; p++)
                            model(i,j,k).divdw(p) = 
                                (model(i+1,j,k).dw(p,0) - model(i-1,j,k).dw(p,0))/2./DX(0) + 
                                (model(i,j+1,k).dw(p,1) - model(i,j-1,k).dw(p,1))/2./DX(1) + 
                                (model(i,j,k+1).dw(p,2) - model(i,j,k-1).dw(p,2))/2./DX(2);

                    });
                }


                amrex::Geometry geom = elastic.Geom(lev);
                for (int i = 0; i < 2; i++)
                {
                    amrex::FabArray<amrex::BaseFab<T>> &mf = a_model[lev];
                    mf.FillBoundary(geom.periodicity());
                    const int ncomp = mf.nComp();
                    const int ng1 = 1;
                    const int ng2 = 2;
                    amrex::FabArray<amrex::BaseFab<T>> tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
                    amrex::Copy(tmpmf, mf, 0, 0, ncomp, ng1);
                    mf.ParallelCopy(tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
                }                

            }
            elastic.SetModel(a_model);
            solveaffine(a_dsol, a_rhs, a_tol_rel, a_tol_abs,false);
            for (int lev = 0; lev < a_dsol.size(); ++lev)
                amrex::MultiFab::Add(*a_sol[lev],*a_dsol[lev],0,0,AMREX_SPACEDIM,2);


            for (int lev = 0; lev < a_rhs.size(); ++lev)
                for (MFIter mfi(a_model[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.tilebox();
    
                    amrex::Array4<Set::Scalar> const &rhs   = a_rhs[lev]->array(mfi);
    
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                    {
                        if (i == bx.loVect()[0] || j == bx.loVect()[1] || k == bx.loVect()[2] ||
                            i == bx.hiVect()[0] || j == bx.hiVect()[1] || k == bx.hiVect()[2] )
                        {
                            rhs(i,j,k,0) = 0.0;
                            rhs(i,j,k,1) = 0.0;
                            rhs(i,j,k,2) = 0.0;
                        }
                   });
                }
            
        }

    }
private:
    int m_nriters = 1;
};
}
}



#endif