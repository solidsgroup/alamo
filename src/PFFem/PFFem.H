
#ifndef _PFFEM_H_
#define _PFFEM_H_

#include <AMReX_LinOp.H>

namespace amrex {


class PFFem
    :
    public LinOp
{
public:
    //
    // constructor;
    //
    PFFem (const BndryData& bd,
               Real             _h);
    //
    // destructor
    //
    virtual ~PFFem() override;
    //
    // Compute flux associated with the op
    //
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
			   int sComp=0, int dComp=0, int nComp=1, int bndComp=0) override;
    
    virtual Real norm (int nm = 0, int level = 0, const bool local = false) override;

protected:
    //
    // compute out=L(in) at level=level
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) override;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) override;
    // //
    // // apply GSRB smoother to improve residual to L(solnL)=rhsL
    // //
    // virtual void Fsmooth (MultiFab&       solnL,
    //                       const MultiFab& rhsL,
    //                       int             level,
    //                       int             rgbflag) override;
    //
    // apply Jacobi smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) override;

private:
  amrex::Real mu, lambda; // Lame parameters 

  inline amrex::Real Phi(int p1, int p2, int n1, int n2, int i, int j);
};

}

#endif /*_PFFEM_H_*/
