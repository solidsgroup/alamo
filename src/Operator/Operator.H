#ifndef OPERATOR_H
#define OPERATOR_H

#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BaseFab.H>
#include "AMReX_BaseFab_f.H"



#include "BC/BC.H"

using namespace amrex;

namespace Operator {

///
/// General abstract class for node-based operators
///
class Operator
	: public amrex::MLNodeLinOp
{


public :

	enum struct CoarseningStrategy { Sigma, RAP };

	Operator () {}
	Operator (const Vector<Geometry>& a_geom,
			 const Vector<BoxArray>& a_grids,
			 const Vector<DistributionMapping>& a_dmap,
			 const LPInfo& a_info = LPInfo(),
			 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
	virtual ~Operator ();

	Operator (const Operator&) = delete;
	Operator (Operator&&) = delete;
	Operator& operator= (const Operator&) = delete;
	Operator& operator= (Operator&&) = delete;

	void define (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info = LPInfo(),
		     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

	void setRZCorrection (bool rz) { m_is_rz = rz; }

	void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
		      const Vector<const MultiFab*>& rhnd,
		      const Vector<MultiFab*>& rhcc);

	void compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& phi,
				     const MultiFab& vold, const MultiFab* rhcc,
				     const BoxArray& fine_grids, const IntVect& ref_ratio);

	void compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
				   const MultiFab* rhcc);

	void setGaussSeidel (bool flag) { m_use_gauss_seidel = flag; }
	void setHarmonicAverage (bool flag) { m_use_harmonic_average = flag; }

	void setCoarseningStrategy (CoarseningStrategy cs) { m_coarsening_strategy = cs; }

	virtual void FApply (int amrlev, int mglev, MultiFab& out, const MultiFab& in)
	{Fapply(amrlev, mglev, out, in);};
	virtual void FSmooth (int amrlev, int mglev, MultiFab& out, const MultiFab& in)
	{Fsmooth(amrlev,mglev,out,in);};
	virtual void Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const;

	virtual void PrepareForSolve () {prepareForSolve();};

	bool VerificationCheck (int amrlev, int mglev, amrex::MultiFab& test) const;

	int NanInfOnBoundary(const amrex::MultiFab &test) const
	{
		int ncomp = test.nComp(), ngrow = test.nGrow();
		if (ngrow == 0) return 0;
	  	for (MFIter mfi(test, true); mfi.isValid(); ++mfi)
	  	{
	  		const Box& bx = mfi.tilebox();
	  		const amrex::FArrayBox &testfab    = test[mfi];
	  		for (int n = 0; n < ncomp; n++)
	  		{
	  			AMREX_D_TERM(for (int m1 = bx.loVect()[0]-ngrow; m1<=bx.hiVect()[0]+ngrow; m1++),
	  				     for (int m2 = bx.loVect()[1]-ngrow; m2<=bx.hiVect()[1]+ngrow; m2++),
	  				     for (int m3 = bx.loVect()[2]-ngrow; m3<=bx.hiVect()[2]+ngrow; m3++))
	  			{
					AMREX_D_TERM(if (m1 >= bx.loVect()[0] || m1 <= bx.hiVect()[0]) continue;,
						     if (m2 >= bx.loVect()[1] || m2 <= bx.hiVect()[1]) continue;,
						     if (m3 >= bx.loVect()[2] || m3 <= bx.hiVect()[2]) continue;)
						amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
	  				if (std::isnan(testfab(m,n)) || std::isnan(testfab(m,n))) Util::Abort(INFO,"nan detected: m=",m," n=",n);
	  			}
	  		}
	  	}
		return 0;
	}


protected:

	virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
	virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const override final;
	virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
					     const MultiFab& fine_sol, const MultiFab& fine_rhs) final;


	virtual void prepareForSolve () final;
	virtual bool isSingular (int amrlev) const final
        {return (amrlev == 0) ? m_is_bottom_singular : false; }
	virtual bool isBottomSingular () const final { return m_is_bottom_singular; }
	virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode,
			      amrex::MLLinOp::StateMode /**/, bool skip_fillboundary=false) const final;
	///
	/// Apply 
	/// 
	virtual void Diagonal (bool recompute=false);
	virtual void Diagonal (int amrlev, int mglev, amrex::MultiFab& diag);

	virtual void Fapply (int amrlev, ///<[in] AMR level
			     int mglev, ///<[in] MG level
			     MultiFab& out, ///<[out] Output of linear operator
			     const MultiFab& in ///<[in] Input of linear operator
			     ) const =0;
	virtual void Fsmooth (int amrlev,
			      int mglev,
			      MultiFab& x,
			      const MultiFab& b
			      ) const;
	virtual void normalize (int amrlev,
				int mglev,
				MultiFab& mf
				) const;
	virtual void reflux (int crse_amrlev,
			     MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const =0;

	virtual Real getAScalar () const final { return  0.0; }
	virtual Real getBScalar () const final { return -1.0; }

	virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

	std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {{AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
								      amrex::IntVect(AMREX_D_DECL(0,1,0)),
								      amrex::IntVect(AMREX_D_DECL(0,0,1)))}};

private:

	int m_is_rz = 0;

	CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::Sigma;

	bool m_use_gauss_seidel = true;
	bool m_use_harmonic_average = false;

	bool m_is_bottom_singular = false;
	bool m_masks_built = false;
	//
	// functions
	//

	void buildMasks ();
public:
	void RegisterNewFab(amrex::Vector<amrex::MultiFab> &input);
	void RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input);
	const amrex::FArrayBox & GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const;

	virtual void averageDownCoeffs () = 0;

protected:
	int m_num_a_fabs = 0;
	bool m_diagonal_computed = false;
	amrex::Vector<amrex::Vector<amrex::Vector<amrex::MultiFab> > > m_a_coeffs;
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > m_diag;
};
}

#endif
