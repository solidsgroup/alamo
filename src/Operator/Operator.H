#ifndef OPERATOR_H
#define OPERATOR_H

#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BaseFab.H>
#include "AMReX_BaseFab_f.H"

#include "BC/BC.H"

#include "Test/Operator/Elastic.H"

#ifdef __cplusplus
extern "C"
{
#endif
	// Defined in AMReX
	void amrex_mg_average (amrex_real* crse,       AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
			       const amrex_real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
			       const int *tlo, const int *thi,
			       const int *nc);

	// Defined in AMReX
	void amrex_mg_interp (amrex_real* fine,       AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
			      const amrex_real* crse, AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
			      const int *tlo, const int *thi,
			      const int *nc);
#ifdef __cplusplus
}
#endif


using namespace amrex;

enum Grid { Cell, Node};

namespace Operator {

template <Grid G> class Operator;

//
//
//  NODE-BASED OPERATOR
//
//

template <>
class Operator<Grid::Node> : public amrex::MLNodeLinOp
{
	friend Test::Operator::Elastic;
public :

	enum struct CoarseningStrategy { Sigma, RAP };

	Operator () {}
	Operator (const amrex::Vector<amrex::Geometry>& a_geom,
		  const amrex::Vector<amrex::BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  const LPInfo& a_info = LPInfo(),
		  const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
	virtual ~Operator ();

	virtual int getNGrow() const override {return 2;}

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

	virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
				       const MultiFab* crse_bcdata=nullptr) override;
	virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
					 BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override;

	void compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& phi,
				     const MultiFab& vold, const MultiFab* rhcc,
				     const BoxArray& fine_grids, const IntVect& ref_ratio);

	void compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
				   const MultiFab* rhcc);

	void setGaussSeidel (bool flag) { m_use_gauss_seidel = flag; }
	void setHarmonicAverage (bool flag) { m_use_harmonic_average = flag; }

	void setCoarseningStrategy (CoarseningStrategy cs) { m_coarsening_strategy = cs; }

	virtual void Reflux (int crse_amrlev,
			     MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
	{reflux(crse_amrlev,res,crse_sol,crse_rhs,fine_res,fine_sol,fine_rhs);}
	virtual void FApply (int amrlev, int mglev, MultiFab& out, const MultiFab& in)
	{Fapply(amrlev, mglev, out, in);};
	virtual void FSmooth (int amrlev, int mglev, MultiFab& out, const MultiFab& in)
	{Fsmooth(amrlev,mglev,out,in);};
	virtual void Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const = 0;

	virtual void PrepareForSolve () {prepareForSolve();};

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

	void BuildMasks() {buildMasks();}
protected:

	void realFillBoundary(MultiFab &phi, const Geometry &geom) const;

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

	virtual void Fapply (int amrlev, 
			     int mglev, 
			     MultiFab& out,
			     const MultiFab& in 
			     ) const override =0;
	virtual void Fsmooth (int amrlev,
			      int mglev,
			      MultiFab& x,
			      const MultiFab& b
			      ) const override;
	virtual void normalize (int amrlev,
				int mglev,
				MultiFab& mf
				) const override;
	virtual void reflux (int crse_amrlev,
			     MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const override;

	virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

	
	static constexpr amrex::IntVect AMREX_D_DECL(dx = {AMREX_D_DECL(1,0,0)},
						     dy = {AMREX_D_DECL(0,1,0)},
						     dz = {AMREX_D_DECL(0,0,1)});

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



//
//
//  CELL-BASED OPERATOR
//
//


template<>
class Operator<Grid::Cell> : public amrex::MLCellLinOp
{
public:

	friend class MLMG;
	friend class MLCGSolver;

	Operator ();
	virtual ~Operator () {};

	Operator (const Operator&) = delete;
	Operator (Operator&&) = delete;
	Operator& operator= (const Operator&) = delete;
	Operator& operator= (Operator&&) = delete;

	void define (amrex::Vector<amrex::Geometry> a_geom,
		     const amrex::Vector<amrex::BoxArray>& a_grids,
		     const amrex::Vector<amrex::DistributionMapping>& a_dmap,
		     BC::BC& a_bc,
		     const amrex::LPInfo& a_info = amrex::LPInfo(),
		     const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory = {});

	// virtual void setLevelBC (int amrlev, const amrex::MultiFab* levelbcdata) override final;

protected:
	
	BC::BC *m_bc;

	amrex::Vector<std::unique_ptr<amrex::MLMGBndry> >   m_bndry_sol;
	amrex::Vector<std::unique_ptr<amrex::BndryRegister> > m_crse_sol_br;

	amrex::Vector<std::unique_ptr<amrex::MLMGBndry> > m_bndry_cor;
	amrex::Vector<std::unique_ptr<amrex::BndryRegister> > m_crse_cor_br;

	// In case of agglomeration, coarse MG grids on amr level 0 are
	// not simply coarsened from fine MG grids.  So we need to build
	// bcond and bcloc for each MG level.
	using RealTuple = std::array<amrex::Real,2*BL_SPACEDIM>;
	using BCTuple   = std::array<amrex::BoundCond,2*BL_SPACEDIM>;
	class BndryCondLoc
	{
	public:
		BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
		void setLOBndryConds (const amrex::Geometry& geom, const amrex::Real* dx,
				      const std::array<BCType,AMREX_SPACEDIM>& lobc,
				      const std::array<BCType,AMREX_SPACEDIM>& hibc,
				      int ratio, const amrex::RealVect& a_loc);
		const BCTuple& bndryConds (const amrex::MFIter& mfi) const {
			return bcond[mfi];
		}
		const RealTuple& bndryLocs (const amrex::MFIter& mfi) const {
			return bcloc[mfi];
		}
	private:
		amrex::LayoutData<BCTuple>   bcond;
		amrex::LayoutData<RealTuple> bcloc;
	};
	amrex::Vector<amrex::Vector<std::unique_ptr<BndryCondLoc> > > m_bcondloc;

	// used to save interpolation coefficients of the first interior cells
	mutable amrex::Vector<amrex::Vector<amrex::BndryRegister> > m_undrrelxr;

	// boundary cell flags for covered, not_covered, outside_domain
	amrex::Vector<amrex::Vector<std::array<amrex::MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;



	//
	// functions
	//

	void updateSolBC (int amrlev, const amrex::MultiFab& crse_bcdata) const;
	void updateCorBC (int amrlev, const amrex::MultiFab& crse_bcdata) const;

	// virtual void applyBC (int amrlev, int mglev, amrex::MultiFab& in, BCMode bc_mode,
	// 		      const amrex::MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const;

	// virtual void restriction (int, int,
	// 			  amrex::MultiFab& crse,
	// 			  amrex::MultiFab& fine) const override final;

	// virtual void interpolation (int amrlev, int fmglev,
	// 			    amrex::MultiFab& fine,
	// 			    const amrex::MultiFab& crse) const override final;

	// virtual void averageDownSolutionRHS (int camrlev,
	// 				     amrex::MultiFab& crse_sol,
	// 				     amrex::MultiFab& crse_rhs,
	// 				     const amrex::MultiFab& fine_sol,
	// 				     const amrex::MultiFab& fine_rhs) override final;

	// virtual void apply (int amrlev, int mglev,
	// 		    amrex::MultiFab& out,
	// 		    amrex::MultiFab& in,
	// 		    BCMode bc_mode,
	// 		    StateMode s_mode,
	// 		    const amrex::MLMGBndry* bndry=nullptr) const override final;


	// virtual void compGrad (int amrlev,
	// 		       const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& grad,
	// 		       amrex::MultiFab& sol, amrex::MLLinOp::Location loc) const override final;

	virtual void prepareForSolve () override;


	// PURE VIRTUAL METHODS

	virtual void Fapply (int amrlev, int mglev, amrex::MultiFab& out, const amrex::MultiFab& in) const = 0;
	virtual void Fsmooth (int amrlev, int mglev, amrex::MultiFab& sol, const amrex::MultiFab& rsh, int redblack) const = 0;
	virtual void FFlux (int amrlev, const MFIter& mfi,
			    const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
			    const FArrayBox& sol, Location loc, const int face_only=0) const = 0;

	void RegisterNewFab(amrex::Vector<amrex::MultiFab> &input);
	void RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input);
	const amrex::FArrayBox & GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) ;

	virtual bool isSingular (int /*amrlev*/) const final { return false; }
	virtual bool isBottomSingular () const final { return false; }

	virtual amrex::Real getAScalar () const final { return 0.0; }
	virtual amrex::Real getBScalar () const final { return 0.0; }

	virtual amrex::MultiFab const* getACoeffs (int /*amrlev*/, int /*mglev*/) const final { return nullptr;}
	virtual std::array<amrex::MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final {
		std::array<amrex::MultiFab const*,AMREX_SPACEDIM> ret;
		AMREX_D_TERM(ret[0] = nullptr;, ret[1] = nullptr;,ret[2] = nullptr;);
		return ret;}

	virtual std::unique_ptr<amrex::MLLinOp> makeNLinOp (int /*grid_size*/) const final {
		Util::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
		return std::unique_ptr<MLLinOp>{};
	}

	void averageDownCoeffs ();
	void averageDownCoeffsSameAmrLevel (amrex::Vector<amrex::MultiFab>&);
	const amrex::FArrayBox & GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const;

	static constexpr amrex::IntVect AMREX_D_DECL(dx = {AMREX_D_DECL(1,0,0)},
						     dy = {AMREX_D_DECL(0,1,0)},
						     dz = {AMREX_D_DECL(0,0,1)});
private:

	// void defineAuxData ();
	// void defineBC ();

	int m_num_a_fabs = 0;
	// m_a_coeffs[index of fab][amrlevel][mglevel]
	amrex::Vector<amrex::Vector<amrex::Vector<amrex::MultiFab> > > m_a_coeffs; 

};



}


#endif
