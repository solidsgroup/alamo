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

using namespace amrex;

enum Grid { Cell, Node};

/// \brief Documentation for operator namespace
namespace Operator
{

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

	//
	// Public: Constructor/Destructor/Operators/defines
	//
public :
	Operator () {}
	Operator (const amrex::Vector<amrex::Geometry>& a_geom,
		  const amrex::Vector<amrex::BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  const LPInfo& a_info = LPInfo(),
		  const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
	virtual ~Operator ();
	Operator (const Operator&) = delete;
	Operator (Operator&&) = delete;
	Operator& operator= (const Operator&) = delete;
	Operator& operator= (Operator&&) = delete;
	void define (const Vector<Geometry>& a_geom, const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info = LPInfo(),
		     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
	const Geometry& Geom (int amr_lev, int mglev=0) const noexcept { return m_geom[amr_lev][mglev]; }

	//
	// Public Utilty functions
	//
public:
	void RegisterNewFab(amrex::Vector<amrex::MultiFab> &input);
	void RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input);
	const amrex::FArrayBox & GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const;
	virtual void SetHomogeneous (bool) {};
	//
	// Pure Virtual: you MUST override these functions
	//
protected:
	virtual void Fapply (int amrlev, int mglev,MultiFab& out,const MultiFab& in) const override =0;
	virtual void averageDownCoeffs () = 0;
	//
	// Virtual: you SHOULD override these functions
	//
protected:
	virtual void Diagonal (bool recompute=false);
	virtual void Diagonal (int amrlev, int mglev, amrex::MultiFab& diag);
	//
	// Virtual: you CAN override these functions (but probably don't need to)
	//
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& x,const MultiFab& b) const override;
	virtual void normalize (int amrlev, int mglev, MultiFab& mf) const override;
	virtual void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const override;
	//
	// Virtual: you CAN'T override these functions (they take care of AMReX business) and you SHOULDN'T
	//
protected:
	virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
	virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const override final;
	virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs, const MultiFab& fine_sol, const MultiFab& fine_rhs) final;
	virtual void prepareForSolve () override;
	virtual bool isSingular (int amrlev) const final {return (amrlev == 0) ? m_is_bottom_singular : false; }
	virtual bool isBottomSingular () const final { return m_is_bottom_singular; }
	virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, amrex::MLLinOp::StateMode /**/, bool skip_fillboundary=false) const final;
	virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;
public:
	virtual int getNGrow() const override final {return 2;}
	virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
				       const MultiFab* crse_bcdata=nullptr) override final;
	virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
					 BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override final;
	//
	// Private: internal machinery that you should never need to mess with
	//
private:
	void buildMasks ();
public:
	static void realFillBoundary(MultiFab &phi, const Geometry &geom);

private:
	bool m_is_bottom_singular = false;
	bool m_masks_built = false;
	/// \todo we need to get rid of this
	// static constexpr amrex::IntVect AMREX_D_DECL(dx = {AMREX_D_DECL(1,0,0)},
	//    					     dy = {AMREX_D_DECL(0,1,0)},
	//    					     dz = {AMREX_D_DECL(0,0,1)});
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
				      const amrex::Array<BCType,AMREX_SPACEDIM>& lobc,
				      const amrex::Array<BCType,AMREX_SPACEDIM>& hibc,
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

	virtual void prepareForSolve () override;


	// PURE VIRTUAL METHODS

	virtual void Fapply (int amrlev, int mglev, amrex::MultiFab& out, const amrex::MultiFab& in) const override = 0;
	virtual void Fsmooth (int amrlev, int mglev, amrex::MultiFab& sol, const amrex::MultiFab& rsh, int redblack) const override = 0;
	virtual void FFlux (int amrlev, const MFIter& mfi,
			    const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
			    const FArrayBox& sol, Location loc, const int face_only=0) const override = 0;

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

	int m_num_a_fabs = 0;
	amrex::Vector<amrex::Vector<amrex::Vector<amrex::MultiFab> > > m_a_coeffs; 

};



}


#endif
