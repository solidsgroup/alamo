#ifndef OPERATOR_ELASTIC_H_
#define OPERATOR_ELASTIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>
#include "Set/Set.H"
#include "Operator/Operator.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/LinearElastic/Isotropic.H"

using namespace amrex;

namespace Operator
{
/// 
/// \todo Remove Elastic derived classes. They have been replaced with the Model construction.
/// 
template<class T>
class Elastic : public Operator
{
	using TArrayBox = amrex::BaseFab<T>;
	using MultiTab  = amrex::FabArray<TArrayBox>;
public:
	enum class BC {Displacement, Traction, Periodic}; 

	enum class Boundary {Lo, Hi, None};

	Elastic () {}
	Elastic (const Vector<Geometry>& a_geom,
		 const Vector<BoxArray>& a_grids,
		 const Vector<DistributionMapping>& a_dmap,
		 const LPInfo& a_info);
	virtual ~Elastic ();
	Elastic (const Elastic&) = delete;
	Elastic (Elastic&&) = delete;
	Elastic& operator= (const Elastic&) = delete;
	Elastic& operator= (Elastic&&) = delete;

	void define (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info = LPInfo(),
		     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});


	void SetModel (int amrlev, const amrex::FabArray<amrex::BaseFab<T> >& a_model);
	void Strain (int amrlev, amrex::MultiFab& epsfab, const amrex::MultiFab& ufab, bool voigt = false) const;
	void Stress (int amrlev, amrex::MultiFab& sigmafab, const amrex::MultiFab& ufab, bool voigt = false) const;
	void Energy (int amrlev, amrex::MultiFab& energy, const amrex::MultiFab& u) const;
	void SetBC(const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_lo,
		   const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_hi)
	{m_bc_lo = a_bc_lo;m_bc_hi = a_bc_hi;};

	void Reflux(int crse_amrlev,
		    MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		    MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs)
	{reflux(crse_amrlev, res, crse_sol, crse_rhs,fine_res, fine_sol, fine_rhs);}
	virtual void FApply (int amrlev, int mglev, MultiFab& out, const MultiFab& in)
	{Fapply(amrlev,mglev,out,in);}

	void Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const;

	inline Set::Vector apply (int amrlev, int mglev,
				  const amrex::FArrayBox &ufab,
				  TArrayBox &C,
				  const amrex::IntVect &m) const;

	void SetTesting(bool a_testing) {m_testing = a_testing;}
protected:

	virtual void Diagonal (int amrlev, int mglev, amrex::MultiFab& diag) override;

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
	//virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const final;
	virtual void FFlux (int amrlev, const MFIter& mfi,
						const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
						const FArrayBox& sol, const int face_only=0) const final;




	virtual int getNComp() const {return AMREX_SPACEDIM;};
	virtual bool isCrossStencil () const { return false; }
	//virtual void Diagonal(bool recompute=false);
	//virtual void normalize (int amrlev, int mglev, MultiFab& mf) const;

private:
	/// \todo: This function needs to be removed. It currently exists only to avoid
	///        breaking other dependencies, but must never be called.
	virtual amrex::Real C(const int , const int , const int , const int , const amrex::IntVect ,
						  const int , const int , const MFIter &) const
	{
		Util::Abort(INFO, "This function is depricated and must not be used.");
		return std::numeric_limits<Set::Scalar>::quiet_NaN();
	};

	/// Simple arrays storing boundary conditions for each component and each face.
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_lo; // m_bc_lo[face][dimension]
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_hi; // m_bc_hi[face][dimension]

	/// This is a multifab-type object containing objects of type
	/// Model::Solid::Elastic::Isotropic::Isotropic
	/// (or some other model type). T is the template argument.
	/// The models contain elastic constants and contain methods for converting strain to stress
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T> > > > > model;



	virtual void averageDownCoeffs ();
	void averageDownCoeffsToCoarseAmrLevel (int flev);
	void averageDownCoeffsSameAmrLevel (int amrlev);

	void FillBoundaryCoeff (amrex::FabArray<amrex::BaseFab<T> >& sigma, const Geometry& geom);

	bool m_testing = false;

};
}

#endif
