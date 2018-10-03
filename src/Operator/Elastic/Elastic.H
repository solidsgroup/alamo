#ifndef OPERATOR_ELASTIC_H_
#define OPERATOR_ELASTIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Set/Set.H"
#include "Operator/Operator.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Elastic/Isotropic/Isotropic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
enum BC	{Displacement, Traction, Periodic}; 

template<class T>
class Elastic : public Operator
{
public:
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
	void Stress (FArrayBox& sigmafab, const FArrayBox& ufab, int amrlev, const MFIter& mfi, bool voigt = false) const;
	void Energy (FArrayBox& energyfab, const FArrayBox& ufab, int amrlev, const MFIter& mfi) const;
	void SetBC(const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_lo,
		   const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_hi)
	{m_bc_lo = a_bc_lo;m_bc_hi = a_bc_hi;};

protected:

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const final;
	virtual void FFlux (int amrlev, const MFIter& mfi,
						const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
						const FArrayBox& sol, const int face_only=0) const final;
	virtual int getNComp() const {return AMREX_SPACEDIM;};
	virtual bool isCrossStencil () const { return false; }
	virtual void normalize (int amrlev, int mglev, MultiFab& mf) const;

private:
	// This is an old depricated method that needs to be eliminated.
	// We keep it around to maintain backward compatibility.
	virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
						  const int amrlev, const int mglev, const MFIter &mfi) const
	{
		Util::Abort("This function is depricated and must not be used.");
	};

	// Simple arrays storing boundary conditions for each component and each face.
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_lo; // m_bc_lo[face][dimension]
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_hi; // m_bc_hi[face][dimension]

	// This is a multifab containing objects of type Model::Solid::Elastic::Isotropic
	// (or some other model type). T is the template argument.
	// The models contain elastic constants and contain methods for converting strain to stress
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T> > > > > model;
};
}
}

#endif
