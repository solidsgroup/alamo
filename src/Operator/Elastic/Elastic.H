#ifndef AMREX_OPERATOR_ELASTIC_H_
#define AMREX_OPERATOR_ELASTIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Set/Set.H"
#include "Operator/Operator.H"
#include "Model/Solid/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
enum BC	{Displacement, Traction, Periodic}; 


class Solid
{
public:
	Solid() {} ;
	Set::Matrix operator () (Set::Matrix &eps);
	{
	};
	Set::Scalar operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps);

	friend Solid operator + (Solid lhs, const Solid rhs)
	{lhs.mu += rhs.mu; lhs.lambda += rhs.lambda;}

	friend Solid operator - (Solid lhs, const Solid rhs)
	{lhs.mu -= rhs.mu; lhs.lambda -= rhs.lambda;}
	
	friend Solid operator / (Solid lhs, const Set::Scalar alpha)
	{lhs.mu /= alpha; lhs.lambda /= alpha;}

	
private:
	Set::Scalar mu, lambda;
};

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


	// void SetEigenstrain(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &);
	// void SetEigenstrain(amrex::Vector<amrex::MultiFab> &);

	// void AddEigenstrainToRHS (FArrayBox& esrhs,
	// 								int amrlev, const MFIter& mfi) const;
	// void AddEigenstrainToRHS (amrex::Vector<amrex::MultiFab>& rhsfab) const;

	void Stress (FArrayBox& sigmafab,
		     const FArrayBox& ufab,
		     int amrlev, const MFIter& mfi,
		     bool voigt = false) const;
	void Energy (FArrayBox& energyfab,
				 const FArrayBox& ufab,
				 int amrlev, const MFIter& mfi) const;
	void SetBC(const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_lo,
		   const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_hi)
	{
		AMREX_D_TERM(m_bc_xlo = a_bc_lo[0];, m_bc_ylo = a_bc_lo[1];, m_bc_zlo = a_bc_lo[2];);
		AMREX_D_TERM(m_bc_xhi = a_bc_hi[0];, m_bc_yhi = a_bc_hi[1];, m_bc_zhi = a_bc_hi[2];);
	};

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
	virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
						  const int amrlev, const int mglev, const MFIter &mfi) const = 0;
	bool usingEigenstrain = false;

	std::array<BC,AMREX_SPACEDIM> AMREX_D_DECL(m_bc_xlo, m_bc_ylo, m_bc_zlo);
	std::array<BC,AMREX_SPACEDIM> AMREX_D_DECL(m_bc_xhi, m_bc_yhi, m_bc_zhi);

	Model::Solid::Elastic::Elastic *model;

	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<Solid> > > > > coeff;

	//std::array<BC,AMREX_SPACEDIM> AMREX_D_DECL{bc_x, bc_y, bc_z};
};
}
}
#endif
