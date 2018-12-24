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
/// **Stress divergence operator definition**
///
/// The stress divergence equation is
/// \f[ f_i = D_{\Omega}(\mathbf{u})_i
///         = \frac{\partial}{\partial x_j}\mathbb{C}_{ijkl}\frac{\partial u_k}{\partial x_l}
///         = \mathbb{C}_{ijkl}u_{k,lj} + \mathbb{C}_{ijkl,j}u_{k,l}\f]
/// Where \f$\mathbf{f}\f$ is the vector of nodal forces, \f$\mathbf{u}\f$ the vector of
/// nodal displacements, and \f$\mathbb{C}_{ijkl}\f$ the fourth-order elasticity tensor.
/// 
/// **Boundary operator**
///
/// If a point is located on the boundary then a boundary operator is used rather than
/// the stress divergence operator.
///
/// *Dirichlet/Essential/Displacement BC Operator*
/// For a Dirichlet point the Dirichlet boundary operator is used, \f$D_{\Omega_1}\f$.
/// This is nothing other than the identity operator, i.e.
/// \f[ D_{\Omega_1}(\mathbf{u}) = \mathbf{u}, \f]
/// that is, it is the identity.
/// The _value_ of the displacement is thus determined by the value of the corresponding point
/// in the right hand side during the solve.
/// 
/// *Natural/Traction BC Operator*
/// This operator is _similar_ (but not the same as) a typical Neumann boundary condition.
/// The operator, denoted \f$D_{\Omega_2}\f$ is defined as
/// \f[ D_{\Omega_2}(\mathbf{u}) = \mathbb{C}\nabla\mathbf{u} \cdot \mathbf{n}
///                              = \sigma\cdot\mathbf{n}\f]
/// The output of this operator is a surface traction. The value of the surface traction is
/// specified in the right hand side during the solve.
///
/// **Solid Model Template Class T**
///
/// This operator is templated with class T.
/// T should be of type Model::Solid::Solid or derived.
/// These models are data structures containing all necessary elastic constants to compute the
/// stress tensor \f$\sigma\f$ given \f$\nabla\mathbf{u}\f$.
/// The models are stored in an amrex::FabArray of amrex::BaseFab<T>.
/// Models have basic arithmetic operators defined, and use the parentheses operator to
/// do the calculation.
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
	
	inline Set::Matrix flux (int amrlev, int mglev,
				 const amrex::FArrayBox &ufab,
				 TArrayBox &C,
				 const amrex::IntVect &m,
				 std::array<Boundary,AMREX_SPACEDIM>) const;

	void SetTesting(bool a_testing) {m_testing = a_testing;}
protected:

	virtual void Diagonal (int amrlev, int mglev, amrex::MultiFab& diag) override;

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
	//virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const final;
	virtual void FFlux (int amrlev, const MFIter& mfi,
						const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
						const FArrayBox& sol, const int face_only=0) const final;



	///
	/// This is a three step process.
	/// **Note** it is implemented in 2D only and all documentation references the 2D implementation
	/// where dimensionality is not general.
	/// 
	/// 1. Restriction.
	///    Updates the residual on the coarse fab with the residual as calculated on the fine fab.
	///    Restriction takes place only on the interior nodes.
	///    In 2D the restriction scheme is:
	///    \f{eqnarray*}{	
	///       res_{c}(I,J)    &=& 4\,res_f(i,j) 
	///                    \\ &+& 2\,\Big[res_f(i+1,j) + res_f(i-1,j) + res_f(i,j+1) + res_f(i,j-1)\Big] 
	///                    \\ &+& \Big[res_f(i+1,j+1) + res_f(i+1,j-1) + res_f(i-1,j+1) + res_f(i-1,j-1)\Big]
	///    \f}
	///

	virtual void reflux (int crse_amrlev,
			     MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const;

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
