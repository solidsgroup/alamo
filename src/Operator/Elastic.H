#ifndef OPERATOR_ELASTIC_H_
#define OPERATOR_ELASTIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>
#include "Set/Set.H"
#include "Operator/Operator.H"
#include "Model/Solid/Solid.H"
#include "Test/Operator/Elastic.H"
#include "BC/Operator/Elastic.H"

using namespace amrex;

namespace Operator
{
/// 
/// \todo Remove Elastic derived classes. They have been replaced with the Model construction.
/// 
template<class T>
class Elastic : public Operator<Grid::Node>
{
	friend Test::Operator::Elastic;
	
	using TArrayBox = amrex::BaseFab<T>;
	using MultiTab  = amrex::FabArray<TArrayBox>;
public:
	enum class BC {Displacement, Traction, Periodic, Neumann}; 
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

	virtual void SetHomogeneous (bool a_homogeneous) override {m_homogeneous = a_homogeneous;}
	void SetModel (T &a_model);
	void SetModel (int amrlev, const amrex::FabArray<amrex::BaseFab<T> >& a_model);
	void SetModel (const amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > & a_model)
	{ for (int ilev = 0; ilev < a_model.size(); ilev++) SetModel(ilev,a_model[ilev]);}
	void SetModel (const Set::Field<T> & a_model)
	{ for (int ilev = 0; ilev < a_model.size(); ilev++) SetModel(ilev,*a_model[ilev]);}

	/// The different types of Boundary Condtiions are listed in the `BC::Operator::Elastic` documentation
	///
	void SetBC (::BC::Operator::Elastic<T> *a_bc) 
	{
		if (this->Geom(0).isAnyPeriodic()) 
		{
			Util::Warning(INFO,"Looks like you're using a periodic domain. \nThat is currently VERY DANGEROUS when using linear elastic solver!!!");
		}
		m_bc = a_bc;
		m_bc_set = true;
	};
	::BC::Operator::Elastic<T> & GetBC()
	{
		return *m_bc;
	}

	/// Compute strain \f$\mathbf{\epsilon}\f$ given the displacement field \f$\mathbf{u}\f$
	/// by
	///   \f[\mathbf{\epsilon}_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})\f]
	/// 
	void Strain (int amrlev, amrex::MultiFab& epsfab, const amrex::MultiFab& ufab, bool voigt = false) const;

	/// Compute stress \f$\mathbf{\sigma}\f$ given the displacement field \f$\mathbf{u}\f$
	/// by
	///   \f[\mathbf{\sigma}_{ij} = \mathbb{C}_{ijkl}\,u_{k,l}\f]
	/// where, \f$\mathbb{C}\f$ is furnished by the templated `model`
	/// 
	void Stress (int amrlev, amrex::MultiFab& sigmafab, const amrex::MultiFab& ufab, bool voigt = false, bool a_homogeneous=false);

	/// Compute energy density \f$\mathbb{W}\f$ given the displacement field \f$\mathbf{u}\f$
	/// by
	///   \f[\mathbb{W} = \mathbb{\Sigma}\mathbb{\Sigma}\mathbf{\sigma}_{ij}\mathbf{\epsilon}_{ij}\f]
	/// where \f$\mathbf{\sigma}\f$ and \f$\mathbf{\epsilon}\f$ are as defined above
	///
	void Energy (int amrlev, amrex::MultiFab& energy, const amrex::MultiFab& u, bool a_homogeneous=false);

	void Energy (int amrlev, amrex::MultiFab& energies, const amrex::MultiFab& u, std::vector<T> models, bool a_homogeneous=false);

	/// This function is depricated and should not be used. Use the other `SetBC` function.
	///
	void SetBC(const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_lo,
		   const std::array<std::array<BC,AMREX_SPACEDIM>,AMREX_SPACEDIM> &a_bc_hi)
	{
		Util::Abort(INFO,"This is now depricated. Use the other SetBC function.");
		m_bc_lo = a_bc_lo;m_bc_hi = a_bc_hi;
	};

	void Reflux(int crse_amrlev,
		    MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		    MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs)
	{reflux(crse_amrlev, res, crse_sol, crse_rhs,fine_res, fine_sol, fine_rhs);}

	void Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const;

	void SetTesting(bool a_testing) {m_testing = a_testing;}
	void SetUniform(bool a_uniform) {m_uniform = a_uniform;}
	
	
protected:

	virtual void Diagonal (int amrlev, int mglev, amrex::MultiFab& diag) override;

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const override final;
	//virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const final;
	virtual void FFlux (int amrlev, const MFIter& mfi,
						const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
						const FArrayBox& sol, const int face_only=0) const final;




	virtual int getNComp() const override {return AMREX_SPACEDIM;};
	virtual bool isCrossStencil () const { return false; }
	//virtual void Diagonal(bool recompute=false);
	//virtual void normalize (int amrlev, int mglev, MultiFab& mf) const;
	virtual void prepareForSolve ()
	{
		if (!m_model_set) Util::Abort(INFO,"Attempting to use operator before calling SetModel!");
		if (!m_bc_set) Util::Warning(INFO,"Attempting to use operator before calling SetBC!");
		Operator<Grid::Node>::prepareForSolve();
	};

private:
	/// Simple arrays storing boundary conditions for each component and each face.
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_lo; // m_bc_lo[face][dimension]
	std::array<std::array<BC,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_hi; // m_bc_hi[face][dimension]

	/// This is a multifab-type object containing objects of type
	/// Model::Solid::Elastic::Isotropic::Isotropic
	/// (or some other model type). T is the template argument.
	/// The models contain elastic constants and contain methods for converting strain to stress
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T> > > > > model;



	virtual void averageDownCoeffs () override;
	void averageDownCoeffsToCoarseAmrLevel (int flev);
	void averageDownCoeffsSameAmrLevel (int amrlev);

	void FillBoundaryCoeff (MultiTab& sigma, const Geometry& geom);

	bool m_testing = false;
	bool m_uniform = false;
	bool m_homogeneous = false;

	::BC::Operator::Elastic<T> *m_bc;

	bool m_model_set = false;
	bool m_bc_set = false;
	
public:
    static void Parse(Elastic<T> & /*value*/, IO::ParmParse & /*pp*/)
    {
    }

};
}

#endif
