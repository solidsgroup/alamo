#ifndef TEST_OPERATOR_ELASTIC_H
#define TEST_OPERATOR_ELASTIC_H

#include <streambuf>
#include <map>
#include <complex>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>

#include"Util/Util.H"

namespace Test
{
namespace Operator
{
class Elastic
{
public:
	enum Grid {None, X, Y, Z, YZ, ZX, XY, XYZ};

public:
	Elastic() {};
	~Elastic() {};

	/// Create a simulation domain where the coarse level has ncells x ncells ( x ncells)
	/// and nlevels.
	/// Each level is identical to the one below but with a quarter of the area and twice
	/// the resolution.
	/// \todo Add option to specify one of multiple refinement geometries
	void Define(const amrex::IntVect _ncells,
		    const int _nlevels,
		    const int _dim = AMREX_SPACEDIM,
		    const Grid config = Grid::XYZ);

	void Define(const int _ncells,
		    const int _nlevels,
		    const int _dim = AMREX_SPACEDIM,
		    const Grid config = Grid::XYZ);

	void Define(const int _ncells,
		    const int _nlevels)
	{
		Define(_ncells,_nlevels,AMREX_SPACEDIM,Grid::XYZ);
	};

	/// Impose a set of prescribed values for \f$\mathbf{u}\f$.
	/// Compute right hand side, then use Operator::Reflux to modify the residual \f$\mathbf{r}=\mathbf{b}-D\mathbf{u}\f$.
	/// Residual should always be zero unless something goes wrong with Reflux.
	/// If this test fails, check Reflux.
	int RefluxTest(int);

	/// Compute the exact solution of the governing equation
	///   \f[C_{ijkl} u_{k,jl} + b_i = 0\f]
	/// Where
	///   \f[C_{ijkl} = \alpha \delta_{ik} \delta_{jl}\f]
	/// (see Model::Solid::LinearElastic::Laplacian for implementation of this model)
	///
	/// Consider a simulation domain \f$\Omega = [0, 1]^2\f$
	///
	/// Expand \f$\mathbf{u}\f$ and \f$\mathbf{b}\f$ out as trigonometric series
	///       \f[ u_i = u_i^{mn} \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f] 
	///       \f[ b_i = b_i^{mn} \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f]
	///
	/// Then the second derivative of displacement is
	///       \f[ u_{i,jj} = - u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f]
	/// and the pointwise force is 
	///       \f[ C_{ijkl} u_{k,jl} = \alpha \delta_{ik} \delta_{jl} u_{k,jl} = \alpha u_{i,jj} \f]
	///                     
	/// Substituting into the governing equation gives
	///       \f[-\alpha u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m) \sin(\pi\ y\ n) 
	///          + b_i^{mn} \sin(\pi\ x\ m)\sin(\pi\ y\ n) 
	///                     = 0\f]
	/// or
	///     \f[\alpha u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m)\sin(\pi y n) = b_i^{mn} \sin(\pi\ x\ m)\sin(\pi\ y\ n)\f]
	///
	/// Using orthognonality of the eigenfunctions gives
	///     \f[\alpha u_i^{mn} \pi^2 (m^2 + n^2) = b_i^{mn}\f]
	/// or, simplifying
	///     \f[u_i^{mn}  = \frac{b_i^{mn}}{(\alpha\ \pi^2\ (m^2 + n^2))}\f]
	/// 
	///
	int TrigTest(int verbose,
		     int component = 0,
		     int n = 1,
		     std::string plotfile = "");
	int mytest() {return 8;};

	// Setter functions
	void setMaxCoarseningLevel(int a_maxCoarseningLevel) {m_maxCoarseningLevel = a_maxCoarseningLevel;}
	void setFixedIter(int a_fixedIter) {m_fixedIter = a_fixedIter;}
	void setMaxIter(int a_maxIter) {m_maxIter = a_maxIter;}
	void setMaxFmgIter(int a_maxFmgIter) {m_maxFmgIter = a_maxFmgIter;}
	void setBottomMaxIter(int a_bottomMaxIter) {m_bottomMaxIter = a_bottomMaxIter;}
	void setBounds(AMREX_D_DECL(Set::Scalar a_x,Set::Scalar a_y, Set::Scalar a_z))
		       {m_bounds = {AMREX_D_DECL(a_x,a_y,a_z)};}
	void setBounds(std::array<Set::Scalar,AMREX_SPACEDIM> a_bounds) {m_bounds = a_bounds;}

private: // Private member functions

	void Generate();

	void WritePlotFile(std::string plotfile, std::vector<int> nghost)
	{
		const int output_comp = varname.size();
		amrex::Vector<amrex::MultiFab> plotmf(nlevels);
		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			amrex::BoxArray plotgrids = ngrids[ilev];
			plotgrids.grow(nghost[ilev]);
			amrex::Box ndomain = geom[ilev].Domain(); ndomain.convert(amrex::IntVect::TheNodeVector());
			plotmf[ilev].define(plotgrids, dmap[ilev], output_comp, nghost[ilev]);
			plotmf[ilev].setVal(0.0);

			for (amrex::MFIter mfi(plotmf[ilev]); mfi.isValid(); ++mfi)
			{
				amrex::Box bx = mfi.validbox();
				bx &= ndomain; // Intersect the box with the domain (that's what &= does)
				
				plotmf[ilev][mfi].copy(solution_exact[ilev][mfi]  ,bx, 0, bx, 0*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(solution_numeric[ilev][mfi],bx, 0, bx, 1*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(solution_error[ilev][mfi]  ,bx, 0, bx, 2*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(rhs_prescribed[ilev][mfi]  ,bx, 0, bx, 3*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(rhs_exact[ilev][mfi]       ,bx, 0, bx, 4*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(rhs_numeric[ilev][mfi]     ,bx, 0, bx, 5*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(res_exact[ilev][mfi]       ,bx, 0, bx, 6*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(res_numeric[ilev][mfi]     ,bx, 0, bx, 7*AMREX_SPACEDIM, AMREX_SPACEDIM); 
				plotmf[ilev][mfi].copy(ghost_force[ilev][mfi]     ,bx, 0, bx, 8*AMREX_SPACEDIM, AMREX_SPACEDIM); 
			}
		}
		amrex::WriteMultiLevelPlotfile(plotfile, nlevels, amrex::GetVecOfConstPtrs(plotmf),
					       varname, geom, 0.0, amrex::Vector<int>(nlevels, 0),
					       amrex::Vector<amrex::IntVect>(nlevels, amrex::IntVect{ref_ratio}));
	}		
	void WritePlotFile(std::string plotfile, int nghost = 0)
	{
		std::vector<int> ngs(nlevels,nghost);
		WritePlotFile(plotfile,ngs);
	}

private: // Private member variables
	int dim = 1;

	int nlevels = 1;
	amrex::IntVect ncells;
	const int ref_ratio = 2;

	Grid m_config = Grid::XYZ;

 	amrex::Vector<amrex::Geometry> 		   geom;
 	amrex::Vector<amrex::BoxArray> 		   cgrids, ngrids;
 	amrex::Vector<amrex::DistributionMapping>  dmap;
	
 	amrex::Vector<amrex::MultiFab>             solution_exact;
 	amrex::Vector<amrex::MultiFab>             solution_numeric;
 	amrex::Vector<amrex::MultiFab>             solution_error;

 	amrex::Vector<amrex::MultiFab>             rhs_prescribed;
 	amrex::Vector<amrex::MultiFab>             rhs_exact;
 	amrex::Vector<amrex::MultiFab>             rhs_numeric;

 	amrex::Vector<amrex::MultiFab>             res_exact;
 	amrex::Vector<amrex::MultiFab>             res_numeric;
 	amrex::Vector<amrex::MultiFab>             ghost_force;

	amrex::Vector<std::string> varname = {AMREX_D_DECL("solution_exact1", "solution_exact2","solution_exact3"),
					      AMREX_D_DECL("solution_numeric1", "solution_numeric2","solution_numeric3"),
					      AMREX_D_DECL("solution_error1", "solution_error2","solution_error3"),
					      AMREX_D_DECL("rhs_prescribed1", "rhs_prescribed2","rhs_prescribed3"),
					      AMREX_D_DECL("rhs_exact1","rhs_exact2","rhs_exact2"),
					      AMREX_D_DECL("rhs_numeric1","rhs_numeric2","rhs_numeric3"),
					      AMREX_D_DECL("res_exact1","res_exact2","res_exact3"),
					      AMREX_D_DECL("res_numeric1","res_numeric2","res_numeric3"),
					      AMREX_D_DECL("ghost_force1", "ghost_force2","ghost_force2")};

	std::array<Set::Scalar,AMREX_SPACEDIM> m_bounds = {AMREX_D_DECL(1.0,1.0,1.0)};
	int m_maxCoarseningLevel = 1000000;
	int m_fixedIter       = 1000000;
	int m_maxIter         = 1000000;
	int m_maxFmgIter      = 1000000;
	int m_bottomMaxIter   = 200;


};
}
}





#endif


