#ifndef TEST_OPERATOR_ELASTIC_H
#define TEST_OPERATOR_ELASTIC_H

#include <AMReX_Vector.H>

namespace Test
{
namespace Operator
{
class Elastic
{
public:
	Elastic() {};

	/// Create a simulation domain where the coarse level has ncells x ncells ( x ncells)
	/// and nlevels.
	/// Each level is identical to the one below but with a quarter of the area and twice
	/// the resolution.
	/// \todo Add option to specify one of multiple refinement geometries
	void Define(int _ncells,
		    int _nlevels);

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
	int TrigTest(bool verbose,
		     int component = 0,
		     int n = 1,
		     std::string plotfile = "");
  
        /// \todo: write documentation for this test
	int UniaxialTest(bool verbose,
		     int component = 0,
		     int n = 1,
		     std::string plotfile = "");

	int SpatiallyVaryingCTest(int);
	
private:


	int nlevels = 1;
	int ncells = 16;
	const int ref_ratio = 2;
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

	amrex::Vector<std::string> varname = {"solution_exact1", "solution_exact2",
				       "solution_numeric1", "solution_numeric2",
				       "solution_error1", "solution_error2",
				       "rhs_prescribed1", "rhs_prescribed2",
				       "rhs_exact1","rhs_exact2",
				       "rhs_numeric1","rhs_numeric2",
				       "res_exact1","res_exact2",
				       "res_numeric1","res_numeric2",
				       "ghost_force1", "ghost_force2"};

};
}
}





#endif


