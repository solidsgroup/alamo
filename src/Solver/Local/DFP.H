#include <eigen3/Eigen/Core>
#include "Set/Set.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
    class DFP
    {
    private:
    
	static constexpr int dim = 8;
    using vector2d = Eigen::Matrix<amrex::Real,dim,1>;
	using matrix22 = Eigen::Matrix<amrex::Real,dim,dim>;

    vector2d getGrad(vector2d x, double dx, const Set::Matrix& sig);
    double secantMethod(double dx, double a1, double a2, double tol, vector2d x_new, vector2d r, const Set::Matrix& sig);
    Set::Scalar f(vector2d x, const Set::Matrix& es);
    vector2d DFPfunc(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig);
    public:
    Set::Matrix relax(const Set::Matrix& sig, double e);
    };
}
}