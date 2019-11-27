#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
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
        static constexpr double dx = 1e-5;

        int rows,cols;
        Set::Matrix sig, esp;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H, Hpre;
        Eigen::Matrix<double, Eigen::Dynamic, 1> a, b, c, d;
        
        //using vector2d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
        //using matrix22 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

        Eigen::Matrix<double, Eigen::Dynamic, 1> getGrad();
        double f(Eigen::Matrix<double, Eigen::Dynamic, 1> x, const Set::Matrix& es);
    public:
        DFP(int _rows, int _cols,double c11, double c12, double c44 )
    :
    cp(c11,c12,c44)
        {
            rows = _rows; cols = _cols;
            H.resize(rows,cols); Hpre.resize(rows, cols);
            a.resize(rows); b.resize(rows); c.resize(rows); d.resize(rows);
        }
        Model::Solid::CrystalPlastic::CrystalPlastic cp;
        };
    }
}