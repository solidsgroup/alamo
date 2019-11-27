#include "DFP.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
	Eigen::Matrix<double, Eigen::Dynamic, 1> DFP::getGrad()
	{
		/*
		a.Zero();
		for (int i = 0; i < rows; i++)
		{
			b = x; 
			c = x;
			b(i) = x(i) + dx;
			c(i) = x(i) - dx;

			a(i) = ( f(b, sig) - f(c, sig) ) / (2 * dx);
		}
		//Util::Message(INFO,"grad = ", grad);
		return a;*/
	}

	double DFP::f(Eigen::Matrix<double, Eigen::Dynamic, 1> x, const Set::Matrix& es)
{/*
	Set::Matrix epsilon = es; 
						 epsilon(0,1) = x(0); epsilon(0,2) = x(1);
	epsilon(1,0) = x(2); epsilon(1,1) = x(3); epsilon(1,2) = x(4);
	epsilon(2,0) = x(5); epsilon(2,1) = x(6); epsilon(2,2) = x(7);
	Set::Matrix temp = epsilon - esp;
	Set::Scalar s = cp.W(temp);
	//Util::Message(INFO,"W = ",s);
	return s;*/
}

}
}