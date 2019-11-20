#include "DFP.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
Eigen::Matrix<amrex::Real,8,1> DFP::getGrad(vector2d x, double dx, const Set::Matrix& sig)
{
	vector2d grad = vector2d::Zero();
	for (int i = 0; i < dim; i++)
	{
		vector2d xp = x; vector2d xm = x;
		xp(i) = x(i) + dx;
		xm(i) = x(i) - dx;

		grad(i) = ( f(xp, sig) - f(xm, sig) ) / (2 * dx);
	}
	//Util::Message(INFO,"grad = ", grad);
	return grad;
}

double DFP::secantMethod(double dx, double a1, double a2, double tol, vector2d x_new, vector2d r, const Set::Matrix& sig)
{
	double a = a1; double b = a2;
	while (abs(a - b) >= tol)
	{
		double df1 = (f(x_new + (b + dx) * r, sig) - f(x_new + (b - dx) * r, sig)) / (2 * dx);
		double df2 = (f(x_new + (a + dx) * r, sig) - f(x_new + (a - dx) * r, sig)) / (2 * dx);
		double temp = b - 0.8*((df1*(b - a)) / (df1 - df2));

		a = b;
		b = temp;
	}
	return b;
}

Eigen::Matrix<amrex::Real,8,1> DFP::DFPfunc(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig)
{
	vector2d xnew = x0;
	vector2d xprev = -xnew; xprev(1) -= 100;
	matrix22 Hnew, Hprev;
	Hnew = matrix22::Identity(); Hprev = matrix22::Identity();

	while (abs(xnew.norm() - xprev.norm()) >= tol)
	{
		
		vector2d R = -Hnew * getGrad(xnew, dx, sig);
		R.normalize();
		
		double temp = secantMethod(dx, alpha1, alpha2, tol, xnew, R, sig);
		//Util::Message(INFO,"temp = ", temp);
		xprev = xnew;
		xnew = xnew + temp * R;

		vector2d gamma = getGrad(xnew, dx, sig) - getGrad(xprev,dx, sig);
		vector2d del = xnew - xprev;

		Hprev = Hnew;
		Hnew = Hprev + (del * del.transpose()) / (del.transpose() * del) 
		- (Hprev * (gamma*gamma.transpose()) * Hprev) / (gamma.transpose() * Hprev * gamma);
		//Util::Message(INFO,"x = ", xnew.transpose());
		
	/*
		vector2d p = -Hnew * getGrad(xnew, dx, sig);

		double temp = secantMethod(dx, alpha1, alpha2, tol, xnew, p, sig);
		vector2d s = temp*p;
		xprev = xnew;
		xnew = xnew + s;
		vector2d y = getGrad(xnew, dx, sig) - getGrad(xprev,dx, sig);

		Hprev = Hnew;
		//Hnew = Hprev + ( (s.transpose() * y + y.transpose() * Hprev * y) * (s*s.transpose()) ) /  (s.transpose()*y*s.transpose()*y);// - (Hprev * y * s.transpose() + s * y.transpose() * Hprev) / (s.transpose() * y); 
    */
	}
	return xnew;
}

Set::Scalar DFP::f(vector2d x, const Set::Matrix& es)
{
	//Set::Matrix epsilon = es; 
	//					 epsilon(0,1) = x(0); epsilon(0,2) = x(1);
	//epsilon(1,0) = x(2); epsilon(1,1) = x(3); epsilon(1,2) = x(4);
	//epsilon(2,0) = x(5); epsilon(2,1) = x(6); epsilon(2,2) = x(7);
	//Set::Matrix temp = epsilon - esp;
	//Set::Scalar s = W(temp);
	//Util::Message(INFO,"W = ",s);
	//return s;
}
Set::Matrix DFP::relax(const Set::Matrix& sig, const double e) 
{
	vector2d x = vector2d::Zero(); x(3) = -e-1e-2; x(7) = -e-1e-2;
	vector2d xprev = vector2d::Zero(); xprev(1) = 100;
	Set::Matrix ep = Set::Matrix::Zero();

	while(abs(x.norm() - xprev.norm()) >= 1e-5)
	{
		xprev = x;
		//x = DFP(x,1e-5,0.1,0.8,1e-5,sig);
		//Util::Message(INFO,"ffff= ", counter);
	}

	ep(0,0) = e;	ep(0,1) = x(0); ep(0,2) = x(1);
	ep(1,0) = x(2); ep(1,1) = x(3); ep(1,2) = x(4);
	ep(2,0) = x(5); ep(2,1) = x(6); ep(2,2) = x(7);
	return ep;
}

}
}