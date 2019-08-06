#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"

namespace Model
{
namespace Solid
{
namespace CrystalPlastic
{

CrystalPlastic::CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
{
	initializeSlip();
	define(C11, C12, C44, R);

	
}
CrystalPlastic::CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	initializeSlip();
	define(C11, C12, C44, phi1, Phi, phi2);
}

void
CrystalPlastic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	Eigen::Matrix3d m;
	m =     Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
	define(C11,C12,C44,m);
}
void CrystalPlastic::initializeSlip()
{
	Set::Vector n1 = {1,1,1};
	Set::Vector n2 = {-1,-1,1};
	Set::Vector n3 = {-1,1,1};
	Set::Vector n4 = {1,-1,1};
	
	Set::Vector s11 = {0,-1,1};
	Set::Vector s12 = {1,0,-1};
	Set::Vector s13 = {-1,1,0};
	Set::Vector s21 = {0,1,1};
	Set::Vector s22 = {-1,0,-1};
	Set::Vector s23 = {1,-1,0};
	Set::Vector s31 = {0,-1,1};
	Set::Vector s32 = {-1,0,-1};
	Set::Vector s33 = {1,1,0};
	Set::Vector s41 = {0,1,1};
	Set::Vector s42 = {1,0,-1};
	Set::Vector s43 = {-1,-1,0};

	this->slp1.n = n1; this->slp1.s = s11;
	this->slp2.n = n1; this->slp2.s = s12;
	this->slp3.n = n1; this->slp3.s = s13;

	this->slp4.n = n2; this->slp4.s = s21;
	this->slp5.n = n2; this->slp5.s = s22;
	this->slp6.n = n2; this->slp6.s = s23;

	this->slp7.n = n3; this->slp7.s = s31;
	this->slp8.n = n3; this->slp8.s = s32;
	this->slp9.n = n3; this->slp9.s = s33;

	this->slp10.n = n4; this->slp10.s = s41;
	this->slp11.n = n4; this->slp11.s = s42;
	this->slp12.n = n4; this->slp12.s = s43;

	slipSystem[0] = slp1;
	slipSystem[1] = slp2;
	slipSystem[2] = slp3; 
	slipSystem[3] = slp4;
	slipSystem[4] = slp5;
	slipSystem[5] = slp6; 
	slipSystem[6] = slp7;
	slipSystem[7] = slp8;
	slipSystem[8] = slp9; 
	slipSystem[9] = slp10;
	slipSystem[10] = slp11; 
	slipSystem[11] = slp12;
	for(int i = 0; i < 12; i++)
	{
		slipSystem[i].on = false;
	}
	
}
double CrystalPlastic::CalcSSigN (const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig) 
{
	double a;
	a = ss.transpose()*sig*nn;
	return a;
}
void CrystalPlastic::GetActivePlains(const Set::Matrix sig)
{
	double a = 0;
	for(int i = 0; i <= 12; i++)
	{
		//calculate a = abs(s*sigma*n)
		a = abs(CalcSSigN(slipSystem[i].s,slipSystem[i].n, sig));
		
		if(a > Tcrss) slipSystem[i].on = true;		
		else slipSystem[i].on = false;
	}	
}

Set::Scalar CrystalPlastic::GetGammaDot(const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig)
{
	double gamma;
	double power = pow( (abs(CalcSSigN( ss,nn,sig )) /Tcrss ) ,n );
	gamma = gammadot0*power;
	return gamma;
}

void CrystalPlastic::AdvanceEsp( const Set::Matrix sig)
{
	GetActivePlains(sig);
	Set::Matrix temp = Set::Matrix::Zero();
	
	for(int i = 0; i < 12; i++)
	{
		if(slipSystem[i].on)
		{
			double a = CalcSSigN(slipSystem[i].s,slipSystem[i].n,sig);
			double gammadot = GetGammaDot(slipSystem[i].s,slipSystem[i].n,sig);
			int sign = sgn(a);
			temp += gammadot*sign*slipSystem[i].s*slipSystem[i].n.transpose();
		}
	}
	// Euler integration
	esp = esp + temp*dt;
}
void CrystalPlastic::update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt)
{
	for(double t = 0.0; t < _dt; t += dt)
	{
		AdvanceEsp(sigma);
		sigma = UpdateSigma(es);
	}
}
Set::Matrix CrystalPlastic::UpdateSigma(const Set::Matrix es)
{
	Set::Matrix temp = (es - esp);
	Set::Matrix sigma = operator()(temp); 
	return sigma;
}
Set::Matrix CrystalPlastic::GetEsp() const
{
	return esp;
}

void CrystalPlastic::Setdt(double _dt)
{
	dt = _dt;
}
//----------DFP Functions-----------//
Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,AMREX_SPACEDIM-1> CrystalPlastic::Jacobi(vector2d x, const Set::Matrix& es, double dx)
{
	matrix22 val;
	val(0,0) = C[ 1]*es(0,0) + C[ 6]*(x(0) + dx) + C[ 7]*x(1) + 2.0 * ( C[ 8]*es(1,2) + C[ 9]*es(2,0) + C[10]*es(0,1));
	val(0,1) = C[ 1]*es(0,0) + C[ 6]*x(0) + C[ 7]*(x(1) + dx) + 2.0 * ( C[ 8]*es(1,2) + C[ 9]*es(2,0) + C[10]*es(0,1));

	val(1,0) = C[ 2]*es(0,0) + C[ 7]*(x(0) + dx) + C[11]*x(1) + 2.0 * ( C[12]*es(1,2) + C[13]*es(2,0) + C[14]*es(0,1));
	val(1,1) = C[ 2]*es(0,0) + C[ 7]*x(0) + C[11]*(x(1) + dx) + 2.0 * ( C[12]*es(1,2) + C[13]*es(2,0) + C[14]*es(0,1));

	return val;
}
Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::newtonfunc(vector2d x, const Set::Matrix& es)
{
	vector2d a;
	a(1) = C[ 1]*es(0,0) + C[ 6]*x(0) + C[ 7]*x(1) + 2.0 * ( C[ 8]*es(1,2) + C[ 9]*es(2,0) + C[10]*es(0,1));
	a(2) = C[ 2]*es(0,0) + C[ 7]*x(0) + C[11]*x(1) + 2.0 * ( C[12]*es(1,2) + C[13]*es(2,0) + C[14]*es(0,1));
	return a;
}
Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::newtonrap(double dx, vector2d x0, const Set::Matrix& es)
{
	matrix22 a;
	double tol = 1e-5;
	vector2d xnew, xpre, temp;
	xnew = x0; xpre = x0; xpre(1) = x0(1) + 10;
	while (abs(xnew.norm() - xpre.norm()) >= tol)
	{
		a = Jacobi( x0, es, dx );
		
		temp = xnew - 0.5*a.inverse()*newtonfunc(xnew,es);
		Util::Message(INFO,"J = ", temp);
		xpre = xnew;
		xnew = temp;
	}
	return xnew;
}

Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::DFP(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig)
{
	vector2d xnew = x0;
	vector2d xprev = -xnew; xprev(1) -= 100;
	matrix22 Hnew, Hprev;
	Hnew = matrix22::Identity(); Hprev = matrix22::Zero();

	while (abs(xnew.norm() - xprev.norm()) >= tol)
	{
		vector2d R = -Hnew * getGrad(xnew, dx, sig);
		R.normalize();
		Util::Message(INFO,"R = ", R.transpose());
		
		double temp = secantMethod(dx, alpha1, alpha2, tol, xnew, R, sig);
		Util::Message(INFO,"temp = ", temp);
		xprev = xnew;
		xnew = xnew + temp * R;

		vector2d gamma = getGrad(xnew, dx, sig) - getGrad(xprev,dx, sig);
		vector2d del = xnew - xprev;

		Hprev = Hnew;
		Hnew = Hprev + (del * del.transpose()) / (del.transpose() * del) - (Hprev * (gamma*gamma.transpose()) * Hprev) / (gamma.transpose() * Hprev * gamma);
		Util::Message(INFO,"x = ", xnew.transpose());
		abort();
	}
	return xnew;
}

double CrystalPlastic::secantMethod(double dx, double a1, double a2, double tol, vector2d x_new, vector2d r, const Set::Matrix& sig)
{
	double a = a1; double b = a2;
	while (abs(a - b) >= tol)
	{
		double df1 = (f(x_new + (b + dx) * r, sig) - f(x_new + (b - dx) * r, sig)) / (2 * dx);
		double df2 = (f(x_new + (a + dx) * r, sig) - f(x_new + (a - dx) * r, sig)) / (2 * dx);
		double temp = b - 0.2*((df1*(b - a)) / (df1 - df2));

		a = b;
		b = temp;
	}
	return b;
}

Set::Scalar CrystalPlastic::f(vector2d x, const Set::Matrix& es)
{
	//double a, b;
	//Set::Matrix sig = (es - esp);
	//a = C[ 1]*sig(0,0) + C[ 6]*x(0) + C[ 7]*x(1) + 2.0 * ( C[ 8]*sig(1,2) + C[ 9]*sig(2,0) + C[10]*sig(0,1));
	//b = C[ 2]*sig(0,0) + C[ 7]*x(0) + C[11]*x(1) + 2.0 * ( C[12]*sig(1,2) + C[13]*sig(2,0) + C[14]*sig(0,1));
	//Set::Scalar val = w1*pow(a,2) + w2*pow(b,2);
	Set::Matrix epsilon = es; epsilon(1,1) = x(0); epsilon(2,2) = x(1);
	Set::Matrix sig = operator()(epsilon);
	Util::Message(INFO,"f = ",0.5 * (sig*es).trace() );
	return 0.5*(sig*es).trace();
}

Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::getGrad(vector2d x, double dx, const Set::Matrix& sig)
{
	vector2d grad = vector2d::Zero();
	for (int i = 0; i < dim; i++)
	{
		vector2d xp = x; vector2d xm = x;
		xp(i) = x(i) + dx;
		xm(i) = x(i) - dx;

		grad(i) = ( f(xp, sig) - f(xm, sig) ) / (2 * dx);
	}
	Util::Message(INFO,"grad = ", grad);
	return grad;
}
//--------------------------------//
Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::reflux(const Set::Matrix& sig, const double e) 
{
	
	vector2d x; x(0) = -e+1; x(1) = -e + 1; w1 = 100.0; w2 = 100.0;
	vector2d xprev = vector2d::Zero(); xprev(1) = 100;
	int counter = 1;
	while(abs(x.norm() - xprev.norm()) >= 1e-6)
	{
		xprev = x;
		x = DFP(x,2e-6,0.8,1.5,2e-6,sig);
		if(counter % 1 == 0)
		{
			w1 = w1 * 1.1;
			w2 = w2 * 1.1;
		}
		counter++;
		//Util::Message(INFO,"ffff= ", counter);
	}
	return x;
	
//vector2d x; x(0) = -e; x(1) = -e;
//	vector2d a = newtonrap(1e-5,x, sig);
//	return a;
}
void
CrystalPlastic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
{
  
	amrex::Real Ctmp[3][3][3][3];
	amrex::Real Crot[3][3][3][3];


	for(int i = 0; i < 3; i++) 
		for(int j = 0; j < 3; j++) 
			for(int k = 0; k < 3; k++) 
				for(int l = 0; l < 3; l++)
				{
					if(i == j && j == k && k == l)  Ctmp[i][j][k][l] = C11;
					else if (i==k && j==l) Ctmp[i][j][k][l] = C44;
					else if (i==j && k==l) Ctmp[i][j][k][l] = C12;
					else Ctmp[i][j][k][l] = 0.0;
				}

	for(int p = 0; p < 3; p++) 
		for(int q = 0; q < 3; q++) 
			for(int s = 0; s < 3; s++) 
				for(int t = 0; t < 3; t++)
				{
					Crot[p][q][s][t] = 0.0;
					for(int i = 0; i < 3; i++) 
						for(int j = 0; j < 3; j++) 
							for(int k = 0; k < 3; k++) 
								for(int l = 0; l < 3; l++) 
									Crot[p][q][s][t] += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
				}

	C[ 0] = Crot[0][0][0][0]; C[ 1] = Crot[0][0][1][1]; C[ 2] = Crot[0][0][2][2]; C[ 3] = Crot[0][0][1][2]; C[ 4] = Crot[0][0][2][0]; C[ 5] = Crot[0][0][0][1];
	/**/                      C[ 6] = Crot[1][1][1][1]; C[ 7] = Crot[1][1][2][2]; C[ 8] = Crot[1][1][1][2]; C[ 9] = Crot[1][1][2][0]; C[10] = Crot[1][1][0][1];
	/**/                                                C[11] = Crot[2][2][2][2]; C[12] = Crot[2][2][1][2]; C[13] = Crot[2][2][2][0]; C[14] = Crot[2][2][0][1];
	/**/                                                                          C[15] = Crot[1][2][1][2]; C[16] = Crot[1][2][2][0]; C[17] = Crot[1][2][0][1];
	/**/                                                                                                    C[18] = Crot[2][0][2][0]; C[19] = Crot[2][0][0][1];
	/**/                                                                                                                              C[20] = Crot[0][1][0][1];
}

void
CrystalPlastic::Randomize()
{
	Set::Scalar C11 = Util::Random();
	Set::Scalar C12 = Util::Random();
	Set::Scalar C44 = Util::Random();

	Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();

	define(C11,C12,C44,phi1,Phi,phi2);
}


#if AMREX_SPACEDIM==2
Set::Matrix
CrystalPlastic::operator () (Set::Matrix &gradu) const
{
	Set::Matrix eps = 0.5*(gradu + gradu.transpose());
	Set::Matrix sig = Set::Matrix::Zero();
		
	sig(0,0) = C[ 0]*eps(0,0) + C[ 1]*eps(1,1) /*+ C[ 2]*eps(2,2)*/ + 2.0*( /*C[ 3]*eps(1,2) + C[ 4]*eps(2,0) +*/ C[ 5]*eps(0,1));
	sig(1,1) = C[ 1]*eps(0,0) + C[ 6]*eps(1,1) /*+ C[ 7]*eps(2,2)*/ + 2.0*( /*C[ 8]*eps(1,2) + C[ 9]*eps(2,0) +*/ C[10]*eps(0,1));
	sig(0,1) = C[ 5]*eps(0,0) + C[10]*eps(1,1) /*+ C[14]*eps(2,2)*/ + 2.0*( /*C[17]*eps(1,2) + C[19]*eps(2,0) +*/ C[20]*eps(0,1));
	sig(1,0) = sig(0,1);

	return sig;
}
Set::Vector
CrystalPlastic::operator () (std::array<Set::Matrix,2> &gradgradu)
{
	std::array<Set::Matrix,2> gradeps;
	gradeps[0](0,0) = gradgradu[0](0,0);
	gradeps[0](0,1) = gradgradu[0](0,1);
	gradeps[0](1,0) = 0.5*(gradgradu[0](1,0) + gradgradu[1](0,0));
	gradeps[0](1,1) = 0.5*(gradgradu[0](1,1) + gradgradu[1](0,1));
	gradeps[1](0,0) = 0.5*(gradgradu[1](0,0) + gradgradu[0](1,0));
	gradeps[1](0,1) = 0.5*(gradgradu[1](0,1) + gradgradu[0](1,1));
	gradeps[1](1,0) = gradgradu[1](1,0);
	gradeps[1](1,1) = gradgradu[1](1,1);
	Set::Vector f = Set::Vector::Zero();

	f(0) =  C[ 0]*gradeps[0](0,0) + C[ 1]*gradeps[1](1,0) + 2.0 * ( C[ 5]*gradeps[0](1,0)) + 
		C[ 5]*gradeps[0](0,1) + C[10]*gradeps[1](1,1) + 2.0 * ( C[20]*gradeps[0](1,1));

	f(1) =  C[ 5]*gradeps[0](0,0) + C[10]*gradeps[1](1,0) + 2.0 * ( C[20]*gradeps[0](1,0)) +
		C[ 1]*gradeps[0](0,1) + C[ 6]*gradeps[1](1,1) + 2.0 * ( C[10]*gradeps[0](1,1));

	return f;
}
#elif AMREX_SPACEDIM==3
Set::Matrix
CrystalPlastic::operator () (Set::Matrix &gradu) const
{
	Set::Matrix eps = 0.5*(gradu + gradu.transpose());
	Set::Matrix sig = Set::Matrix::Zero();
	
	sig(0,0) = C[ 0]*eps(0,0) + C[ 1]*eps(1,1) + C[ 2]*eps(2,2) + 2.0 * ( C[ 3]*eps(1,2) + C[ 4]*eps(2,0) + C[ 5]*eps(0,1));
	sig(1,1) = C[ 1]*eps(0,0) + C[ 6]*eps(1,1) + C[ 7]*eps(2,2) + 2.0 * ( C[ 8]*eps(1,2) + C[ 9]*eps(2,0) + C[10]*eps(0,1));
	sig(2,2) = C[ 2]*eps(0,0) + C[ 7]*eps(1,1) + C[11]*eps(2,2) + 2.0 * ( C[12]*eps(1,2) + C[13]*eps(2,0) + C[14]*eps(0,1));
	sig(1,2) = C[ 3]*eps(0,0) + C[ 8]*eps(1,1) + C[12]*eps(2,2) + 2.0 * ( C[15]*eps(1,2) + C[16]*eps(2,0) + C[17]*eps(0,1));
	sig(2,0) = C[ 4]*eps(0,0) + C[ 9]*eps(1,1) + C[13]*eps(2,2) + 2.0 * ( C[16]*eps(1,2) + C[18]*eps(2,0) + C[19]*eps(0,1));
	sig(0,1) = C[ 5]*eps(0,0) + C[10]*eps(1,1) + C[14]*eps(2,2) + 2.0 * ( C[17]*eps(1,2) + C[19]*eps(2,0) + C[20]*eps(0,1));

	sig(2,1) = sig(1,2);
	sig(0,2) = sig(2,0);
	sig(1,0) = sig(0,1);

	return sig;
}
Set::Vector
CrystalPlastic::operator () (std::array<Set::Matrix,3> &gradgradu)
{
	std::array<Set::Matrix,3> gradeps;
	gradeps[0](0,0) = gradgradu[0](0,0);
	gradeps[0](0,1) = gradgradu[0](0,1);
	gradeps[0](0,2) = gradgradu[0](0,2);
	gradeps[0](1,0) = 0.5*(gradgradu[0](1,0) + gradgradu[1](0,0));
	gradeps[0](1,1) = 0.5*(gradgradu[0](1,1) + gradgradu[1](0,1));
	gradeps[0](1,2) = 0.5*(gradgradu[0](1,2) + gradgradu[1](0,2));
	gradeps[0](2,0) = 0.5*(gradgradu[0](2,0) + gradgradu[2](0,0));
	gradeps[0](2,1) = 0.5*(gradgradu[0](2,1) + gradgradu[2](0,1));
	gradeps[0](2,2) = 0.5*(gradgradu[0](2,2) + gradgradu[2](0,2));
	gradeps[1](0,0) = 0.5*(gradgradu[1](0,0) + gradgradu[0](1,0));
	gradeps[1](0,1) = 0.5*(gradgradu[1](0,1) + gradgradu[0](1,1));
	gradeps[1](0,2) = 0.5*(gradgradu[1](0,2) + gradgradu[0](1,2));
	gradeps[1](1,0) = gradgradu[1](1,0);
	gradeps[1](1,1) = gradgradu[1](1,1);
	gradeps[1](1,2) = gradgradu[1](1,2);
	gradeps[1](2,0) = 0.5*(gradgradu[1](2,0) + gradgradu[2](1,0));
	gradeps[1](2,1) = 0.5*(gradgradu[1](2,1) + gradgradu[2](1,1));
	gradeps[1](2,2) = 0.5*(gradgradu[1](2,2) + gradgradu[2](1,2));
	gradeps[2](0,0) = 0.5*(gradgradu[2](0,0) + gradgradu[0](2,0));
	gradeps[2](0,1) = 0.5*(gradgradu[2](0,1) + gradgradu[0](2,1));
	gradeps[2](0,2) = 0.5*(gradgradu[2](0,2) + gradgradu[0](2,2));
	gradeps[2](1,0) = 0.5*(gradgradu[2](1,0) + gradgradu[1](2,0));
	gradeps[2](1,1) = 0.5*(gradgradu[2](1,1) + gradgradu[1](2,1));
	gradeps[2](1,2) = 0.5*(gradgradu[2](1,2) + gradgradu[1](2,2));
	gradeps[2](2,0) = gradgradu[2](2,0);
	gradeps[2](2,1) = gradgradu[2](2,1);
	gradeps[2](2,2) = gradgradu[2](2,2);

	Set::Vector f = Set::Vector::Zero();

	f(0) =  C[ 0]*gradeps[0](0,0) + C[ 1]*gradeps[1](1,0) + C[ 2]*gradeps[2](2,0) + 2.0 * ( C[ 3]*gradeps[1](2,0) + C[ 4]*gradeps[2](0,0) + C[ 5]*gradeps[0](1,0)) + 
		C[ 5]*gradeps[0](0,1) + C[10]*gradeps[1](1,1) + C[14]*gradeps[2](2,1) + 2.0 * ( C[17]*gradeps[1](2,1) + C[19]*gradeps[2](0,1) + C[20]*gradeps[0](1,1)) + 
		C[ 4]*gradeps[0](0,2) + C[ 9]*gradeps[1](1,2) + C[13]*gradeps[2](2,2) + 2.0 * ( C[16]*gradeps[1](2,2) + C[18]*gradeps[2](0,2) + C[19]*gradeps[0](1,2));

	f(1) =  C[ 5]*gradeps[0](0,0) + C[10]*gradeps[1](1,0) + C[14]*gradeps[2](2,0) + 2.0 * ( C[17]*gradeps[1](2,0) + C[19]*gradeps[2](0,0) + C[20]*gradeps[0](1,0)) +
		C[ 1]*gradeps[0](0,1) + C[ 6]*gradeps[1](1,1) + C[ 7]*gradeps[2](2,1) + 2.0 * ( C[ 8]*gradeps[1](2,1) + C[ 9]*gradeps[2](0,1) + C[10]*gradeps[0](1,1)) +
		C[ 3]*gradeps[0](0,2) + C[ 8]*gradeps[1](1,2) + C[12]*gradeps[2](2,2) + 2.0 * ( C[15]*gradeps[1](2,2) + C[16]*gradeps[2](0,2) + C[17]*gradeps[0](1,2));

	f(2) =  C[ 4]*gradeps[0](0,0) + C[ 9]*gradeps[1](1,0) + C[13]*gradeps[2](2,0) + 2.0 * ( C[16]*gradeps[1](2,0) + C[18]*gradeps[2](0,0) + C[19]*gradeps[0](1,0)) +
		C[ 3]*gradeps[0](0,1) + C[ 8]*gradeps[1](1,1) + C[12]*gradeps[2](2,1) + 2.0 * ( C[15]*gradeps[1](2,1) + C[16]*gradeps[2](0,1) + C[17]*gradeps[0](1,1)) + 
		C[ 2]*gradeps[0](0,2) + C[ 7]*gradeps[1](1,2) + C[11]*gradeps[2](2,2) + 2.0 * ( C[12]*gradeps[1](2,2) + C[13]*gradeps[2](0,2) + C[14]*gradeps[0](1,2));
	return f;
}
#endif

}
}
}
