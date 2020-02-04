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
	initializeSlip(R);
	F = setF(); 
	define(C11, C12, C44, R);
}
CrystalPlastic::CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	initializeSlip(phi1, Phi, phi2);
	F = setF(); 
	define(C11, C12, C44, phi1, Phi, phi2);
}
void CrystalPlastic::initializeSlip(Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
	Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
	initializeSlip(m);
}
Eigen::Matrix<double,12,12> CrystalPlastic::setF()
{
	Eigen::Matrix<double,12,12> f = Eigen::Matrix<double,12,12>::Zero();
	Eigen::Matrix<double,12,12> ggg = Eigen::Matrix<double,12,12>::Zero();
	f(0,1) = cc; f(1,2) = cc; f(2,3) = hh; f(3,4) = cc; f(4,5) = cc; f(5,6) = hh; f(6,7) = cc; f(7,8) = cc; f(8,9) = gg; f(9,10) = cc; f(10,11) = cc;
	f(0,2) = cc; f(1,3) = gg; f(2,4) = gg; f(3,5) = cc; f(4,6) = ss; f(5,7) = gg; f(6,8) = cc; f(7,9) = ss; f(8,10) = nn; f(9,11) = cc;
	f(0,3) = ss; f(1,4) = nn; f(2,5) = ss; f(3,6) = gg; f(4,7) = gg; f(5,8) = ss; f(6,9) = hh; f(7,10) = gg; f(8,11) = gg;
	f(0,4) = gg; f(1,5) = gg; f(2,6) = gg; f(3,7) = nn; f(4,8) = hh; f(5,9) = nn; f(6,10) = gg; f(7,11) = hh;
	f(0,5) = hh; f(1,6) = gg; f(2,7) = hh; f(3,8) = gg; f(4,9) = gg; f(5,10) = gg; f(6,11) = ss;
	f(0,6) = nn; f(1,7) = ss; f(2,8) = ss; f(3,9) = gg; f(4,10) = hh; f(5,11) = gg;
	f(0,7) = gg; f(1,8) = hh; f(2,9) = gg; f(3,10) = ss; f(4,11) = ss;
	f(0,8) = gg; f(1,9) = ss; f(2,10) = gg; f(3,11) = hh;
	f(0,9) = hh; f(1,10) = hh; f(2,11) = nn;
	f(0,10) = ss; f(1,11) = gg;
	f(0,11) = gg;
	ggg = f.transpose() + f;
	//Util::Message(INFO,g);
	return ggg;
}
void CrystalPlastic::initializeSlip(Set::Matrix R)
{
	double n = 1/sqrt(3);
	double s = 1/sqrt(2);
	Set::Vector n1 = {n,n,-n}; //{1,1,1};
	Set::Vector n2 = {n,-n,n}; //{-1,-1,1};
	Set::Vector n3 = {n,-n,-n}; //{-1,1,1};
	Set::Vector n4 = {n,n,n}; //{1,-1,1};
	
	Set::Vector s11 = {s,0,s}; //{0,-1,1};
	Set::Vector s12 = {0,s,s}; //{1,0,-1};
	Set::Vector s13 = {s,-s,0}; //{-1,1,0};
	Set::Vector s21 = {s,s,0}; //{0,1,1};
	Set::Vector s22 = {0,s,s}; //{-1,0,-1};
	Set::Vector s23 = {s,0,-s}; //{1,-1,0};
	Set::Vector s31 = {s,0,s}; //{0,-1,1};
	Set::Vector s32 = {s,s,0}; //{-1,0,-1};
	Set::Vector s33 = {0,s,-s}; //{1,1,0};
	Set::Vector s41 = {s,0,-s}; //{0,1,1};
	Set::Vector s42 = {0,-s,s}; //{1,0,-1};
	Set::Vector s43 = {s,-s,0}; //{-1,-1,0};

	this->slp1.n = R*n1; this->slp1.s = R*s11;
	this->slp2.n = R*n1; this->slp2.s = R*s12;
	this->slp3.n = R*n1; this->slp3.s = R*s13;

	this->slp4.n = R*n2; this->slp4.s = R*s21;
	this->slp5.n = R*n2; this->slp5.s = R*s22;
	this->slp6.n = R*n2; this->slp6.s = R*s23;

	this->slp7.n = R*n3; this->slp7.s = R*s31;
	this->slp8.n = R*n3; this->slp8.s = R*s32;
	this->slp9.n = R*n3; this->slp9.s = R*s33;

	this->slp10.n = R*n4; this->slp10.s = R*s41;
	this->slp11.n = R*n4; this->slp11.s = R*s42;
	this->slp12.n = R*n4; this->slp12.s = R*s43;

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

	//Util::Message(INFO,"R = ", R );

	for(int i = 0; i < 12; i++)
	{
		slipSystem[i].on = false;
		slipSystem[i].galphad = gammadot0;
		slipSystem[i].Tcrss = tcrss;
		slipSystem[i].gam = 0;
		//Util::Message(INFO,"n = ", slipSystem[i].n.transpose());
		//Util::Message(INFO,"s = ", slipSystem[i].s.transpose());
		//Util::Message(INFO,"galpha = ", slipSystem[i].galpha);
		//Util::Message(INFO,"Tcrss = ", slipSystem[i].Tcrss);
	}
}

double CrystalPlastic::CalcSSigN (const Set::Vector ss, const Set::Vector nn, const Set::Matrix& sig) 
{
	double a;
	a = ss.transpose()*sig*nn;
	return a;
}
std::array<double, 12> CrystalPlastic::StressSlipSystem(const Set::Matrix& sig)
{
	std::array<double,12> a;
	for(int i = 0; i < 12; i++)
	{
		a[i] = slipSystem[i].s.transpose()*sig*slipSystem[i].n;
	}
	return a;
}
void CrystalPlastic::GetActivePlains(const Set::Matrix& sig)
{
	double a = 0;
	for(int i = 0; i < 12; i++)
	{
		a = abs(slipSystem[i].s.transpose()*sig*slipSystem[i].n);
		//if(Time>.0035) Util::Message(INFO,"a = ", a);
		//Util::Message(INFO,"sig = ", sig);
		
		if(a > slipSystem[i].Tcrss) slipSystem[i].on = true;		
		else slipSystem[i].on = false;
	}	
}

void CrystalPlastic::AdvanceEsp( const Set::Matrix& sig)
{
	GetActivePlains(sig);
	Set::Matrix temp = Set::Matrix::Zero();
	
	for(int i = 0; i < 12; i++)
	{ 
		//double a = slipSystem[i].s.transpose()*sig*slipSystem[i].n;
		//if(Time >.0035) Util::Message(INFO,"tcrss = ", slipSystem[i].Tcrss," a = ", a, "\n sig = ",sig);
		if(slipSystem[i].on)
		{
			double a = slipSystem[i].s.transpose()*sig*slipSystem[i].n;
			int sign = sgn(a);
			slipSystem[i].galphad = (double)sign*gammadot0*pow( abs(a)/slipSystem[i].Tcrss, n);
			
			temp += slipSystem[i].galphad*(double)sign*slipSystem[i].s*slipSystem[i].n.transpose();
		}
	}
	// Euler integration
	esp = esp + temp*dt;
}
Eigen::Matrix<double,12,1> CrystalPlastic::G()
{
	Eigen::Matrix<double,12,1> k = Eigen::Matrix<double,12,1>::Zero();
	double l;
	for(int a = 0; a < 12; a++)
	{
		l =  hs + (h0 - hs) * 1/pow(cosh( ((h0 - hs)/( ts - t0)) * slipSystem[a].gam ) , 2);
		for(int b = 0; b < 12; b++)
		{
			if (b == a) continue;
			k(a) += F(a,b)*tanh(slipSystem[b].gam/gammadot0);
		}
		k(a) += 1.0;
		k(a) = k(a)*l;
	}
	//Util::Message(INFO, k,"\n");
	return k;
}

void CrystalPlastic::LatentHardening()
{
	Eigen::Matrix<double,12,12> H = Eigen::Matrix<double,12,12>::Identity();
	Eigen::Matrix<double,12,1> haa = Eigen::Matrix<double,12,1>::Zero();
	for(int i = 0; i < 12; i++) 
	{
		slipSystem[i].gam = slipSystem[i].gam + slipSystem[i].galphad * dt; 
	}
	haa = G();
	for(int a = 0; a < 12; a++)
	{
		double temp = 0;
		for(int b = 0; b < 12; b++)
		{
			if(a == b) H(a,b) = haa(b);
			else H(a,b) = q * haa(b);
			temp += H(a,b) * abs(slipSystem[b].galphad);
		}
		slipSystem[a].Tcrss = slipSystem[a].Tcrss + temp*dt;
	}
	//Util::Message(INFO, haa.transpose() ,"\n");
	//Util::Message(INFO, H ,"\n");
}

void CrystalPlastic::update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt)
{
	//for(double t = 0.0; t < _dt; t += dt)
	//{
		Time += dt;
		AdvanceEsp(sigma);
		LatentHardening();
		sigma = UpdateSigma(es);
	//}
}
Set::Matrix CrystalPlastic::UpdateSigma(const Set::Matrix& es)
{
	Set::Matrix temp = (es - esp);
	return operator()(temp); 
}
Set::Matrix CrystalPlastic::GetEsp() const
{
	return esp;
}
double CrystalPlastic::getGamma(int index) const
{
	//return slipSystem[index].galpha;
}
void CrystalPlastic::Setdt(double _dt)
{
	dt = _dt;
}
//----------DFP Functions-----------//

Eigen::Matrix<amrex::Real,9,1> CrystalPlastic::DFP(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig)
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
		double temp = b - 0.8*((df1*(b - a)) / (df1 - df2));

		a = b;
		b = temp;
	}
	return b;
}

Set::Scalar CrystalPlastic::f(vector2d x, const Set::Matrix& es)
{
	Set::Matrix epsilon = es; 
	for(int i= 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(Mask(i,j)) epsilon(i,j) = es(i,j); 
			else epsilon(i,j) = x(i*3+j);
		}
	}

	//					 epsilon(0,1) = x(0); epsilon(0,2) = x(1);
	//epsilon(1,0) = x(2); epsilon(1,1) = x(3); epsilon(1,2) = x(4);
	//epsilon(2,0) = x(5); epsilon(2,1) = x(6); epsilon(2,2) = x(7);
	Set::Matrix temp = epsilon - esp;
	Set::Scalar s = W(temp);
	//Util::Message(INFO,"W = ",s);
	return s;
}

Eigen::Matrix<amrex::Real,9,1> CrystalPlastic::getGrad(vector2d x, double dx, const Set::Matrix& sig)
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
//--------------------------------//
Set::Matrix CrystalPlastic::relax(const Set::Matrix& _es, const double e, const Set::iMatrix& _mask) 
{
	Mask = _mask;
	vector2d x;
	for(int i= 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(Mask(i,j)) x(i*3+j) = _es(i,j); 
			else x(i*3+j) = 0;
		}
	}
	x(4) = -e-1e-2; x(8) = -e-1e-2; //for extension 
	//x(3) = e-1e-2; //x(7) = -e-1e-2; // for shear
	vector2d xprev = vector2d::Zero(); xprev(1) = 100;
	Set::Matrix ep = Set::Matrix::Zero();

	while(abs(x.norm() - xprev.norm()) >= 1e-5)
	{
		xprev = x;
		x = DFP(x,1e-5,0.1,0.8,1e-5,_es);
		//Util::Message(INFO,"ffff= ", counter);
	}
	for(int i= 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(Mask(i,j)) ep(i,j) = _es(i,j); 
			else ep(i,j) = x(i*3+j);
		}
	}
	//ep(0,0) = e;	ep(0,1) = x(0); ep(0,2) = x(1);
	//ep(1,0) = x(2); ep(1,1) = x(3); ep(1,2) = x(4);
	//(2,0) = x(5); ep(2,1) = x(6); ep(2,2) = x(7);
	return ep;
}
void
CrystalPlastic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
	define(C11,C12,C44,m);
}
void
CrystalPlastic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
{
  
	amrex::Real Ctmp[3][3][3][3];
	amrex::Real Crot[3][3][3][3];
	C = Set::Matrix4<3,Set::Sym::MajorMinor>::Zero();

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
					C(p,q,s,t) = 0.0;
					for(int i = 0; i < 3; i++) 
						for(int j = 0; j < 3; j++) 
							for(int k = 0; k < 3; k++) 
								for(int l = 0; l < 3; l++) 
									C(p,q,s,t) += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
				}
}
void
CrystalPlastic::Randomize()
{
	Set::Scalar C11 = 0.5 + 0.5*Util::Random();
	Set::Scalar C12 = 0.5 + 0.5*Util::Random();
	Set::Scalar C44 = 0.5 + 0.5*Util::Random();

	Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();

	define(C11,C12,C44,phi1,Phi,phi2);
}
Set::Scalar 
CrystalPlastic::W(Set::Matrix &gradu) const
{
	Set::Matrix sig = C*gradu;
	return 0.5 * (sig*gradu).trace();
}
Set::Matrix CrystalPlastic::operator () (Set::Matrix &gradu,bool) const
{
	return C*gradu;
}
Set::Matrix
CrystalPlastic::DW (Set::Matrix &gradu) const
{
	return (*this)(gradu);
}
Set::Vector
CrystalPlastic::operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu,bool)
{
	Set::Vector ret = Set::Vector::Zero();
	for (int i = 0; i < AMREX_SPACEDIM; i++)
		for (int j = 0; j < AMREX_SPACEDIM; j++)
			for (int k = 0; k < AMREX_SPACEDIM; k++)
				for (int l = 0; l < AMREX_SPACEDIM; l++)
					ret(i) += C(i,j,k,l)*gradgradu[k](l,j);
	return ret;
}
Set::Vector
CrystalPlastic::DW (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu)
{
	return (*this)(gradgradu);
}
Set::Matrix4<3,Set::Sym::MajorMinor>
CrystalPlastic::DDW(Set::Matrix &/*gradu*/) const
{
	return C;
}
}
}
}