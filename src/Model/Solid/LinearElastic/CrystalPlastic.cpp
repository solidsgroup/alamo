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
	//Q = coplanar();
	F = setF(); 
	define(C11, C12, C44, R);
}
CrystalPlastic::CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	initializeSlip(phi1, Phi, phi2);
	//Q = coplanar();
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
	for(int i = 0; i < 12; i++)
	{
		slipSystem[i].on = false;
		slipSystem[i].galpha = 0;
		slipSystem[i].Tcrss = tcrss;
		//Util::Message(INFO,"n = ", slipSystem[i].n.transpose());
		//Util::Message(INFO,"s = ", slipSystem[i].s.transpose());
		//Util::Message(INFO,"galpha = ", slipSystem[i].galpha);
		//Util::Message(INFO,"Tcrss = ", slipSystem[i].Tcrss);
	}
}
Eigen::Matrix<double,12,12> CrystalPlastic::coplanar()
{
	Eigen::Matrix<double,12,12> D;
	for(int i = 0; i < 12; i++)
	{
		Set::Vector h = slipSystem[i].n;
		for(int j = 0; j < 12; j++)
		{
			Set::Vector k = slipSystem[j].n;
			Set::Vector a = h.cross(k);

			if(a.norm() <= 1e-8) D(i,j) = 1.0;
			else D(i,j) = q;
		}
	}
	//Util::Message(INFO, D);
	return D;
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
		a[i] = CalcSSigN(slipSystem[i].s, slipSystem[i].n, sig);
	}
	return a;
}
void CrystalPlastic::GetActivePlains(const Set::Matrix& sig)
{
	double a = 0;
	for(int i = 0; i < 12; i++)
	{
		a = abs(CalcSSigN(slipSystem[i].s,slipSystem[i].n, sig));
		//if(Time > 0.01) Util::Message(INFO, a ,"\n");
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
		if(slipSystem[i].on)
		{
			double a = CalcSSigN(slipSystem[i].s,slipSystem[i].n,sig);
			int sign = sgn(a);
			slipSystem[i].galpha = (double)sign*gammadot0*pow( abs(a)/slipSystem[i].Tcrss, n);
			//if(Time > 0.01) Util::Message(INFO, sgn(a));
			temp += slipSystem[i].galpha*(double)sign*slipSystem[i].s*slipSystem[i].n.transpose();
		}
	}
	// Euler integration
	esp = esp + temp*dt;
}
Eigen::Matrix<double,12,1> CrystalPlastic::G()
{
	Eigen::Matrix<double,12,1> a = Eigen::Matrix<double,12,1>::Ones();;
	for(int i = 0; i < 12; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			if (i == j) continue;
			a(i) += F(i,j)*tanh(gam/gammadot0);
		}
		a(i) += 1.0;
	}
	//Util::Message(INFO, a,"\n");
	return a;
}
void CrystalPlastic::LatentHardening()
{
	Eigen::Matrix<double,12,12> H = Eigen::Matrix<double,12,12>::Identity();
	Eigen::Matrix<double,12,1> haa = Eigen::Matrix<double,12,1>::Zero();
	double gammatemp = 0;
	for(int i = 0; i < 12; i++) 
	{
		gammatemp += abs(slipSystem[i].galpha);
	}
	gam = gammatemp;// gam + gammatemp*dt;
	/*
	double h = hs + (h0 - hs) * 1/pow(cosh( ((h0 - hs)/( ts - t0)) * gam ) , 2);
	H = Q*h;
	//if(Time > 0.0028) Util::Message(INFO, h,"\n");
	for(int a = 0; a < 12; a++)
	{
		double temp = 0;
		for(int b = 0; b < 12; b++)
		{
			temp += H(a,b) * abs(slipSystem[b].galpha);
		}
		slipSystem[a].Tcrss = slipSystem[a].Tcrss + temp*dt;
		//if(Time > 0.0034) Util::Message(INFO, slipSystem[a].Tcrss ,"\n");
		if(slipSystem[a].Tcrss < tcrss)
		{
			//slipSystem[a].Tcrss = tcrss;
		}
	}
	*/
/*
	for(int i = 0; i < 12; i++)
	{
		double G = 0; 
		double a =  hs + (h0 - hs) * 1/pow(cosh( ((h0 - hs)/( ts - t0)) * slipSystem[i].galpha ) , 2);
		for(int j = 0; j < 12; j++)
		{
			if(i == j) continue;
			G = F(i,j) * tanh(slipSystem[j].galpha/gammadot0);
		}
		G += 1;
		haa(i) = a*G;
	}
*/

	haa = G();
	double a =  hs + (h0 - hs) * 1/pow(cosh( ((h0 - hs)/( ts - t0)) * gam ) , 2);
	haa = a* haa;
	for(int i = 0; i < 12; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			if(i == j) H(i,j) = haa(i);

			else H(i,j) = q * haa(i);
		}
	}
	for(int a = 0; a < 12; a++)
	{
		double temp = 0;
		for(int b = 0; b < 12; b++)
		{
			temp += H(a,b) * abs(slipSystem[b].galpha);
		}
		//if(temp < tcrss) temp = tcrss;
		slipSystem[a].Tcrss = slipSystem[a].Tcrss + temp*dt;
	}
	//Util::Message(INFO, haa.transpose() ,"\n");
	//Util::Message(INFO, H ,"\n");
}

void CrystalPlastic::update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt)
{
	for(double t = 0.0; t < _dt; t += dt)
	{
		Time += dt;
		AdvanceEsp(sigma);
		LatentHardening();
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
double CrystalPlastic::getGamma(int index) const
{
	return slipSystem[index].galpha;
}
void CrystalPlastic::Setdt(double _dt)
{
	dt = _dt;
}
//----------DFP Functions-----------//

Eigen::Matrix<amrex::Real,8,1> CrystalPlastic::DFP(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig)
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
						 epsilon(0,1) = x(0); epsilon(0,2) = x(1);
	epsilon(1,0) = x(2); epsilon(1,1) = x(3); epsilon(1,2) = x(4);
	epsilon(2,0) = x(5); epsilon(2,1) = x(6); epsilon(2,2) = x(7);
	Set::Matrix temp = epsilon - esp;
	Set::Scalar s = W(temp);
	//Util::Message(INFO,"W = ",s);
	return s;
}

Eigen::Matrix<amrex::Real,8,1> CrystalPlastic::getGrad(vector2d x, double dx, const Set::Matrix& sig)
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
Set::Matrix CrystalPlastic::relax(const Set::Matrix& sig, const double e) 
{
	vector2d x = vector2d::Zero(); x(3) = -e-1e-2; x(7) = -e-1e-2;
	vector2d xprev = vector2d::Zero(); xprev(1) = 100;
	Set::Matrix ep = Set::Matrix::Zero();

	while(abs(x.norm() - xprev.norm()) >= 1e-5)
	{
		xprev = x;
		x = DFP(x,1e-5,0.1,0.8,1e-5,sig);
		//Util::Message(INFO,"ffff= ", counter);
	}

	ep(0,0) = e;	ep(0,1) = x(0); ep(0,2) = x(1);
	ep(1,0) = x(2); ep(1,1) = x(3); ep(1,2) = x(4);
	ep(2,0) = x(5); ep(2,1) = x(6); ep(2,2) = x(7);
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

Set::Matrix CrystalPlastic::operator () (Set::Matrix &gradu) const
{
	return C*gradu;
}
Set::Matrix CrystalPlastic::DW (Set::Matrix &gradu) const
{
	return (*this)(gradu);
}

Set::Vector CrystalPlastic::operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu)
{
	Set::Vector ret = Set::Vector::Zero();
	for (int i = 0; i < AMREX_SPACEDIM; i++)
		for (int j = 0; j < AMREX_SPACEDIM; j++)
			for (int k = 0; k < AMREX_SPACEDIM; k++)
				for (int l = 0; l < AMREX_SPACEDIM; l++)
					ret(i) += C(i,j,k,l)*gradgradu[k](l,j);
	return ret;
}
Set::Vector CrystalPlastic::DW (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu)
{
	return (*this)(gradgradu);
}
Set::Matrix4<3,Set::Sym::MajorMinor> CrystalPlastic::DDW(Set::Matrix &gradu) const
{
	return C;
}

}
}
}