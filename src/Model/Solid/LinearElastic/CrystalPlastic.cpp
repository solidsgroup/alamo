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
	define(C11, C12, C44, R);
}
CrystalPlastic::CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	initializeSlip(phi1, Phi, phi2);
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
void CrystalPlastic::initializeSlip(Set::Matrix R)
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
		slipSystem[i].Tcrss = 0.5;
		Util::Message(INFO,"n = ", slipSystem[i].n.transpose());
		Util::Message(INFO,"s = ", slipSystem[i].s.transpose());
		Util::Message(INFO,"galpha = ", slipSystem[i].galpha);
		Util::Message(INFO,"Tcrss = ", slipSystem[i].Tcrss);
	}
}

double CrystalPlastic::CalcSSigN (const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig) 
{
	double a;
	a = ss.transpose()*sig*nn;
	return a;
}
void CrystalPlastic::GetActivePlains(const Set::Matrix& sig)
{
	double a = 0;
	for(int i = 0; i <= 12; i++)
	{
		//calculate a = abs(s*sigma*n)
		a = abs(CalcSSigN(slipSystem[i].s,slipSystem[i].n, sig));
		
		if(a > slipSystem[i].Tcrss) slipSystem[i].on = true;		
		else slipSystem[i].on = false;
	}	
}

Set::Scalar CrystalPlastic::GetGammaDot(const Set::Vector ss, const Set::Vector nn, const Set::Matrix& sig)
{
	//double gamma;
	//double power = pow( (abs(CalcSSigN( ss,nn,sig )) /slipSystem[i].Tcrss ) ,n );
	//gamma = gammadot0*power;
	//return gamma;
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
			double gammadot = pow( abs(a)/slipSystem[i].Tcrss, n);
			slipSystem[i].galpha = gammadot;
			int sign = sgn(a);
			temp += gammadot*sign*slipSystem[i].s*slipSystem[i].n.transpose();
		}
	}
	// Euler integration
	esp = esp + temp*dt;
}
void CrystalPlastic::LatentHardening()
{
	std::array<double,12> hab;
	Set::Matrix H = Set::Matrix::Zero();
	for(int i = 0; i < 12; i++)
	{
		double g = 0; 
		for(auto f : slipSystem) 
		{
			g += f.galpha;
		}
		g = g*Time;
		if(coplanar(slipSystem[i].n, i))
		{
			hab[i] = (hs + (hs - h0) * pow( 1/cosh( (h0-hs)*g / (ts-t0)) ,2));
		}
		else
		{
			hab[i] = q*(hs + (hs - h0) * pow( 1/cosh( (h0-hs)*g / (ts-t0) ) ,2));
		}
		slipSystem[i].Tcrss = hab[i]*abs(slipSystem[i].galpha);
		
		if(slipSystem[i].Tcrss == 0)
		{
			slipSystem[i].Tcrss = 0.1;
		}
		Util::Message(INFO,"Tcrss = ", slipSystem[i].Tcrss);
		
	}
	
}
bool CrystalPlastic::coplanar(Set::Vector nn, int i) const
{
	Set::Vector a = slipSystem[i].n.cross(nn);
	//Util::Message(INFO,"a = ", a.norm());
	if(a.norm() <= 1e-8)
	{
		return true;
	}
	else 
	{
		return false;
	}
}
void CrystalPlastic::update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt)
{
	for(double t = 0.0; t < _dt; t += dt)
	{
		Time += dt;
		//LatentHardening();
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
double CrystalPlastic::getGamma(int index) const
{
	return slipSystem[index].galpha;
}
void CrystalPlastic::Setdt(double _dt)
{
	dt = _dt;
}
//----------DFP Functions-----------//

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
		
		
		double temp = secantMethod(dx, alpha1, alpha2, tol, xnew, R, sig);
		//Util::Message(INFO,"temp = ", temp);
		xprev = xnew;
		xnew = xnew + temp * R;

		vector2d gamma = getGrad(xnew, dx, sig) - getGrad(xprev,dx, sig);
		vector2d del = xnew - xprev;

		Hprev = Hnew;
		Hnew = Hprev + (del * del.transpose()) / (del.transpose() * del) - (Hprev * (gamma*gamma.transpose()) * Hprev) / (gamma.transpose() * Hprev * gamma);
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
	Set::Matrix epsilon = es; epsilon(1,1) = x(0); epsilon(2,2) = x(1);
	Set::Matrix temp = epsilon - esp;
	Set::Scalar s = W(temp);
	//Util::Message(INFO,"W = ",s);
	return s;
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
	//Util::Message(INFO,"grad = ", grad);
	return grad;
}
//--------------------------------//
Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> CrystalPlastic::relax(const Set::Matrix& sig, const double e) 
{
	
	vector2d x; x(0) = -e+1; x(1) = -e + 1; //w1 = 10; w2 = 10;
	vector2d xprev = vector2d::Zero(); xprev(1) = 100;
	int counter = 1;
	while(abs(x.norm() - xprev.norm()) >= 1e-5)
	{
		xprev = x;
		x = DFP(x,1e-5,0.1,0.8,1e-5,sig);
		if(counter % 1 == 0)
		{
			//w1 = w1 * 1.1;
			//w2 = w2 * 1.1;
		}
		counter++;
		//Util::Message(INFO,"ffff= ", counter);
	}
	return x;
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