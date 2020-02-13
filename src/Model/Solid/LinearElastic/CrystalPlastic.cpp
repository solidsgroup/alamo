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
		//if(Time>.000001) 
		//{
		//Util::Message(INFO,"a = ", a, "\n");
		//Util::Message(INFO,"sig = ", sig);
		//}
		
		if(a > slipSystem[i].Tcrss) slipSystem[i].on = true;		
		else slipSystem[i].on = false;
	}	
}

void CrystalPlastic::AdvanceEsp( const Set::Matrix& sig)
{
	//GetActivePlains(sig);
	Set::Matrix Dp = Set::Matrix::Zero();
	
	for(int i = 0; i < 12; i++)
	{ 
		double a = slipSystem[i].s.transpose()*sig*slipSystem[i].n;
		//if(Time >.00035) Util::Message(INFO,"tcrss = ", slipSystem[i].Tcrss," a = ", a, "\n sig = ",sig);
		if(abs(a) > slipSystem[i].Tcrss)
		{
			//double a = slipSystem[i].s.transpose()*sig*slipSystem[i].n;
			int sign = sgn(a);
			slipSystem[i].galphad = /* (double)sign* */gammadot0*pow( abs(a)/slipSystem[i].Tcrss, n);
			Dp += slipSystem[i].galphad*(double)sign*slipSystem[i].s*slipSystem[i].n.transpose();
			//Dp += (slipSystem[i].s*slipSystem[i].n.transpose() + slipSystem[i].n*slipSystem[i].s.transpose())*slipSystem[i].galphad/2;
			slipSystem[i].on = true;
		}
		else slipSystem[i].on = false;
	}
	esp = esp + Dp*dt;
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
		if(slipSystem[i].on)
		{
			slipSystem[i].gam = slipSystem[i].gam + slipSystem[i].galphad * dt; 
		}
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
		if(slipSystem[a].on)
		{
			slipSystem[a].Tcrss = slipSystem[a].Tcrss + temp*dt;
		}
	}
	//Util::Message(INFO, haa.transpose() ,"\n");
	//Util::Message(INFO, H ,"\n");
}
void CrystalPlastic::update(const Set::Matrix es, Set::Matrix& sigma)
{
	Time += dt;
	AdvanceEsp(sigma);
	//if(Time >.00035) Util::Message(INFO,"esp = ", esp);
	sigma = UpdateSigma(es);
	LatentHardening();
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
	return slipSystem[index].gam;
}
void CrystalPlastic::Setdt(double _dt)
{
	dt = _dt;
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
	ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();

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
	for(int p = 0; p < AMREX_SPACEDIM; p++) 
		for(int q = 0; q < AMREX_SPACEDIM; q++) 
			for(int s = 0; s < AMREX_SPACEDIM; s++) 
				for(int t = 0; t < AMREX_SPACEDIM; t++)
				{
					ddw(p,q,s,t) = 0.0;
					for(int i = 0; i < 3; i++) 
						for(int j = 0; j < 3; j++) 
							for(int k = 0; k < 3; k++) 
								for(int l = 0; l < 3; l++) 
									ddw(p,q,s,t) += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
}
/*
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
*/
}
}
}
}