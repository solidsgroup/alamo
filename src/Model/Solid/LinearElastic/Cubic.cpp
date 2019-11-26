#include "Util/Util.H"
#include "Cubic.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{

Cubic::Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
{
	define(C11, C12, C44, R);
}
Cubic::Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	define(C11, C12, C44, phi1, Phi, phi2);
}

void
Cubic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
{
	Eigen::Matrix3d m;
	m =     Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
	define(C11,C12,C44,m);
}
void
Cubic::define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
{
  
	amrex::Real Ctmp[3][3][3][3];
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
Cubic::Randomize()
{
	Set::Scalar C11 = 0.5 + 0.5*Util::Random();
	Set::Scalar C12 = 0.5 + 0.5*Util::Random();
	Set::Scalar C44 = 0.5 + 0.5*Util::Random();
	Randomize(C11,C12,C44);
}
void
Cubic::Randomize(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
{
	Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
	define(C11,C12,C44,phi1,Phi,phi2);
}

Set::Scalar 
Cubic::W(Set::Matrix &gradu) const
{
	Set::Matrix sig = C*gradu;
	return 0.5 * (sig*gradu).trace();
}

Set::Matrix
Cubic::operator () (Set::Matrix &gradu,bool) const
{
	return C*gradu;
}
Set::Matrix
Cubic::DW (Set::Matrix &gradu) const
{
	return (*this)(gradu);
}

Set::Vector
Cubic::operator () (Set::Matrix3 &gradgradu,bool)
{
	Set::Vector ret = Set::Vector::Zero();
	for (int i = 0; i < AMREX_SPACEDIM; i++)
		for (int j = 0; j < AMREX_SPACEDIM; j++)
			for (int k = 0; k < AMREX_SPACEDIM; k++)
				for (int l = 0; l < AMREX_SPACEDIM; l++)
					ret(i) += C(i,j,k,l)*gradgradu(k,l,j);
	return ret;
}
Set::Vector
Cubic::DW (Set::Matrix3 &gradgradu)
{
	return (*this)(gradgradu);
}
Set::Matrix4<3,Set::Sym::MajorMinor>
Cubic::DDW(Set::Matrix &/*gradu*/) const
{
	return C;
}

}
}
}
