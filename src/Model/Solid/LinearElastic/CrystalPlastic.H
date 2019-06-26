#ifndef MODEL_SOLID_LINEARELASTIC_CrystalPlastic_H_
#define MODEL_SOLID_LINEARELASTIC_CrystalPlastic_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"


namespace Model
{
namespace Solid
{
namespace CrystalPlastic
{

class CrystalPlastic : public Model::Solid::LinearElastic::LinearElastic
{
public:
	CrystalPlastic() 
	{
		initializeSlip();
	} 
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	~CrystalPlastic() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps) const;
#if AMREX_SPACEDIM==2
	virtual Set::Vector operator () (std::array<Set::Matrix,2> &gradeps);
#elif AMREX_SPACEDIM==3
	virtual Set::Vector operator () (std::array<Set::Matrix,3> &gradeps);
#endif
//----------------------//
	template <typename T> inline constexpr
	int sgn(T a)
	{
		if (a > 0) return 1;
		if (a < 0) return -1;
		return 0;
	}
	double CalcSSigN (Set::Vector ss, Set::Vector nn);

	void UpdateSigma();
	Set::Scalar GetGammaDot(Set::Vector ss, Set::Vector nn);
	Set::Matrix AdvanceEsp();
	void SetEs(Set::Matrix _es); 
	void Setdt(double _dt);
	Set::Matrix GetSigma();
	Set::Matrix GetEsp();
//-----------------------//
	void Randomize();
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	
	CrystalPlastic operator + (const CrystalPlastic &rhs) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] + rhs.C[i]; return ret;}

	void operator += (const CrystalPlastic &rhs)
	{for (unsigned int i = 0; i < C.size(); i++) C[i] += rhs.C[i];}

	CrystalPlastic operator - (const CrystalPlastic &rhs) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] - rhs.C[i]; return ret;}

	CrystalPlastic operator * (const CrystalPlastic &rhs) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] * rhs.C[i]; return ret;}

	CrystalPlastic operator / (const CrystalPlastic &rhs) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] / rhs.C[i]; return ret;}

	CrystalPlastic operator / (const Set::Scalar alpha) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i]/alpha; return ret;}

	CrystalPlastic operator * (const Set::Scalar alpha) const
	{CrystalPlastic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i]*alpha; return ret;}

private:
	void GetActivePlains(); //used in update esp function
	void initializeSlip();
	struct slip
{
	Set::Vector n;
	Set::Vector s;
};
	std::array<Set::Scalar,21> C;
	//********Matrial Slip Data**********//
	Set::Vector n1 = {1,1,1};
	Set::Vector n2 = {-1,-1,1};
	Set::Vector n3 = {-1,1,1};
	Set::Vector n4 = {1,-1,1};
	
	Set::Vector s11 = {0,-1,1};
	Set::Vector s12 = {1,0,1};
	Set::Vector s13 = {-1,1,0};
	Set::Vector s21 = {0,1,1};
	Set::Vector s22 = {-1,0,-1};
	Set::Vector s23 = {1,-1,0};
	Set::Vector s31 = {0,-1,1};
	Set::Vector s32 = {-1,0,-1};
	Set::Vector s33 = {1,1,0};
	Set::Vector s41 = {0,-1,1};
	Set::Vector s42 = {1,0,-1};
	Set::Vector s43 = {-1,-1,0};

	std::array<slip,12> slipSystem;
	std::array<bool,12> Systems;
	slip slp1;
	slip slp2;
	slip slp3;
	slip slp4;
	slip slp5;
	slip slp6;
	slip slp7;
	slip slp8;
	slip slp9;
	slip slp10;
	slip slp11;
	slip slp12;
	//********Matrial Matrices**********//
	Set::Matrix sigma;
	Set::Matrix esp = Set::Matrix::Zero();
	Set::Matrix es = Set::Matrix::Zero();

	//********Matrial Constants**********//
	Set::Scalar Tcrss = .1;
	Set::Scalar E = 1; 
	Set::Scalar n = 2.1;
	Set::Scalar gammadot0 = .5;
	double gamma0 = 0.1;
	Set::Scalar dt = 1e-5;

	

};

}
}
}

#endif 


