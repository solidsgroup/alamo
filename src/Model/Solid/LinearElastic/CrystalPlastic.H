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
private:
	static constexpr int dim = AMREX_SPACEDIM-1;
	using vector2d = Eigen::Matrix<amrex::Real,dim,1>;
	using matrix22 = Eigen::Matrix<amrex::Real,dim,dim>;

	void initializeSlip();
	
	struct slip
	{
	Set::Vector n;
	Set::Vector s;
	bool on;
	};
	
	std::array<Set::Scalar,21> C;
	//********Matrial Slip Data**********//

	std::array<slip,12> slipSystem;
	slip slp1; slip slp2;
	slip slp3; slip slp4;
	slip slp5; slip slp6; 
	slip slp7; slip slp8;
	slip slp9; slip slp10;
	slip slp11; slip slp12;
	//********Matrial Matrix**********//
	
	Set::Matrix esp = Set::Matrix::Zero();
	//********Matrial Constants**********//

	static Set::Scalar constexpr Tcrss = 0.5;
	static Set::Scalar constexpr n = 2.0;
	static Set::Scalar constexpr gammadot0 = 0.5;
	Set::Scalar dt = 1e-5;
	//***********DFP************//

	double w1 = 10.0;
	double w2 = 10.0;

public:
	CrystalPlastic() 
	{
		initializeSlip();
	} 
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	~CrystalPlastic() {};
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
	double CalcSSigN (const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig);
	void GetActivePlains(const Set::Matrix sig); //used in update esp function
	Set::Matrix UpdateSigma(const Set::Matrix es);
	Set::Scalar GetGammaDot(const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig);
	void update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt);
	void AdvanceEsp(const Set::Matrix sig);
	void Setdt(const double _dt);
	vector2d reflux(const Set::Matrix& sig, double e);
	Set::Matrix GetEsp() const;
//DFP
	double secantMethod(double dx, double a1, double a2, double tol, vector2d x_new, vector2d r, const Set::Matrix& sig);
	Set::Scalar f(vector2d x, const Set::Matrix& es);
	vector2d getGrad(vector2d x, double dx, const Set::Matrix& sig);
	vector2d DFP(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig);

	matrix22 Jacobi(vector2d x, const Set::Matrix& es, double dx);
	vector2d newtonrap(double dx, vector2d x0, const Set::Matrix& es);
	vector2d newtonfunc(vector2d x, const Set::Matrix& es);
//
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


};

}
}
}

#endif 