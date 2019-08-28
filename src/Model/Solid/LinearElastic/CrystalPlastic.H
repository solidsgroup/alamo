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
	static constexpr int dim = 8;
	using vector2d = Eigen::Matrix<amrex::Real,dim,1>;
	using matrix22 = Eigen::Matrix<amrex::Real,dim,dim>;

	void initializeSlip(Set::Matrix R = Set::Matrix::Identity());
	void initializeSlip(Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	Eigen::Matrix<double,12,12> setF();
	Eigen::Matrix<double,12,1> G();
	
	struct slip
	{
		Set::Vector n;
		Set::Vector s;
		Set::Scalar Tcrss;
		Set::Scalar galpha;
		bool on;
	};
	//********Matrial Slip Data**********//

	std::array<slip,12> slipSystem;
	slip slp1; slip slp2;
	slip slp3; slip slp4;
	slip slp5; slip slp6; 
	slip slp7; slip slp8;
	slip slp9; slip slp10;
	slip slp11; slip slp12;
	//********Matrial Matrix**********//

	Set::Matrix4<3,Set::Sym::MajorMinor> C;
	Set::Matrix esp = Set::Matrix::Zero();
	Eigen::Matrix<double,12,12> Q;
	Eigen::Matrix<double,12,12> F;
	//********Matrial Constants**********//

	static Set::Scalar constexpr tcrss = 2.0e-3;
	static Set::Scalar constexpr n = 1/0.02; 
	static Set::Scalar constexpr gammadot0 = 0.001; // 0.001 1/s
	static Set::Scalar constexpr t0 = tcrss; 	// MPa
	static Set::Scalar constexpr hs = 1.5*t0; 	// MPa
	static Set::Scalar constexpr h0 = 90*t0; 	// MPa
	static Set::Scalar constexpr ts = 1.3*t0; 	// MPa

	static Set::Scalar constexpr cc = 200.0; 
	static Set::Scalar constexpr gg = 8.0;
	static Set::Scalar constexpr hh = 15.0; 
	static Set::Scalar constexpr ss = 20.0; 
	static Set::Scalar constexpr nn = 0.0;
	
	static Set::Scalar constexpr q = 0; // 1
	Set::Scalar gam = gammadot0;

	Set::Scalar dt = 1.0e-5;
	double Time = 0.0;
	//***********DFP************//
	//double w1 = 10.0;
	//double w2 = 10.0;

public:
	CrystalPlastic() 
	{
		initializeSlip();
	} 
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	~CrystalPlastic() {};

	virtual Set::Matrix operator () (Set::Matrix &eps) const;
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps);
	virtual Set::Scalar    W (Set::Matrix &eps) const;
	virtual Set::Matrix   DW (Set::Matrix &eps) const;
	virtual Set::Vector   DW (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps);
	virtual Set::Matrix4<3,Set::Sym::MajorMinor> DDW (Set::Matrix &eps) const;
//----------------------//
	template <typename T> inline constexpr
	int sgn(T a)
	{
		if (a > 0) return 1;
		if (a < 0) return -1;
		return 0;
	}
	double CalcSSigN (const Set::Vector ss, const Set::Vector nn, const Set::Matrix sig);
	void GetActivePlains(const Set::Matrix& sig); //used in update esp function
	Set::Matrix UpdateSigma(const Set::Matrix es);
	void update(const Set::Matrix es, Set::Matrix& sigma, const Set::Scalar _dt);
	void AdvanceEsp(const Set::Matrix& sig);
	void Setdt(const double _dt);
	void LatentHardening();
	Eigen::Matrix<double,12,12> coplanar();
	Set::Matrix relax(const Set::Matrix& sig, double e);
	Set::Matrix GetEsp() const;
	double getGamma(int index) const;
//DFP
	double secantMethod(double dx, double a1, double a2, double tol, vector2d x_new, vector2d r, const Set::Matrix& sig);
	Set::Scalar f(vector2d x, const Set::Matrix& es);
	vector2d getGrad(vector2d x, double dx, const Set::Matrix& sig);
	vector2d DFP(vector2d x0, double tol, double alpha1, double alpha2, double dx, const Set::Matrix& sig);
//
	void Randomize();
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	
	CrystalPlastic operator + (const CrystalPlastic &rhs) const
	{CrystalPlastic ret = *this; ret.C += rhs.C; return ret;}

	CrystalPlastic operator - (const CrystalPlastic &rhs) const
	{CrystalPlastic ret = *this; ret.C -= rhs.C; return ret;}

	CrystalPlastic operator * (const CrystalPlastic &rhs) const
	{CrystalPlastic ret = *this; ret.C *= rhs.C; return ret;}

	CrystalPlastic operator / (const CrystalPlastic &rhs) const
	{CrystalPlastic ret = *this; ret.C /= rhs.C; return ret;}

	CrystalPlastic operator * (const Set::Scalar alpha) const
	{CrystalPlastic ret = *this; ret.C *= alpha; return ret;}

	CrystalPlastic operator / (const Set::Scalar alpha) const
	{CrystalPlastic ret = *this; ret.C /= alpha; return ret;}

	void operator += (const CrystalPlastic &rhs) {C += rhs.C;}


};

}
}
}

#endif 