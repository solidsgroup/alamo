#ifndef MODEL_SOLID_LINEARELASTIC_CrystalPlastic_H_
#define MODEL_SOLID_LINEARELASTIC_CrystalPlastic_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace CrystalPlastic
{

class CrystalPlastic : public Solid<Set::Sym::MajorMinor>
{
private:

	struct slip
	{
		Set::Vector n;
		Set::Vector s;
		Set::Scalar Tcrss;
		Set::Scalar galphad;
		Set::Scalar gam;
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
	Eigen::Matrix<double,12,12> F;
	//********Matrial Constants**********//

	static Set::Scalar constexpr tcrss = 3.75e-3;//1e-3;
	static Set::Scalar constexpr n = 50.0; 
	static Set::Scalar constexpr gammadot0 = 1e-3; // 0.001 1/s
	static Set::Scalar constexpr t0 = tcrss; 	// MPa
	static Set::Scalar constexpr hs = 1.5*t0; 	// MPa
	static Set::Scalar constexpr h0 = 90*t0; 	// MPa
	static Set::Scalar constexpr ts = 1.3*t0; 	// MPa

	static Set::Scalar constexpr cc = 8.0; 
	static Set::Scalar constexpr gg = 15.0;
	static Set::Scalar constexpr hh = 8.0; 
	static Set::Scalar constexpr ss = 20.0; 
	static Set::Scalar constexpr nn = 8.0;
	
	static Set::Scalar constexpr q = 0; 

	Set::Scalar dt = 1.0e-5;
	double Time = 0.0;
	//******** Functions **********//

	void initializeSlip(Set::Matrix R = Set::Matrix::Identity());
	void initializeSlip(Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	template <typename T> inline constexpr
	int sgn(T a)
	{
		if (a > 0) return 1;
		if (a < 0) return -1;
		return 0;
	}
	Eigen::Matrix<double,12,12> setF();
	Eigen::Matrix<double,12,1> G();
	void GetActivePlains(const Set::Matrix& sig); 
	Set::Matrix UpdateSigma(const Set::Matrix& es);
	void AdvanceEsp(const Set::Matrix& sig);
	void LatentHardening();

public:
	static const KinematicVariable kinvar = KinematicVariable::gradu;
	Set::Matrix Epsp = Set::Matrix::Zero(); // Eshelby
	std::array<Set::Matrix,AMREX_SPACEDIM> gradEpsp;
	CrystalPlastic(Solid<Set::Sym::MajorMinor> base) : Solid<Set::Sym::MajorMinor>(base) {};
	CrystalPlastic() 
	{
		initializeSlip();
	} 
	CrystalPlastic(Set::Matrix Fgb)
	{
		esp = Fgb;
	}
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	CrystalPlastic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	virtual ~CrystalPlastic() {};

	//virtual Set::Matrix operator () (Set::Matrix &eps,bool =true) const;
	//virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps,bool =true);
	Set::Scalar W(Set::Matrix gradu)
	{
		return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
	}
	Set::Matrix DW(Set::Matrix gradu)
	{
		return ddw*gradu;
	}
	Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(Set::Matrix /*gradu*/)
	{
		return ddw;
	}
//----------------------//

	std::array<double,12> StressSlipSystem(const Set::Matrix& sig);
	void Setdt(const double _dt);
	Set::Matrix GetEsp() const;
	double getGamma(int index) const;
	void update(const Set::Matrix es, Set::Matrix& sigma);

	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	
	static CrystalPlastic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static CrystalPlastic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        CrystalPlastic ret;
	    Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        ret.define(C11,C12,C44,phi1,Phi,phi2);
        return ret;
    }
/*
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
	*/
};

}
}
}

#endif 