#ifndef MODEL_SOLID_LINEARELASTIC_CUBIC_H_
#define MODEL_SOLID_LINEARELASTIC_CUBIC_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"


namespace Model
{
namespace Solid
{
namespace LinearElastic
{

class Cubic : public Model::Solid::LinearElastic::LinearElastic
{
public:
	Cubic() {} ;
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	~Cubic() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps) const;
#if AMREX_SPACEDIM==2
	virtual Set::Vector operator () (std::array<Set::Matrix,2> &gradeps);
#elif AMREX_SPACEDIM==3
	virtual Set::Vector operator () (std::array<Set::Matrix,3> &gradeps);
#endif

	void Randomize();
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	
	Cubic operator + (const Cubic &rhs) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] + rhs.C[i]; return ret;}

	void operator += (const Cubic &rhs)
	{for (unsigned int i = 0; i < C.size(); i++) C[i] += rhs.C[i];}

	Cubic operator - (const Cubic &rhs) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] - rhs.C[i]; return ret;}

	Cubic operator * (const Cubic &rhs) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] * rhs.C[i]; return ret;}

	Cubic operator / (const Cubic &rhs) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i] / rhs.C[i]; return ret;}

	Cubic operator / (const Set::Scalar alpha) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i]/alpha; return ret;}

	Cubic operator * (const Set::Scalar alpha) const
	{Cubic ret; for (unsigned int i = 0; i < C.size(); i++) ret.C[i] = C[i]*alpha; return ret;}
private:
	std::array<Set::Scalar,21> C;

};

}
}
}

#endif 


