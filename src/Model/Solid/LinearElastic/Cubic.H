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
///
/// \class Cubic
///
/// Cubic implements the following linear elastic model, expressed in Voigt notation
///
/// \f[\begin{bmatrix}\hat{\sigma}_{11}\\\hat{\sigma}_{22}\\\hat{\sigma}_{33}\\\hat{\sigma}_{23}\\\hat{\sigma}_{31}\\\hat{\sigma}_{12}\end{bmatrix} = 
///    \begin{bmatrix} C_{11} & C_{12} & C_{12} \\ %
///                    C_{12} & C_{11} & C_{12} \\ %
///                    C_{12} & C_{12} & C_{11} \\ %
///                    &&& C_{44}   \\ %
///                    &&&& C_{44}  \\ %
///                    &&&&& C_{44} \\ %
///    \end{bmatrix}
///    \begin{bmatrix}\hat{\varepsilon}_{11}\\\hat{\varepsilon}_{22}\\\hat{\varepsilon}_{33}\\\hat{\varepsilon}_{23}\\\hat{\varepsilon}_{31}\\\hat{\varepsilon}_{12}\end{bmatrix} 
/// \f]
/// where \f$\{\hat{\sigma}\}\f$ are the components of \f$\hat{\mathbf{\sigma}}\f$ and
/// \f$\{\hat{\varepsilon}\}\f$ are the components of \f$\hat{\mathbf{\varepsilon}}\f$.
/// \f$\hat{\mathbf{\sigma}}\f$, \f$\hat{\mathbf{\varepsilon}}\f$ are rotations of the input \f$\mathbf{\varepsilon} = \frac{1}{2}(\nabla\mathbf{u} + \nabla\mathbf{u}^T)\f$
/// and the output \f$\mathbf{\sigma}\f$; i.e.
/// \f[ \hat{\mathbf{\sigma}} = R \mathbf{\sigma} \f]
/// \f[ \hat{\mathbf{\varepsilon}} = R \mathbf{\varepsilon} \f]
/// where \f$R\in SO(3)\f$.
/// Note that \f$R\f$ is _always_ three dimensional, even if the code is compiled in 2D.
/// Consequently, the eigen Matrix3d type is used instead of the usual Set::Matrix.
///
/// Rotations can be specified in one of two ways:
///    1. Passing a rotation matrix explicitly
///    2. Passing in Bunge Euler angles
/// If no matrix is specified, the Identity is used by default.
/// 
class Cubic : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::MajorMinor>
{
public:
	Cubic() {} ;
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	~Cubic() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps,bool =true) const;
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps,bool =true);

	virtual Set::Scalar    W (Set::Matrix &eps) const;
	virtual Set::Matrix   DW (Set::Matrix &eps) const;
	virtual Set::Vector   DW (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps);
	virtual Set::Matrix4<3,Set::Sym::MajorMinor> DDW (Set::Matrix &eps) const;

	void Randomize();
	void Randomize(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44);
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity());
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2);
	
	Cubic operator + (const Cubic &rhs) const
	{Cubic ret = *this; ret.C += rhs.C; return ret;}

	Cubic operator - (const Cubic &rhs) const
	{Cubic ret = *this; ret.C -= rhs.C; return ret;}

	Cubic operator * (const Cubic &rhs) const
	{Cubic ret = *this; ret.C *= rhs.C; return ret;}

	Cubic operator / (const Cubic &rhs) const
	{Cubic ret = *this; ret.C /= rhs.C; return ret;}

	Cubic operator * (const Set::Scalar alpha) const
	{Cubic ret = *this; ret.C *= alpha; return ret;}

	Cubic operator / (const Set::Scalar alpha) const
	{Cubic ret = *this; ret.C /= alpha; return ret;}

	void operator += (const Cubic &rhs) {C += rhs.C;}

private:
	Set::Matrix4<3,Set::Sym::MajorMinor> C;
	//std::array<Set::Scalar,21> C;

};

}
}
}

#endif 


