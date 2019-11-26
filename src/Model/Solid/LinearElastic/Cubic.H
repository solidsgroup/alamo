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
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity())
	{
		define(C11, C12, C44, R);
	}
	Cubic(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
	{
		define(C11, C12, C44, phi1, Phi, phi2);
	}
	~Cubic() {} ;

	AMREX_FORCE_INLINE
	virtual Set::Matrix operator () (Set::Matrix &a_gradu,bool a_homogeneous=true) const
	{
		Set::Matrix gradu = a_gradu;
		if (!a_homogeneous) gradu -= F0;
		return C*gradu;
	}

	AMREX_FORCE_INLINE
	virtual Set::Vector operator () (Set::Matrix3 &a_gradgradu,bool a_homogeneous=true)
	{
		Set::Matrix3 gradgradu = a_gradgradu;
		if (!a_homogeneous) gradgradu -= gradF0;
		Set::Vector ret = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j = 0; j < AMREX_SPACEDIM; j++)
				for (int k = 0; k < AMREX_SPACEDIM; k++)
					for (int l = 0; l < AMREX_SPACEDIM; l++)
						ret(i) += C(i,j,k,l)*gradgradu(k,l,j);
		return ret;
	}

	virtual Set::Scalar    W (Set::Matrix &gradu) const
	{
		Util::Abort(INFO,"Depricated");
		Set::Matrix sig = C*gradu;
		return 0.5 * (sig*gradu).trace();
	}

	virtual Set::Matrix   DW (Set::Matrix &gradu) const
	{
		Util::Abort(INFO,"Depricated");
		return (*this)(gradu);
	}

	virtual Set::Vector   DW (Set::Matrix3 &gradgradu)
	{
		Util::Abort(INFO,"Depricated");
		return (*this)(gradgradu);
	}

	virtual Set::Matrix4<3,Set::Sym::MajorMinor> DDW (Set::Matrix & /*gradu*/) const
	{
		Util::Abort(INFO,"Depricated");
		return C;
	}

	void Randomize()
	{
		Set::Scalar C11 = 0.5 + 0.5*Util::Random();
		Set::Scalar C12 = 0.5 + 0.5*Util::Random();
		Set::Scalar C44 = 0.5 + 0.5*Util::Random();
		Randomize(C11,C12,C44);
	}
	void Randomize(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
	{
		Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
		Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
		Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
		define(C11,C12,C44,phi1,Phi,phi2);
	}

	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R = Eigen::Matrix3d::Identity())
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
	void define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
	{
		Eigen::Matrix3d m;
		m =     Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
			Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
		 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
		define(C11,C12,C44,m);
	}
	
	AMREX_FORCE_INLINE
	Cubic operator + (const Cubic &rhs) const
	{Cubic ret = *this; ret.C += rhs.C; ret.F0 += rhs.F0; ret.gradF0 += rhs.gradF0; return ret;}

	AMREX_FORCE_INLINE
	Cubic operator - (const Cubic &rhs) const
	{Cubic ret = *this; ret.C -= rhs.C; ret.F0 -= rhs.F0; ret.gradF0 -= rhs.gradF0; return ret;}

	AMREX_FORCE_INLINE
	Cubic operator * (const Set::Scalar alpha) const
	{Cubic ret = *this; ret.C *= alpha; ret.F0 *= alpha; ret.gradF0 *= alpha; return ret;}

	AMREX_FORCE_INLINE
	Cubic operator / (const Set::Scalar alpha) const
	{Cubic ret = *this; ret.C /= alpha; ret.F0 /= alpha; ret.gradF0 /= alpha; return ret;}

	AMREX_FORCE_INLINE
	void operator += (const Cubic &rhs) 
	{C += rhs.C; F0 += rhs.F0; gradF0 += rhs.gradF0;}

private:
	Set::Matrix4<3,Set::Sym::MajorMinor> C;

public:
	Set::Matrix F0      = Set::Matrix::Zero();
	Set::Matrix3 gradF0 = Set::Matrix3::Zero();

};

}
}
}

#endif 


