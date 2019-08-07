#include "Set.H"
namespace Set
{
	Eigen::Matrix<amrex::Real,3,3> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,3,3> b)
	{
    	Eigen::Matrix<amrex::Real,3,3> ret = Eigen::Matrix<amrex::Real,3,3>::Zero();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						ret(i,j) += a(i,j,k,l)*b(k,l);
		return ret;
	}
	Eigen::Matrix<amrex::Real,2,2> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,2,2> b)
	{
    	Eigen::Matrix<amrex::Real,2,2> ret = Eigen::Matrix<amrex::Real,2,2>::Zero();
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					for (int l = 0; l < 2; l++)
						ret(i,j) += a(i,j,k,l)*b(k,l);
		return ret;
	}
}

namespace Util
{
Set::Scalar Random()
{
	return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}
}