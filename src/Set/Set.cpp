#include "Set.H"
namespace Set
{
//	AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
//	Set::Matrix operator * (Matrix4<3,Sym::Isotropic> a, Set::Matrix b)
//	{
//		Set::Matrix ret;
//		
//		ret(0,0) = (a.lambda + 2.*a.mu) * b(0,0) +       a.lambda      *b(1,1) +       a.lambda      *b(2,2);
//		ret(1,1) =        a.lambda      * b(0,0) + (a.lambda + 2.*a.mu)*b(1,1) +       a.lambda      *b(2,2);
//		ret(2,2) =        a.lambda      * b(0,0) +       a.lambda      *b(1,1) + (a.lambda + 2.*a.mu)*b(2,2);
//
//		ret(1,2) = a.mu*(b(1,2) + b(2,1)); ret(2,1) = ret(1,2);
//		ret(2,0) = a.mu*(b(2,0) + b(0,2)); ret(0,2) = ret(2,0);
//		ret(0,1) = a.mu*(b(0,1) + b(1,0)); ret(1,0) = ret(0,1);
//		
//		return ret;
//	}

	AMREX_GPU_HOST_DEVICE
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
	AMREX_GPU_HOST_DEVICE
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
	AMREX_GPU_HOST_DEVICE 
	Matrix4<3,Sym::Major> operator - (Matrix4<3,Sym::Major> a, Matrix4<3,Sym::Major> b)
	{
		Matrix4<3,Sym::Major> ret;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						ret(i,j,k,l) = a(i,j,k,l) - b(i,j,k,l);
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
