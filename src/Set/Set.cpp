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


}

namespace Util
{
Set::Scalar Random()
{
	return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}
}