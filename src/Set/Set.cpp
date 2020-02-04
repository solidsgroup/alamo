#include "Set.H"
#include <random>
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
Set::Scalar Gaussian(amrex::Real mean,amrex::Real std_deviation)
{
	std::random_device randomness_device{};
    std::mt19937 pseudorandom_generator{randomness_device()};
    std::normal_distribution<double> distribution{mean, std_deviation};
	auto sample = distribution(pseudorandom_generator);
	return sample;
}

}
