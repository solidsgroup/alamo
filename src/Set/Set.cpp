#include "Set.H"
#include "AMReX_RandomEngine.H"
#include <random>
namespace Set
{
}

namespace Util
{
Set::Scalar Random()
{
    return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}
// AMREX_GPU_HOST_DEVICE
// Set::Scalar Random(int i, int j, int k, int n, int lev)
// {
//     amrex::Long key1    = (static_cast<amrex::Long>(i) << 32) | static_cast<unsigned int>(j);
//     amrex::Long key2    = (static_cast<amrex::Long>(k) << 32) | static_cast<unsigned int>(n);
//     amrex::Long counter =  static_cast<amrex::Long>(lev);

//     //auto engine = amrex::RandomEngine(key1, key2, static_cast<amrex::Long>(lev));

//     return amrex::Random();// (key1, key2, counter);
// }
Set::Scalar Gaussian(amrex::Real mean,amrex::Real std_deviation)
{
    std::random_device randomness_device{};
    std::mt19937 pseudorandom_generator{randomness_device()};
    std::normal_distribution<double> distribution{mean, std_deviation};
    auto sample = distribution(pseudorandom_generator);
    return sample;
}

}
