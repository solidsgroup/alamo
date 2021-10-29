#include "Set.H"
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
Set::Scalar Gaussian(amrex::Real mean,amrex::Real std_deviation)
{
    std::random_device randomness_device{};
    std::mt19937 pseudorandom_generator{randomness_device()};
    std::normal_distribution<double> distribution{mean, std_deviation};
    auto sample = distribution(pseudorandom_generator);
    return sample;
}

}
