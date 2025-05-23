#include "Set.H"
#include "AMReX_ParallelDescriptor.H"
#include <random>
namespace Set
{
}

namespace Util
{
Set::Scalar Random(bool sync)
{
    if (sync)
    {
        Set::Scalar ret = NAN;

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            ret = amrex::Random();
        }
        amrex::ParallelDescriptor::Bcast(&ret, 1);
        return ret;
    }

    return amrex::Random();
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
