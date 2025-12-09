#include "Set.H"
#include <random>
namespace Set
{
}

namespace Set
{
namespace Constant
{

    // Pi
    const Set::Scalar Pi  = 3.14159265359;
    Set::Scalar Rg        = 0.0;
    Set::Scalar Boltzmann = 0.0;
    Set::Scalar Avogadro  = 0.0;

    // Function to update constants to match system level units
    void SetGlobalConstants()
    {
        // Gas constant
        Rg        = Unit::Parse("8.31446261815_J/mol/K").normalized_value();
        // Boltzmann constant
        Boltzmann = Unit::Parse("1.380649e-23_J/K").normalized_value();
        // Avogadro's number
        Avogadro  = Unit::Parse("6.02214076e23_1/mol").normalized_value();
    }

}
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
