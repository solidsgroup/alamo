#include "Integrator/NarrowBandLevelset.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    srand(2);

    {
    Integrator::NarrowBandLevelset integrator(pp);
    }
    
    Util::Finalize();
} 
