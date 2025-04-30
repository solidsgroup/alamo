#include "IO/ParmParse.H"
#include "Integrator/HeatConduction.H"
int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);
    srand(2);
    Integrator::Integrator *integrator = nullptr;
    IO::ParmParse pp;
    pp.select_only<Integrator::HeatConduction>(integrator);
    integrator->InitData();
    integrator->Evolve();
    delete integrator;
    Util::Finalize();
} 
