#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/ThermoElastic.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    srand(2);

    Integrator::Integrator *integrator;
    pp.select_only<Integrator::ThermoElastic>(integrator);
    integrator->InitData();
    integrator->Evolve();
    delete integrator;
    
    Util::Finalize();
} 
