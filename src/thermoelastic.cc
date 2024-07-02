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

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp_query("alamo.program",program);
    srand(2);

    Integrator::ThermoElastic integrator(pp);
    integrator.InitData();
    integrator.Evolve();
    
    Util::Finalize();
} 
