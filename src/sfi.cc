#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/SFI.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    std::string program;
    pp.query("alamo.program",program);
    srand(2);

    Integrator::SFI integrator(pp);
    integrator.InitData();
    integrator.Evolve();
    
    Util::Finalize();
} 
