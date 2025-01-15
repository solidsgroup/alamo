#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/ScimitarX.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    srand(2);

    Integrator::ScimitarX integrator(pp);
    integrator.InitData();
    integrator.Evolve();
    
    Util::Finalize();
} 
