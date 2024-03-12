#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Mechanics.H"
#include "Model/Solid/Elastic/CrystalPlastic.H"
#include "Model/Solid/Elastic/PseudoLinearCubic.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp.query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if (program == "mechanics")
    {
        std::string model = "linear.isotropic";
        pp.query("alamo.program.mechanics.model",model);
        if (model == "elastic.pseudolinearcubic")        integrator = new Integrator::Mechanics<Model::Solid::Elastic::PseudoLinearCubic>(pp);
        else if (model == "elastic.crystalplastic")        integrator = new Integrator::Mechanics<Model::Solid::Elastic::CrystalPlastic>(pp);
        else Util::Abort(INFO,model," is not a valid model");
    }
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
