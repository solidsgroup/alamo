#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/LagrangianShock.H"

// For now, use a simple isotropic linear elastic model as placeholder.
// The MODEL is not yet coupled to stress computation (EOS only),
// but the template is required by Base::Mechanics.
#include "Model/Solid/Linear/Isotropic.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    pp.query_validate("alamo.program",program,{"lagrangian_shock"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    if (program == "lagrangian_shock")
    {
        pp.select_only<Integrator::LagrangianShock<Model::Solid::Linear::Isotropic>>(integrator);
    }
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}