#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Finite/CrystalPlastic.H"

#include "Integrator/LagrangianShockCP.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    pp.query_validate("alamo.program", program, {"lagrangian_shock_cp"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    if (program == "lagrangian_shock_cp")
    {
        std::string model = "finite.crystalplastic";
        pp.query_default("alamo.program.mechanics.model", model, "finite.crystalplastic");

        if (model == "finite.crystalplastic")
            pp.select_only<Integrator::LagrangianShockCP<Model::Solid::Finite::CrystalPlastic>>(integrator);
        else if (model == "linear.isotropic")
            pp.select_only<Integrator::LagrangianShockCP<Model::Solid::Linear::Isotropic>>(integrator);
        else
            Util::Abort(INFO, model, " is not a valid model for lagrangian_shock_cp");
    }
    else
        Util::Abort(INFO, "Error: \"", program, "\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}