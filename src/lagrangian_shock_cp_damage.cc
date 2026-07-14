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

#include "Integrator/LagrangianShockCPDamage.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    pp.query_validate("alamo.program", program, {"lagrangian_shock_cp_damage"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    if (program == "lagrangian_shock_cp_damage")
    {
        std::string model = "finite.crystalplastic";
        pp.query_default("alamo.program.mechanics.model", model, "finite.crystalplastic");

        if (model == "finite.crystalplastic")
            pp.select_only<Integrator::LagrangianShockCPDamage<Model::Solid::Finite::CrystalPlastic>>(integrator);
        else if (model == "linear.isotropic")
            pp.select_only<Integrator::LagrangianShockCPDamage<Model::Solid::Linear::Isotropic>>(integrator);
        else
            Util::Abort(INFO, model, " is not a valid model for lagrangian_shock_cp_damage");
    }
    else
        Util::Abort(INFO, "Error: \"", program, "\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}