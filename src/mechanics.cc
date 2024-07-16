#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

//#include "Model/Solid/Linear/Isotropic.H"
//#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Linear/Cubic.H"
//#include "Model/Solid/Affine/Cubic.H"
//#include "Model/Solid/Finite/NeoHookean.H"
//#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
//#include "Model/Solid/Finite/PseudoLinearCubic.H"
//#include "Model/Solid/Finite/PseudoLinearCubicPredeformed.H"
//#include "Model/Solid/Linear/Laplacian.H"
//#include "Model/Solid/Affine/J2.H"
//#include "Model/Solid/Affine/Hexagonal.H"

#include "Integrator/Mechanics.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "mechanics";
    IO::ParmParse pp;
    pp_query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if (program == "mechanics")
    {
        std::string model = "linear.isotropic";
        pp_query("alamo.program.mechanics.model",model);
        // if (model == "linear.isotropic")        integrator = new Integrator::Mechanics<Model::Solid::Linear::Isotropic>(pp);
        // else
        if (model == "linear.cubic")       integrator = new Integrator::Mechanics<Model::Solid::Linear::Cubic>(pp);
        // else if (model == "affine.cubic")       integrator = new Integrator::Mechanics<Model::Solid::Affine::Cubic>(pp);
        // else if (model == "affine.hexagonal")   integrator = new Integrator::Mechanics<Model::Solid::Affine::Hexagonal>(pp);
        // else if (model == "affine.isotropic")   integrator = new Integrator::Mechanics<Model::Solid::Affine::Isotropic>(pp);
        // else if (model == "linear.laplacian")   integrator = new Integrator::Mechanics<Model::Solid::Linear::Laplacian>(pp);
        // else if (model == "finite.neohookean") integrator = new Integrator::Mechanics<Model::Solid::Finite::NeoHookean>(pp);
        // else if (model == "finite.neohookeanpre") integrator = new Integrator::Mechanics<Model::Solid::Finite::NeoHookeanPredeformed>(pp);
        // else if (model == "finite.pseudolinearcubic") integrator = new Integrator::Mechanics<Model::Solid::Finite::PseudoLinearCubic>(pp);
        // else if (model == "finite.pseudolinearcubicpredeformed") integrator = new Integrator::Mechanics<Model::Solid::Finite::PseudoLinearCubicPredeformed>(pp);
        // else if (model == "affine.j2")          integrator = new Integrator::Mechanics<Model::Solid::Affine::J2>(pp);
        else Util::Abort(INFO,model," is not a valid model");
    }
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
