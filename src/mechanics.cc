#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Isotropic.H"
//#include "Model/Solid/Elastic/NeoHookean.H"
//#include "Model/Solid/Elastic/PseudoLinearCubic.H"
//#include "Model/Solid/Elastic/PseudoLinearCubicPredeformed.H"
//#include "Model/Solid/Linear/Laplacian.H"
#include "Model/Solid/Affine/J2.H"
//#include "Model/Solid/Affine/J2Plastic.H"
//#include "Model/Solid/Affine/Hexagonal.H"

#include "Integrator/Mechanics.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp.query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    std::string model = "linear.isotropic";
    pp.query("alamo.program.mechanics.model",model);
    //if (model == "linear.isotropic")        integrator = new Integrator::Mechanics<Model::Solid::Linear::Isotropic>(pp);
    //else if (model == "linear.cubic")       integrator = new Integrator::Mechanics<Model::Solid::Linear::Cubic>(pp);
    //else if (model == "affine.cubic")       integrator = new Integrator::Mechanics<Model::Solid::Affine::Cubic>(pp);
    //else if (model == "affine.hexagonal")   integrator = new Integrator::Mechanics<Model::Solid::Affine::Hexagonal>(pp);
    if (model == "affine.isotropic")   
        integrator = new Integrator::Mechanics<Model::Solid::Affine::Isotropic>(pp);
    //else if (model == "linear.laplacian")   integrator = new Integrator::Mechanics<Model::Solid::Linear::Laplacian>(pp);
    //else if (model == "elastic.neohookean") integrator = new Integrator::Mechanics<Model::Solid::Elastic::NeoHookean>(pp);
    //else if (model == "elastic.neohookeanpre") integrator = new Integrator::Mechanics<Model::Solid::Elastic::NeoHookeanPredeformed>(pp);
    //else if (model == "elastic.pseudolinearcubic") integrator = new Integrator::Mechanics<Model::Solid::Elastic::PseudoLinearCubic>(pp);
    //else if (model == "elastic.pseudolinearcubicpredeformed") integrator = new Integrator::Mechanics<Model::Solid::Elastic::PseudoLinearCubicPredeformed>(pp);
    else if (model == "affine.j2")          
        integrator = new Integrator::Mechanics<Model::Solid::Affine::J2>(pp);
    //else if (model == "affine.j2plastic")          
    //    integrator = new Integrator::Mechanics<Model::Solid::Affine::J2Plastic>(pp);
    //else Util::Abort(INFO,model," is not a valid model");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
