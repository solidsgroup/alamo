#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Finite/Adhesion.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "Model/Solid/Finite/PseudoLinear/Cubic.H"
#include "Model/Solid/Finite/PseudoAffine/Cubic.H"
#include "Model/Solid/Linear/Laplacian.H"
#include "Model/Solid/Affine/J2.H"
#include "Model/Solid/Affine/Hexagonal.H"

#include "Integrator/Mechanics.H"
#include "Model/Solid/Finite/CrystalPlastic.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    // which integrator to use (can only be mechanics)
    pp.query_validate("alamo.program",program,{"mechanics"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    if (program == "mechanics")
    {
        std::string model = "linear.isotropic";
        // which mechanics model to use
        pp.query_default("alamo.program.mechanics.model",model,"linear.isotropic");
        if (model == "linear.isotropic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Isotropic>>(integrator);
        else if (model == "linear.cubic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Cubic>>(integrator);
        else if (model == "affine.cubic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Cubic>>(integrator);
        else if (model == "affine.hexagonal")
            pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Hexagonal>>(integrator);
        else if (model == "affine.isotropic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Isotropic>>(integrator);
        else if (model == "linear.laplacian")
            pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Laplacian>>(integrator);
        else if (model == "finite.neohookean")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::NeoHookean>>(integrator);
        else if (model == "finite.neohookeanpre")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::NeoHookeanPredeformed>>(integrator);
        else if (model == "finite.adhesion")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::Adhesion>>(integrator);
        else if (model == "finite.pseudolinear.cubic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::PseudoLinear::Cubic>>(integrator);
        else if (model == "finite.pseudoaffine.cubic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::PseudoAffine::Cubic>>(integrator);
        else if (model == "affine.j2")
            pp.select_only<Integrator::Mechanics<Model::Solid::Affine::J2>>(integrator);
        else if (model == "finite.crystalplastic")
            pp.select_only<Integrator::Mechanics<Model::Solid::Finite::CrystalPlastic>>(integrator);
        else Util::Abort(INFO,model," is not a valid model");
    }
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
