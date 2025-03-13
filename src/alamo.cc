#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX_ParallelDescriptor.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Affine/Hexagonal.H"
#include "Model/Solid/Finite/PseudoAffine/Cubic.H"

#include "Integrator/AllenCahn.H"
#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mechanics.H"
#include "Integrator/Flame.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"
#include "Integrator/ThermoElastic.H"
#include "Integrator/Dendrite.H"
#include "Integrator/PFC.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    // The integrator to select
    pp.query_default("alamo.program",program,"microstructure");
    srand(2 + amrex::ParallelDescriptor::MyProc());

    Integrator::Integrator *integrator = nullptr;
    if (program == "microstructure")
    {
        std::string model;
        // which model to use (if using PFM with elasticity)
        pp.query_default("alamo.program.microstructure.model",model,"affine.cubic");
        if (model == "affine.cubic")
            pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>>(integrator);
        else if (model == "affine.hexagonal")
            pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Hexagonal>>(integrator);
        else if (model == "finite.pseudoaffine.cubic")
            pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Finite::PseudoAffine::Cubic>>(integrator);
        else Util::Abort(INFO,model," is not a valid model");
    }
    else if (program == "flame")                pp.select_only<Integrator::Flame>(integrator);
    else if (program == "heat")                 pp.select_only<Integrator::HeatConduction>(integrator);
    else if (program == "thermoelastic")        pp.select_only<Integrator::ThermoElastic>(integrator);
    //else if (program == "fracture")             pp.select_only<Integrator::Fracture>();
    else if (program == "dendrite")             pp.select_only<Integrator::Dendrite>(integrator);
    else if (program == "allencahn")            pp.select_only<Integrator::AllenCahn>(integrator);
    else if (program == "cahnhilliard")         pp.select_only<Integrator::CahnHilliard>(integrator);
    else if (program == "pfc")                  pp.select_only<Integrator::PFC>(integrator);
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
