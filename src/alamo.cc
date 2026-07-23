// This is the main entry point for alamo and is a general-purpose launcher for
// many of the main integrators.
// Check the possible values for :code:`alamo.program` below to see the possible
// integrators that can be launched.
//

#include <iostream>
#include <fstream>
#include <iomanip>

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
    // This input determines which integrator is used.
    pp.query_validate(  "alamo.program", program,
                        {"microstructure", "flame", "heat", "dendrite","allencahn","cahnhilliard","pfc"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    pp.query_switch("alamo.program",
                    pp.forward_args(
                        "microstructure",  [&] ()
                        {
                            std::string model;
                            // This input determines which elastic model is used - only if using
                            // the PhaseFieldMicrostructure integrator.
                            pp.query_validate(  "alamo.program.microstructure.model",model,
                                                {"affine.cubic","affine.hexagonal","finite.pseudoaffine.cubic"});

                            pp.query_switch("alamo.program.microstructure.model",
                                            pp.forward_args(
                                                "affine.cubic", [&]()
                                                {
                                                    pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>>(integrator);
                                                },
                                                "affine.hexagonal", [&]()
                                                {
                                                    pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Hexagonal>>(integrator);
                                                },
                                                "finite.pseudoaffine.cubic", [&]()
                                                {
                                                    pp.select_only<Integrator::PhaseFieldMicrostructure<Model::Solid::Finite::PseudoAffine::Cubic>>(integrator);
                                                }));
                        },
                        "flame", [&] ()        {pp.select_only<Integrator::Flame>(integrator);} ,
                        "heat", [&] ()         {pp.select_only<Integrator::HeatConduction>(integrator);},
                        "fracture", [&] ()     {pp.select_only<Integrator::Fracture>(integrator);},
                        "dendrite", [&] ()     {pp.select_only<Integrator::Dendrite>(integrator);},
                        "allencahn", [&] ()    {pp.select_only<Integrator::AllenCahn>(integrator);},
                        "cahnhilliard", [&] () {pp.select_only<Integrator::CahnHilliard>(integrator);},
                        "pfc", [&] ()          {pp.select_only<Integrator::PFC>(integrator); }
                    ));


    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
