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
#include "Model/Solid/Finite/NeoHookean.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "Model/Solid/Linear/Transverse.H"
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

    pp.query_switch("alamo.program.mechanics.model",{
            {"linear.isotropic",          [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Isotropic>>(integrator);             }},
            { "linear.cubic",             [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Cubic>>(integrator);                 }},
            { "affine.cubic",             [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Cubic>>(integrator);                 }},
            { "affine.hexagonal",         [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Hexagonal>>(integrator);             }},
            { "affine.isotropic",         [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Affine::Isotropic>>(integrator);             }},
            { "linear.laplacian",         [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Laplacian>>(integrator);             }},
            { "finite.neohookean",        [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Finite::NeoHookean>>(integrator);            }},
            { "finite.neohookeanpre",     [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Finite::NeoHookeanPredeformed>>(integrator); }},
            { "linear.transverse",        [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Linear::Transverse> >(integrator);           }},
            { "finite.pseudolinear.cubic",[&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Finite::PseudoLinear::Cubic>>(integrator);   }},
            { "finite.pseudoaffine.cubic",[&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Finite::PseudoAffine::Cubic>>(integrator);   }},
            { "affine.j2",                [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Affine::J2>>(integrator);                    }},
            { "finite.crystalplastic",    [&]() {pp.select_only<Integrator::Mechanics<Model::Solid::Finite::CrystalPlastic>>(integrator);        }}
        });

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
