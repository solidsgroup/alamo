///
/// \file PhaseFieldMicrostructure.H
///
#ifndef INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H
#define INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"
#include "Integrator/Mechanics.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Step.H"
#include "IC/Constant.H"
#include "IC/TabulatedInterface.H"
#include "IC/PerturbedInterface.H"
#include "IC/Voronoi.H"
#include "IC/Sphere.H"
#include "IC/Expression.H"

#include "Model/Interface/GB/GB.H"
#include "Model/Interface/GB/Sin.H"
#include "Model/Interface/GB/AbsSin.H"
#include "Model/Interface/GB/Read.H"
#include "Model/Interface/GB/SH.H"

#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Operator/Elastic.H"

namespace Integrator
{
using model_type = Model::Solid::Affine::Cubic;

enum RegularizationType{
    Wilmore = 1,
    K12 = 2};

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public Base::Mechanics<model_type>
{
public:
    PhaseFieldMicrostructure() : Integrator() {};
    PhaseFieldMicrostructure(IO::ParmParse &pp) : Integrator() 
    {pp.queryclass(*this);}

    static void Parse(PhaseFieldMicrostructure &value, IO::ParmParse &pp)
    {
        BL_PROFILE("PhaseFieldMicrostructure::Parse");
        
        pp.query("pf.number_of_grains", value.number_of_grains);     // Number of grain fields (may be more if using different IC)
        pp.query("pf.M", value.pf.M);                                // Mobility
        if (pp.contains("pf.mu")) pp.query("pf.mu", value.pf.mu);       // Phase field :math:`\mu`
        pp.query("pf.gamma", value.pf.gamma);                        // Phase field :math:`\gamma`
        pp.query("pf.sigma0", value.pf.sigma0);                      // Initial GB energy if not using GB anisotropy
        pp.query("pf.l_gb", value.pf.l_gb);                          // Mobility
        pp.query("pf.elastic_mult",value.pf.elastic_mult);           // Multiplier of elastic energy
        pp.query("pf.elastic_threshold",value.pf.elastic_threshold); // Elastic threshold (:math:`\phi_0`)
        value.pf.L = (4./3.)*value.pf.M / value.pf.l_gb;
    
        pp.query("amr.max_level", value.max_level);            // Maximum AMR level
        pp.query("amr.ref_threshold", value.ref_threshold);    // Phase field refinement threshold

        // Elasticity
        pp.query("mechanics.tstart",value.mechanics.tstart);
        //value.mechanics.model.resize(value.number_of_grains);

        pp.queryclass("mechanics.model",value.mechanics.model);
        value.m_type = Base::Mechanics<model_type>::Disable; // Turn mechanics off by default
        pp.queryclass<Base::Mechanics<model_type>>("mechanics",&value);
        if (value.m_type == Base::Mechanics<model_type>::Type::Static) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,3);
        


        // Lagrange multiplier method for enforcing volumes
        pp.query("lagrange.on", value.lagrange.on);               // Turn on
        if (value.lagrange.on)
        {
            pp.query("lagrange.lambda", value.lagrange.lambda);   // Lagrange multiplier value
            pp.query("lagrange.vol0", value.lagrange.vol0);       // Prescribed volume
            pp.query("lagrange.tstart", value.lagrange.tstart);   // Time to start enforcing Lagrange multipler
            value.SetThermoInt(1);
        }

        pp.query("sdf.on",value.sdf.on);
        if(value.sdf.on)
        {
            std::vector<std::string> vals;
            pp.queryarr("val",vals);
            if (vals.size() == 1)
                for (int i = 0; i < value.number_of_grains; i++)
                    value.sdf.val.push_back(Numeric::Interpolator::Linear<Set::Scalar>(vals[0]));
            else if (vals.size() == value.number_of_grains)
                for (int i = 0; i < value.number_of_grains; i++)
                    value.sdf.val.push_back(Numeric::Interpolator::Linear<Set::Scalar>(vals[i]));
            else
                Util::Abort(INFO,"sdf.val received ", vals.size(), " but requires 1 or ", value.number_of_grains);

            pp.query("tstart",value.sdf.tstart);
        }

        // Anisotropic grain boundary energy parameters
        pp.query("anisotropy.on", value.anisotropy.on);                            // Turn on
        pp.query("anisotropy.beta", value.anisotropy.beta);                        // Regularization param 
        pp.query("anisotropy.tstart", value.anisotropy.tstart);                    // Time to turn on anisotropy
        value.anisotropy.timestep = value.timestep;
        pp.query("anisotropy.timestep", value.anisotropy.timestep);                // Modify timestep when turned on
        value.anisotropy.plot_int = value.plot_int;
        pp.query("anisotropy.plot_int", value.anisotropy.plot_int);                // Modify plot_int when turned on
        value.anisotropy.plot_dt = value.plot_dt;
        pp.query("anisotropy.plot_dt", value.anisotropy.plot_dt);                  // Modify plot_dt when turned on
        pp.query("anisotropy.thermo_int", value.anisotropy.thermo_int);            // Modify thermo int when turned on
        pp.query("anisotropy.thermo_plot_int", value.anisotropy.thermo_plot_int);  // Modify thermo plot int when turned on
        pp.query("anisotropy.elastic_int",value.anisotropy.elastic_int);           // Frequency of elastic calculation
        if (value.anisotropy.on) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,2);

        // Determine the kind of regularization to use
        std::map<std::string, RegularizationType> regularization_type;
        regularization_type["wilmore"] = RegularizationType::Wilmore;
        regularization_type["k12"] = RegularizationType::K12;
        std::string regularization_type_input = "k12";
        pp.query("anisotropy.regularization", regularization_type_input);    // Type of regularization to use  
        value.regularization = regularization_type[regularization_type_input];

        // Set the anisotropic GB model
        std::string gb_type = "";
        pp.query("anisotropy.gb_type", gb_type); // Type of GB energy to use
        if      (gb_type == "abssin") value.boundary = new Model::Interface::GB::AbsSin(pp,"anisotropy.abssin");
        else if (gb_type == "sin")    value.boundary = new Model::Interface::GB::Sin(pp,"anisotropy.sin");
        else if (gb_type == "read")   value.boundary = new Model::Interface::GB::Read(pp,"anisotropy.read");
        else if (gb_type == "sh")     value.boundary = new Model::Interface::GB::SH(pp,"anisotropy.sh");
        else if (value.anisotropy.on)       Util::Abort(INFO,"A GB model must be specified");

        // Thermal disconnection fluctuations
        pp.query("fluctuation.on",value.fluctuation.on);
        pp.query("fluctuation.amp",value.fluctuation.amp);
        pp.query("fluctuation.sd",value.fluctuation.sd);
        pp.query("fluctuation.tstart", value.fluctuation.tstart);
        value.fluctuation.norm_dist = std::normal_distribution<double>(0.0,value.fluctuation.sd);

        // Disconnection generation
        pp.query("disconnection.on",value.disconnection.on);
        pp.query("disconnection.tstart", value.disconnection.tstart);
        pp.query("disconnection.nucleation_energy",value.disconnection.nucleation_energy);
        pp.query("disconnection.tau_vol",value.disconnection.tau_vol);
        pp.query("disconnection.temp",value.disconnection.temp);
        pp.query("disconnection.box_size",value.disconnection.box_size);
        pp.query("disconnection.interval",value.disconnection.interval);
        pp.query("disconnectionepsilon",value.disconnection.epsilon);
        pp.query("disconnection.fixed.on",value.disconnection.fixed.on);
        if (value.disconnection.fixed.on)
        {
            pp.queryarr("disconnection.fixed.sitex",value.disconnection.fixed.sitex);
            pp.queryarr("disconnection.fixed.sitey",value.disconnection.fixed.sitey);
            pp.queryarr("disconnection.fixed.phases",value.disconnection.fixed.phases);
            pp.queryarr("disconnection.fixed.time",value.disconnection.fixed.time);
            Util::Assert(INFO,TEST(value.disconnection.fixed.sitex.size() == value.disconnection.fixed.sitey.size()));
            Util::Assert(INFO,TEST(value.disconnection.fixed.sitex.size() == value.disconnection.fixed.phases.size()));
            Util::Assert(INFO,TEST(value.disconnection.fixed.sitex.size() == value.disconnection.fixed.time.size()));
            value.disconnection.fixed.done.resize(value.disconnection.fixed.sitex.size(),false);
        }
        else
        {
            value.disconnection.unif_dist = std::uniform_real_distribution<double>(0.0,1.0);
            value.disconnection.int_dist = std::uniform_int_distribution<int>(0,1);
            value.disconnection.rand_num_gen.seed(amrex::ParallelDescriptor::MyProc());
        }

        // Boundary conditions for Etas
        std::string bc_type = "constant";
        pp.query("bc.eta.type",bc_type);  // Type (constnat)
        if (bc_type == "constant")  value.mybc = new BC::Constant(value.number_of_grains,pp,"bc.eta");
        else if (bc_type == "step") value.mybc = new BC::Step(pp,"bc.eta");

        // Initial conditions for Etas
        std::string ic_type = "constant";
        pp.query("ic.type", ic_type); // IC Type
        if (ic_type == "constant")                 value.ic = new IC::Constant(value.geom,pp,"ic");
        else if (ic_type == "perturbed_interface") value.ic = new IC::PerturbedInterface(value.geom,pp,"ic");
        else if (ic_type == "tabulated_interface") value.ic = new IC::TabulatedInterface(value.geom,pp,"ic");
        else if (ic_type == "voronoi")             value.ic = new IC::Voronoi(value.geom, pp, "ic.voronoi");
        else if (ic_type == "expression")          value.ic = new IC::Expression(value.geom,pp,"ic.expression");
        else if (ic_type == "sphere")              value.ic = new IC::Sphere(value.geom,pp,"ic.sphere"); 
        else Util::Abort(INFO, "No valid initial condition specified");


        value.RegisterNewFab(value.eta_new_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta",true);
        value.RegisterNewFab(value.eta_old_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta old",false);
        value.RegisterNewFab(value.totaldf_mf, new BC::Nothing(), value.number_of_grains,0, "totaldf",true);
        value.RegisterNewFab(value.disc_mf, new BC::Nothing(), 1, value.number_of_ghost_cells, "disc",true);  // see box
        value.RegisterNodalFab(value.elasticdf_mf,value.number_of_grains,value.number_of_ghost_cells,"elasticdf",true);

        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.area, "area");
        value.RegisterIntegratedVariable(&value.gbenergy, "gbenergy");
        value.RegisterIntegratedVariable(&value.realgbenergy, "realgbenergy");
        value.RegisterIntegratedVariable(&value.regenergy, "regenergy");
        //value.RegisterIntegratedVariable(&elastic.strainenergy, "strainenergy");
    }

    
protected:

    /// \fn    Advance
    /// \brief Evolve phase field in time
    void Advance (int lev, Real time, Real dt) override;
    void Initialize (int lev) override;

    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void TimeStepBegin(amrex::Real time, int iter) override;
    void TimeStepComplete(amrex::Real time, int iter) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;


    virtual void UpdateModel(int /*a_step*/);

private:

    int number_of_grains = 2;
    int number_of_ghost_cells = 3;
    Set::Scalar ref_threshold = 0.1;

    // Cell fab
    Set::Field<Set::Scalar> eta_new_mf; // Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
    Set::Field<Set::Scalar> eta_old_mf; // Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
    Set::Field<Set::Scalar> fluct_mf;
    Set::Field<Set::Scalar> disc_mf; //see box
    // Node fab
    Set::Field<Set::Scalar> elasticdf_mf;
    Set::Field<Set::Scalar> totaldf_mf;

    BC::BC<Set::Scalar> *mybc;

    //amrex::Real M, mu, gamma, sigma0, l_gb, beta;
    RegularizationType regularization = RegularizationType::K12;
    struct {
        Set::Scalar M;
        Set::Scalar L;
        Set::Scalar mu;
        Set::Scalar gamma;
        Set::Scalar sigma0;
        Set::Scalar l_gb;
        Set::Scalar elastic_mult = 1.0;
        Set::Scalar elastic_threshold = 0.0;
    } pf;

    struct {
        int on = 0;
        Set::Scalar beta;
        Set::Scalar timestep;
        Set::Scalar tstart;
        int plot_int = -1;
        Set::Scalar plot_dt = -1.0;
        int thermo_int = -1, thermo_plot_int = -1;
        Set::Scalar theta0,sigma0,sigma1;
        Set::Scalar phi0 = 0.0;
        int elastic_int = -1;
    } anisotropy;
    
    struct { 
        int on = 0;
        Set::Scalar tstart = 0.0;
        Set::Scalar vol0 = 0.0;
        Set::Scalar lambda = 0.0;
    } lagrange;

    struct {
        int on = 0;
        std::vector<Numeric::Interpolator::Linear<Set::Scalar>> val;
        Set::Scalar tstart = 0.0;
    } sdf;

    struct {
        int on = 0;
        Set::Scalar amp = 0.0;
        Set::Scalar sd = 0.0;
        Set::Scalar tstart = 0.0;
        std::normal_distribution<double> norm_dist;
        std::default_random_engine rand_num_gen;
    } fluctuation;

    struct {
        int on = 0;
        //int phase = 0;
        int interval = -1;
        bool nucleate = false;
        Set::Scalar tstart = 0.0;
        Set::Scalar fixed_site = 4.0;
        Set::Scalar range = 0.06;
        Set::Scalar nucleation_energy = 0.0;
        Set::Scalar tau_vol = 1.0;
        Set::Scalar temp = 0.0;
        Set::Scalar K_b = 8.617333262145e-5; // eV/K
        Set::Scalar box_size = 0.0;
        Set::Scalar epsilon = 1E-20;
        std::uniform_real_distribution<double> unif_dist;
        std::uniform_int_distribution<int> int_dist;
        std::default_random_engine rand_num_gen;
        std::vector<Set::Scalar> sitex;
        std::vector<Set::Scalar> sitey;
        std::vector<int> phases;
        struct {
            int on = 0;
            std::vector<Set::Scalar> sitex;
            std::vector<Set::Scalar> sitey;
            std::vector<int> phases;
            std::vector<Set::Scalar> time;
            std::vector<bool> done;
        } fixed;
    } disconnection;

    Model::Interface::GB::GB *boundary;

    IC::IC *ic;

    Set::Scalar volume = 5;
    Set::Scalar area = 0.0;
    Set::Scalar gbenergy = 0.0;
    Set::Scalar realgbenergy = 0.0;
    Set::Scalar regenergy = 0.0;

    struct
    {
        Set::Scalar tstart = 0.0;
        model_type model;
    } mechanics;
    

};
}
#endif
