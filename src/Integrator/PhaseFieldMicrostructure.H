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
#include "IC/Ellipse.H"

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
//using model_type = Model::Solid::Affine::Cubic;

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
template<class model_type>
class PhaseFieldMicrostructure : public Base::Mechanics<model_type>
{
public:
    PhaseFieldMicrostructure() : Integrator() {};
    PhaseFieldMicrostructure(IO::ParmParse &pp) : Integrator() 
    {pp_queryclass(*this);}
    virtual ~PhaseFieldMicrostructure()
    {
        delete boundary;
        delete ic;
        delete mybc;
    }
    static void Parse(PhaseFieldMicrostructure &value, IO::ParmParse &pp)
    {
        BL_PROFILE("PhaseFieldMicrostructure::Parse");
        
        // Number of grain fields (may be more if using different IC)
        pp_query_default("pf.number_of_grains", value.number_of_grains,2);     
        pp_query_required("pf.M", value.pf.M);                         // Mobility
        pp_query_required("pf.mu", value.pf.mu);                       // Phase field :math:`\mu`
        pp_query_required("pf.gamma", value.pf.gamma);                 // Phase field :math:`\gamma`
        pp_query_required("pf.sigma0", value.pf.sigma0);               // Initial GB energy if not using  anisotropy
        pp_query_required("pf.l_gb", value.pf.l_gb);                   // Mobility
        pp_query_default("pf.elastic_df",value.pf.elastic_df,false);   // Determine whether to use elastic driving force
        pp_query_default("pf.elastic_mult",value.pf.elastic_mult,1.0); // Multiplier of elastic energy

        pp_query_default("pf.threshold.value",value.pf.threshold.value,0.0);          // Value used for thresholding kinetic relation
        pp_query_default("pf.threshold.chempot",value.pf.threshold.chempot,false);    // Whether to include chemical potential in threshold
        pp_query_default("pf.threshold.boundary",value.pf.threshold.boundary,false);  // Whether to include boundary energy in threshold
        pp_query_default("pf.threshold.corner",value.pf.threshold.corner,false);      // Whether to include corner regularization in threshold
        pp_query_default("pf.threshold.lagrange",value.pf.threshold.lagrange,false);  // Whether to include lagrange multiplier in threshold
        pp_query_default("pf.threshold.mechanics",value.pf.threshold.mechanics,false);// Whether to include mechanical driving force in threshold
        pp_query_default("pf.threshold.sdf",value.pf.threshold.sdf,false);// Whether to include mechanical driving force in threshold
        {
            value.pf.threshold.on =
                value.pf.threshold.chempot || value.pf.threshold.boundary ||
                value.pf.threshold.corner || value.pf.threshold.lagrange ||
                value.pf.threshold.mechanics;

            std::string type = "continuous";
            pp_query_validate("pf.threshold.type",type,{"continuous","chop"}); // Type of thresholding to use
            if (type == "continuous") value.pf.threshold.type = ThresholdType::Continuous;
            else if (type == "chop") value.pf.threshold.type = ThresholdType::Chop;
        }
        
        value.pf.L = (4./3.)*value.pf.M / value.pf.l_gb;
    
        pp_query("amr.max_level", value.max_level);            // Maximum AMR level
        pp_query("amr.ref_threshold", value.ref_threshold);    // Phase field refinement threshold

        // Elasticity
        pp_query("mechanics.tstart",value.mechanics.tstart);
        value.mechanics.model.clear();
        value.mechanics.model.resize(value.number_of_grains,model_type::Zero());
        for (int i=0; i < value.number_of_grains; i++)
        {
            // By default, read in the model specified by "mechanics.model"
            if (pp.getEntries("mechanics.model").size())
            pp_queryclass("mechanics.model", value.mechanics.model[i]);

            // If individual models are specified, read those in and overwrite
            std::string name = "mechanics.model" + std::to_string(i+1);
            if (pp.getEntries(name).size())
                pp_queryclass(std::string(name.data()), value.mechanics.model[i]);
        }

        value.m_type = Base::Mechanics<model_type>::Disable; // Turn mechanics off by default
        pp_queryclass("mechanics",static_cast<Base::Mechanics<model_type>&>(value));
        if (value.m_type == Base::Mechanics<model_type>::Type::Static) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,3);


        // Lagrange multiplier method for enforcing volumes
        pp_query_default("lagrange.on", value.lagrange.on,false);
        if (value.lagrange.on)
        {
            pp_query_required("lagrange.lambda", value.lagrange.lambda);      // Lagrange multiplier value
            pp_query_required("lagrange.vol0", value.lagrange.vol0);          // Prescribed volume
            pp_query_default("lagrange.tstart", value.lagrange.tstart,0.0);   // Lagrange multipler start time
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

        pp_query_default("anisotropy.on", value.anisotropy.on,false);              // Turn on
        if (value.anisotropy.on)
        {
            pp_query_required("anisotropy.beta", value.anisotropy.beta);               // Regularization param 
            pp_query_required("anisotropy.tstart", value.anisotropy.tstart);           // Time to turn on anisotropy
        value.anisotropy.timestep = value.timestep;
        pp_query_required("anisotropy.timestep", value.anisotropy.timestep);       // Modify timestep when turned on
        value.anisotropy.plot_int = value.plot_int;
        pp_query("anisotropy.plot_int", value.anisotropy.plot_int);                // Modify plot_int when turned on
        value.anisotropy.plot_dt = value.plot_dt;
        pp_query("anisotropy.plot_dt", value.anisotropy.plot_dt);                  // Modify plot_dt when turned on
        pp_query("anisotropy.thermo_int", value.anisotropy.thermo_int);            // Modify thermo int when turned on
        pp_query("anisotropy.thermo_plot_int", value.anisotropy.thermo_plot_int);  // Modify thermo plot int when turned on
        pp_query("anisotropy.elastic_int",value.anisotropy.elastic_int);           // Frequency of elastic calculation
        if (value.anisotropy.on) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,2);

        // Determine the kind of regularization to use
        std::map<std::string, RegularizationType> regularization_type;
        regularization_type["wilmore"] = RegularizationType::Wilmore;
        regularization_type["k12"] = RegularizationType::K12;

            std::string regularization_type_input;
            // Type of regularization to use  
            pp_query_validate("anisotropy.regularization", regularization_type_input,{"k12","wilmore"});    
        value.regularization = regularization_type[regularization_type_input];

            // Type of GB to use
            pp_query_validate("anisotropy.gb_type", value.gb_type,{"abssin","sin","read","sh"}); 
        if      (value.gb_type == "abssin") value.boundary = new Model::Interface::GB::AbsSin(pp,"anisotropy.abssin");
        else if (value.gb_type == "sin")    value.boundary = new Model::Interface::GB::Sin(pp,"anisotropy.sin");
        else if (value.gb_type == "read")   value.boundary = new Model::Interface::GB::Read(pp,"anisotropy.read");
        else if (value.gb_type == "sh")     value.boundary = new Model::Interface::GB::SH(pp,"anisotropy.sh");
            else if (value.anisotropy.on)       Util::Abort(INFO,"A GB model must be specified");
        }

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

        // Shear coupling matrices
        pp_query_default("shearcouple.on",value.shearcouple.on,false);
        if (value.shearcouple.on)
        {
            value.shearcouple.Fgb.resize(value.number_of_grains * value.number_of_grains, Set::Matrix::Zero());
            for (int i = 0 ; i < value.number_of_grains ; i++)
                for (int j = 0 ; j < value.number_of_grains ; j++)
                {
                    std::string name    = "shearcouple.Fgb."+std::to_string(i)+"."+std::to_string(j);
                    std::string namerev = "shearcouple.Fgb."+std::to_string(j)+"."+std::to_string(i);
                    if ( i==j && pp.contains(name.data())) 
                        Util::Abort(INFO,"Cannot specify self FGB ", name);
                    if (pp.contains(name.data()) && pp.contains(namerev.data())) 
                        Util::Abort(INFO,"Cannot specify both ",name," and ",namerev);
                    if (pp.contains(name.data())) 
                    {
                        pp.queryarr(name.data(),value.shearcouple.Fgb[i*value.number_of_grains + j]);
                        value.shearcouple.Fgb[j*value.number_of_grains + i] = - value.shearcouple.Fgb[i*value.number_of_grains + j];
                    }
                }
        }

        std::string bc_type;
        // Type of boundary condition to use for eta
        pp_query_validate("bc.eta.type",bc_type,{"constant","step"});  
        if (bc_type == "constant")  value.mybc = new BC::Constant(value.number_of_grains,pp,"bc.eta");
        else if (bc_type == "step") value.mybc = new BC::Step(pp,"bc.eta");

        std::string ic_type;
        // Initial condition to use for eta
        pp_query_validate("ic.type", ic_type,{"constant","perturbed_interface","tabulated_interface","voronoi","expression","sphere","ellipse"});  
        // [switch ic.type] Initial conditions for Etas
        if (ic_type == "constant")                 value.ic = new IC::Constant(value.geom,pp,"ic");
        else if (ic_type == "perturbed_interface") value.ic = new IC::PerturbedInterface(value.geom,pp,"ic");
        else if (ic_type == "tabulated_interface") value.ic = new IC::TabulatedInterface(value.geom,pp,"ic");
        else if (ic_type == "voronoi")             value.ic = new IC::Voronoi(value.geom, pp, "ic.voronoi");
        else if (ic_type == "expression")          value.ic = new IC::Expression(value.geom,pp,"ic.expression");
        else if (ic_type == "sphere")              value.ic = new IC::Sphere(value.geom,pp,"ic.sphere"); 
        else if (ic_type == "ellipse")             value.ic = new IC::Ellipse(value.geom,pp,"ic.ellipse");
        else Util::Abort(INFO, "No valid initial condition specified");
        // [end switch]

        // Anisotropic mobility
        pp_query("anisotropic_kinetics.on",value.anisotropic_kinetics.on);
        if (value.anisotropic_kinetics.on)
        {
            // simulation time when anisotropic kinetics is activated
            pp_query("anisotropic_kinetics.tstart",value.anisotropic_kinetics.tstart);
            std::string mobility_filename, threshold_filename;
            // file containing anisotropic mobility data
            pp_query_file("anisotropic_kinetics.mobility",mobility_filename); 
            value.anisotropic_kinetics.mobility = Numeric::Interpolator::Linear<Set::Scalar>::Read(mobility_filename);
            // file containing anisotropic mobility data
            pp_query_file("anisotropic_kinetics.threshold",threshold_filename); 
            value.anisotropic_kinetics.threshold = Numeric::Interpolator::Linear<Set::Scalar>::Read(threshold_filename);
            value.RegisterNewFab(value.anisotropic_kinetics.L_mf, value.mybc, value.number_of_grains, 0, "mobility",true);
            value.RegisterNewFab(value.anisotropic_kinetics.threshold_mf, value.mybc, value.number_of_grains, 0, "theshold",true);
        }




        value.RegisterNewFab(value.eta_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta",true);
        if (value.shearcouple.on)
        {
            value.RegisterNewFab(value.eta_old_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "EtaOld", false);
            value.RegisterNodalFab(value.elasticdf_mf,value.number_of_grains,value.number_of_ghost_cells,"elasticdf",true);
        }
        value.RegisterNewFab(value.driving_force_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForce",false);
        if (value.pf.threshold.on)
            value.RegisterNewFab(value.driving_force_threshold_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForceThreshold",false);
        value.RegisterNewFab(value.disc_mf, new BC::Nothing(), 1, value.number_of_ghost_cells, "disc",true);  // see box

        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.area, "area");
        value.RegisterIntegratedVariable(&value.gbenergy, "gbenergy");
        value.RegisterIntegratedVariable(&value.realgbenergy, "realgbenergy");
        value.RegisterIntegratedVariable(&value.regenergy, "regenergy");


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

    void UpdateEigenstrain(int lev);

    virtual void UpdateModel(int /*a_step*/, Set::Scalar /*a_time*/) override;

private:

    int number_of_grains = -1;
    int number_of_ghost_cells = 1;
    Set::Scalar ref_threshold = 0.1;

    // Cell fab
    Set::Field<Set::Scalar> eta_mf; // Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> driving_force_mf;
    Set::Field<Set::Scalar> driving_force_threshold_mf;
    Set::Field<Set::Scalar> fluct_mf;
    Set::Field<Set::Scalar> disc_mf; //see box
    // Node fab
    Set::Field<Set::Scalar> elasticdf_mf;
    Set::Field<Set::Scalar> totaldf_mf;

    BC::BC<Set::Scalar> *mybc = nullptr;

    //amrex::Real M, mu, gamma, sigma0, l_gb, beta;
    RegularizationType regularization = RegularizationType::K12;
    enum ThresholdType {
        Continuous = 0, Chop = 1
    };
    struct {
        Set::Scalar M = NAN;
        Set::Scalar L = NAN;
        Set::Scalar mu = NAN;
        Set::Scalar gamma = NAN;
        Set::Scalar sigma0 = NAN;
        Set::Scalar l_gb = NAN;
        bool elastic_df = false;
        Set::Scalar elastic_mult = 1.0;
        struct {
            bool on = false;
            bool chempot = false;
            bool boundary = false;
            bool corner = false;
            bool lagrange = false;
            bool mechanics = false;
            bool sdf = false;
            Set::Scalar value = NAN;
            ThresholdType type = ThresholdType::Continuous;
        } threshold;
    } pf;

    struct {
        int on = 0;
        Set::Scalar tstart = 0.0;
        Numeric::Interpolator::Linear<Set::Scalar> mobility;
        Numeric::Interpolator::Linear<Set::Scalar> threshold;
        Set::Field<Set::Scalar> L_mf;
        Set::Field<Set::Scalar> threshold_mf;
    } anisotropic_kinetics;

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
        bool on = 0;
        Set::Scalar tstart = NAN;
        Set::Scalar vol0 = NAN;
        Set::Scalar lambda = NAN;
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

    struct {
        int on = 0;
        std::vector<Set::Matrix> Fgb;
    } shearcouple;

    std::string gb_type, filename;

    Model::Interface::GB::GB *boundary = nullptr;

    IC::IC *ic = nullptr;

    Set::Scalar volume = 5;
    Set::Scalar area = 0.0;
    Set::Scalar gbenergy = 0.0;
    Set::Scalar realgbenergy = 0.0;
    Set::Scalar regenergy = 0.0;

    struct
    {
        Set::Scalar tstart = 0.0;
        std::vector<model_type> model;
    } mechanics;
    

};
}
#endif
