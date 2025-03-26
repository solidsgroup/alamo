///
/// \file PhaseFieldMicrostructure.H
///
#ifndef INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H
#define INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H

#include <random>

#include "AMReX_ParmParse.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Constant.H"
#include "IC/Constant.H"
#include "IC/PerturbedInterface.H"
#include "IC/Voronoi.H"
#include "IC/Sphere.H"
#include "IC/Expression.H"
#include "IC/Ellipse.H"
#include "IC/Random.H"
#include "Model/Interface/GB/GB.H"
#include "Model/Interface/GB/Sin.H"
#include "Model/Interface/GB/AbsSin.H"
#include "Model/Interface/GB/Read.H"
#include "Model/Interface/GB/SH.H"
#include "Model/Defect/Disconnection.H"
#include "Base/Mechanics.H"

namespace Integrator
{

enum RegularizationType{
    Wilmore = 1,
    K12 = 2};

///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
template<class model_type>
class PhaseFieldMicrostructure : public Base::Mechanics<model_type>
{
public:
    const std::string name = "phasefieldmicrostructure." + std::string(model_type::name);

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
    
        pp_query_required("amr.max_level", value.max_level);            // Maximum AMR level
        pp_query_default("amr.ref_threshold", value.ref_threshold, 0.1);    // Phase field refinement threshold


        std::string type_str;
        // Reading this is redundant but necessary because of the way the code was originally structured
        // (need to fix eventually)
        pp.query_validate("mechanics.type", type_str, {"disable","static","dynamic"}); 
        value.m_type = Base::Mechanics<model_type>::Disable; // Turn mechanics off by default

        if (type_str != "disable") // do this only if mechanics is activated
        {
            // Elasticity
            pp_query_default("mechanics.tstart",value.mechanics.tstart, 0.0);

            // Read in models - either one model for all grains, or
            // individual models, one for each grain.
            pp.queryclass_enumerate("mechanics.model",value.mechanics.model, value.number_of_grains);
        
            // Mixing order 
            pp_query_validate("mechanics.mix_order",value.mechanics.model_mix_order,{1,2});
            // Force Neumann BCs on the model
            pp_query_default("mechanics.model_neuman_boundary",value.mechanics.model_neumann_boundary,false);

            pp.queryclass<Base::Mechanics<model_type>>("mechanics",value);
            if (value.m_type == Base::Mechanics<model_type>::Type::Static)
                value.number_of_ghost_cells = std::max(value.number_of_ghost_cells, 3);
        }


        // Lagrange multiplier method for enforcing volumes
        pp_query_default("lagrange.on", value.lagrange.on,false);
        if (value.lagrange.on)
        {
            pp_query_required("lagrange.lambda", value.lagrange.lambda);      // Lagrange multiplier value
            pp_query_required("lagrange.vol0", value.lagrange.vol0);          // Prescribed volume
            pp_query_default("lagrange.tstart", value.lagrange.tstart,0.0);   // Lagrange multipler start time
            value.SetThermoInt(1);
        }

        pp_query_default("sdf.on",value.sdf.on,false); // synthetic driving force (SDF)
        if(value.sdf.on)
        {
            std::vector<std::string> vals;
            pp_queryarr("sdf.val",vals);  // value of SDF for each grain
            int nvals = static_cast<int>(vals.size());
            if (nvals == 1)
                for (int i = 0; i < value.number_of_grains; i++)
                    value.sdf.val.push_back(Numeric::Interpolator::Linear<Set::Scalar>(vals[0]));
            else if (nvals == value.number_of_grains)
                for (int i = 0; i < value.number_of_grains; i++)
                    value.sdf.val.push_back(Numeric::Interpolator::Linear<Set::Scalar>(vals[i]));
            else
                Util::Abort(INFO,"sdf.val received ", vals.size(), " but requires 1 or ", value.number_of_grains);

            pp_query_default("sdf.tstart",value.sdf.tstart,0.0); // time to begin applying SDF
        }

        // Anisotropic grain boundary energy parameters

        pp_query_default("anisotropy.on", value.anisotropy.on,false);              // Turn on
        if (value.anisotropy.on)
        {
            // Regularization para m
            pp_query_required("anisotropy.beta", value.anisotropy.beta);
            // Time to turn on anisotropy
            pp_query_required("anisotropy.tstart", value.anisotropy.tstart);
            value.anisotropy.timestep = value.timestep;
            // Modify timestep when turned on
            pp_query_required("anisotropy.timestep", value.anisotropy.timestep);
            value.anisotropy.plot_int = value.plot_int;
            // Modify plot_int when turned on
            pp_query_default("anisotropy.plot_int", value.anisotropy.plot_int, -1);
            value.anisotropy.plot_dt = value.plot_dt;
            // Modify plot_dt when turned on
            pp_query_default("anisotropy.plot_dt", value.anisotropy.plot_dt, -1.0);
            // Modify thermo int when turned on
            pp_query_default("anisotropy.thermo_int", value.anisotropy.thermo_int, -1);
            // Modify thermo plot int when turned on
            pp_query_default("anisotropy.thermo_plot_int", value.anisotropy.thermo_plot_int, -1);
            // Frequency of elastic calculation
            pp_query_default("anisotropy.elastic_int",value.anisotropy.elastic_int, -1);           
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

            pp_forbid("anisotropy.gb_type"," --> anisotropy.type");
            // Type of GB to use
            pp.select<Model::Interface::GB::AbsSin, Model::Interface::GB::Sin, Model::Interface::GB::Read, Model::Interface::GB::SH>("anisotropy",value.boundary); 
        }
        // Thermal fluctuations
        pp.query("fluctuation.on",value.fluctuation.on);
        if (value.fluctuation.on)
        {
            pp.query("fluctuation.amp",value.fluctuation.amp); // fluctuation amplitude
            pp.query("fluctuation.sd",value.fluctuation.sd); // fluctuation stadard deviation
            pp.query("fluctuation.tstart", value.fluctuation.tstart); // time to start applying fluctuation
            value.fluctuation.norm_dist = std::normal_distribution<double>(0.0,value.fluctuation.sd);
        }

        // Disconnection generation
        pp.query("disconnection.on",value.disconnection.on);
        if (value.disconnection.on)
        {
            // Read in nucleation parameters from disconnection class
            pp.queryclass<Model::Defect::Disconnection>("disconnection",value.disconnection.model);
        }

        // Shear coupling matrices
        pp_query_default("shearcouple.on",value.shearcouple.on,false);
        if (value.shearcouple.on)
        {
            Util::AssertException(INFO,TEST(value.m_time_evolving == true), " mechanics.time_evolving must be true when using shearcouple");
            Util::AssertException(INFO,TEST(value.mechanics.model_mix_order == 2), " mechanics.model.mix_order must be 2 when using shearcouple");
            
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

        // Boundary condition for eta
        pp.select<BC::Constant>("bc.eta",value.mybc,value.number_of_grains);

        // Initial condition for the order parameter eta
        pp.select<IC::Constant,IC::PerturbedInterface,IC::Voronoi,IC::Expression,IC::Sphere,IC::Ellipse,IC::Random>("ic",value.ic,value.geom);

        // Anisotropic mobility
        pp_query_default("anisotropic_kinetics.on",value.anisotropic_kinetics.on, 0);
        if (value.anisotropic_kinetics.on)
        {
            // simulation time when anisotropic kinetics is activated
            pp_query_default("anisotropic_kinetics.tstart",value.anisotropic_kinetics.tstart, 0.0);
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
        value.RegisterNewFab(value.eta_old_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "EtaOld",false);
        //value.RegisterNewFab(value.driving_force_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForce",false);
        // if (value.pf.threshold.on)
        //     value.RegisterNewFab(value.driving_force_threshold_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForceThreshold",false);
        // if (value.disconnection.on)
        //     value.RegisterNewFab(value.disc_mf, new BC::Nothing(), 1, value.number_of_ghost_cells, "disc",true);  // see box

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

    virtual void TimeStepBegin(amrex::Real time, int iter) override;
    virtual void TimeStepComplete(amrex::Real time, int iter) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;

    virtual void UpdateEigenstrain(int lev);
    virtual void UpdateEigenstrain()
    {
        for (int lev = 0; lev <= this->max_level; lev++)
            UpdateEigenstrain(lev);
    }

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
    //Set::Field<Set::Scalar> disc_mf; //see box
    // Node fab
    //Set::Field<Set::Scalar> elasticdf_mf;
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
        Set::Scalar elastic_mult = NAN;
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
        Model::Defect::Disconnection model;
    } disconnection;

    struct {
        int on = 0;
        std::vector<Set::Matrix> Fgb;
    } shearcouple;

    std::string gb_type, filename;

    Model::Interface::GB::GB *boundary = nullptr;

    IC::IC<Set::Scalar> *ic = nullptr;

    Set::Scalar volume = 5;
    Set::Scalar area = 0.0;
    Set::Scalar gbenergy = 0.0;
    Set::Scalar realgbenergy = 0.0;
    Set::Scalar regenergy = 0.0;

    struct
    {
        Set::Scalar tstart = 0.0;
        std::vector<model_type> model;
        int model_mix_order;
        bool model_neumann_boundary = false;
    } mechanics;
    
    using Base::Mechanics<model_type>::model_mf;
    using Base::Mechanics<model_type>::stress_mf;

};
}
#endif
