#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Constant.H"
#include "IC/Voronoi.H"
#include "IC/IC.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/PSRead.H"
#include "IC/Laminate.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"



namespace Integrator
{
class Flame : public Integrator::Integrator
{

    using model_type = Model::Solid::Affine::Isotropic;

public:
    Flame ();
    Flame (IO::ParmParse &pp) : Flame()
    {Parse(*this,pp);}
    static void Parse(Flame &value, IO::ParmParse &pp)
    {
        {
            // These are the phase field method parameters
            // that you use to inform the phase field method.
            pp.query("pf.eps",    value.pf.eps); // Burn width thickness
            pp.query("pf.kappa",  value.pf.kappa); // Interface energy param
            pp.query("pf.gamma",  value.pf.gamma); // Scaling factor for mobility
            pp.query("pf.lambda", value.pf.lambda); // Chemical potential multiplier
            pp.query("pf.w1",     value.pf.w1); // Unburned rest energy
            pp.query("pf.w12",    value.pf.w12);  // Barrier energy
            pp.query("pf.w0",     value.pf.w0);    // Burned rest energy
            pp.query("pf.P",      value.pf.P);             // Pressure [UNITS?]
            pp.query("pf.r_ap",   value.pf.r_ap);       // AP Power law multiplier
            pp.query("pf.n_ap",   value.pf.n_ap);       // AP Power law exponent
            pp.query("pf.r_htpb", value.pf.r_htpb);   // HTPB Power law multiplier
            pp.query("pf.n_htpb", value.pf.n_htpb);   // HTPB Power law exponent
            pp.query("pf.r_comb", value.pf.r_comb);   // Combination power law multiplier
            pp.query("pf.n_comb", value.pf.n_comb);   // Combination power law exponent

            value.EtaBC = new BC::Constant(1);
            pp.queryclass("pf.eta.bc", *static_cast<BC::Constant *>(value.EtaBC)); // See :ref:`BC::Constant`

            value.RegisterNewFab(value.Eta_mf, value.EtaBC, 1, 1, "Eta", true);
            value.RegisterNewFab(value.Eta_old_mf, value.EtaBC, 1, 1, "Eta_old", false);
        }

        {
            // Read in parameters to determine the IC for eta
            std::string type = "constant";
            pp.query("eta.ic.type", type); // IC type - [packedspheres,laminate] - see classes for more information
            if (type == "constant") value.EtaIC = new IC::Constant(value.geom,pp,"eta.ic.constant");
            else if (type == "expression") value.EtaIC  = new IC::Expression(value.geom,pp,"eta.ic.expression");
            else Util::Abort(INFO,"Invalid eta.ic ", type);
        }

        {
            // These parameters are for the **Thermal transport model**
            pp.query("thermal.on",   value.thermal.on);       // Whether to use the thermal model
            pp.query("thermal.rho1", value.thermal.rho1); // Density (before)
            pp.query("thermal.rho0", value.thermal.rho0); // Density (after)
            pp.query("thermal.ka",   value.thermal.ka); // Thermal conductivity (before and after)
            pp.query("thermal.kh",   value.thermal.kh); // Thermal conductivity (before and after)
            pp.query("thermal.k0",   value.thermal.k0); // Thermal conductivity (before and after)
            pp.query("thermal.cp1",  value.thermal.cp1); // Specific heat (before and after)
            pp.query("thermal.cp0",  value.thermal.cp0); // Specific heat (before and after)
            pp.query("thermal.delA", value.thermal.delA); // Thermal flux of each material
            pp.query("thermal.delH", value.thermal.delH); // Thermal flux of each material
            if (value.thermal.on)
            {
                value.TempBC = new BC::Constant(1);
                pp.queryclass("thermal.bc", *static_cast<BC::Constant *>(value.TempBC)); // See :ref:`BC::Constant`
                value.RegisterNewFab(value.Temp_mf, value.TempBC, 1, 1, "Temp", true);
                value.RegisterNewFab(value.Temp_old_mf, value.TempBC, 1, 1, "Temp_old", false);
            }
        }

        {
            pp.query("amr.refinement_criterion", value.m_refinement_criterion); // Refinement criterion for eta field
            pp.query("amr.refinement_criterion_temp", value.t_refinement_criterion); // Refinement criterion for temperature field
        }

        std::vector<Set::Scalar> fs(value.fs_number, 1);
        
        {
            // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
            // specified using these parameters. 
            //IO::ParmParse pp("phi.ic");
            std::string type = "packedspheres";
            pp.query("phi.ic.type", type); // IC type - [packedspheres,laminate] - see classes for more information
            if (type == "psread") value.PhiIC = new IC::PSRead(value.geom,pp,"phi.ic.psread"); 
            else if (type == "laminate") value.PhiIC = new IC::Laminate(value.geom,pp,"phi.ic.laminate");
            else if (type == "constant") value.PhiIC = new IC::Constant(value.geom,pp,"phi.ic.constant");
            else Util::Abort(INFO,"Invalid phi.ic.type = ", type);
            
            value.RegisterNewFab(value.phi_mf, value.EtaBC, 1, 1, "phi", true);
        }

        //
        // Register integrated variables
        //
        value.RegisterIntegratedVariable(&value.volume,"volume");
        value.RegisterIntegratedVariable(&value.area,"area");
    }


protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void Advance (int lev, amrex::Real time, amrex::Real dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;
private:

    Set::Field<Set::Scalar> Temp_mf;
    Set::Field<Set::Scalar> Temp_old_mf;
    Set::Field<Set::Scalar> Eta_mf;
    Set::Field<Set::Scalar> Eta_old_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    BC::BC<Set::Scalar> *TempBC;
    BC::BC<Set::Scalar> *EtaBC;
    IC::IC *EtaIC;
    IC::IC *PhiIC;


    int fs_number = 1;
 
   
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;
    
    struct {
       // Set::Scalar M = 0.0;
        Set::Scalar gamma = 1.0;
        Set::Scalar eps = 0.0;
        Set::Scalar lambda = 0.0;
        Set::Scalar kappa = 0.0;
        Set::Scalar w1 = 0.0, w12 = 0.0, w0 = 0.0;
        Set::Scalar P = 0.0;
        Set::Scalar r_ap = 0.0;
        Set::Scalar n_ap = 0.0;
        Set::Scalar r_htpb = 0.0;
        Set::Scalar n_htpb = 0.0;
        Set::Scalar r_comb = 0.0;
        Set::Scalar n_comb = 0.0;
    } pf;

    struct {
        int on = 0;
        Set::Scalar rho1, rho0;         
        Set::Scalar ka, kh, k0;             
        Set::Scalar cp1, cp0;           
        Set::Scalar delA, delH;       
    } thermal;

    struct {
        Set::Field<Set::Scalar> disp_mf;
        Set::Field<Set::Scalar> rhs_mf;
        Set::Field<Set::Scalar> res_mf;
        Set::Field<Set::Scalar> stress_mf;
        Set::Field<model_type>  model_mf;
        model_type model_ap, model_htpb, model_void;
        int interval = 0;
        BC::Operator::Elastic::Constant bc;
        Solver::Nonlocal::Newton<model_type> *solver;//(elastic.op);
    } elastic;


    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;
};
}

#endif
