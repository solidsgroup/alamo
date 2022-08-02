///
/// \file PolymerDegradation.H
///
#ifndef POLYMERDEGRADATION_H
#define POLYMERDEGRADATION_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Operator/Elastic/Constant.H"
#include "IC/IC.H"
#include "IC/Constant.H"
#include "Numeric/Interpolator/Linear.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"


#include "Operator/Elastic.H"
#include "Model/Solid/Linear/IsotropicDegradableTanh.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"


namespace Integrator
{
//using pd_model_type = Model::Solid::Linear::IsotropicDegradableTanh;
using pd_model_type = Model::Solid::Linear::IsotropicDegradable;

///
/// \class PolymerDegradation
/// \brief Phase field degradation model for hydro-thermal induced degradation of polymers.
///
/// Solve damage evolution laws for damage variable \f$ \eta \f$ or tensor \f$ \mathbf{F}_d \f$,
/// degrade material modulus tensor based on damage variable, and perform elasticity tests.
///
class PolymerDegradation : public Integrator::Integrator
{
public:
    PolymerDegradation();
    ~PolymerDegradation();

protected:

    /// \fn    Advance
    /// \brief Evolve phase field in time
    void Advance (int lev, Real time, Real dt);

    void Initialize (int lev);

    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

    void TimeStepBegin(amrex::Real time, int iter);

    void TimeStepComplete(amrex::Real time, int iter);

    void DegradeMaterial(int lev,amrex::FabArray<amrex::BaseFab<pd_model_type> > &model);

private:

    int number_of_ghost_cells = 3;
    int number_of_ghost_nodes = 2;

    // Degradation variable
    Set::Field<Set::Scalar> eta_new;         ///< Degradation variable for the __current__timestep
    Set::Field<Set::Scalar> eta_old;        ///< Degradation variable for the __previous__timestep

    // For water induced degradation
    Set::Field<Set::Scalar> water_conc;        ///< Water concentration for the __current__timestep
    Set::Field<Set::Scalar> water_conc_old;    ///< Water concentration for the __previous__timestep
    Set::Field<Set::Scalar> damage_start_time;

    // For temperature induced degradation
    Set::Field<Set::Scalar> Temp;            ///< Temperature for the __current__timestep
    Set::Field<Set::Scalar> Temp_old;        ///< Temperature for the __previous__timestep

    // For mechanical degradation to be studied later
    Set::Field<Set::Scalar>  displacement; ///< Multicomponent field variable storing pointwise displacements
    Set::Field<Set::Scalar>  strain; ///< Multicomponent field variable storing pointwise strains (6 components)
    Set::Field<Set::Scalar>  stress; ///< Multicomponent field variable storing pointwise strains (6 components)
    Set::Field<Set::Scalar>  stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
    Set::Field<Set::Scalar>  rhs; ///< Multicomponent field variable storing pointwise displacements
    Set::Field<Set::Scalar>  residual;
    Set::Field<Set::Scalar>  energy; ///< Single component variable storing strain energy
    Set::Field<Set::Scalar>  energies; ///< Multicomponent field variable storing individual strain energies

    int nlevels;

    std::string ic_type, gb_type, filename;


    amrex::Vector<amrex::BoxArray> ngrids,cgrids;
    amrex::Vector<amrex::DistributionMapping> ndmap;

    // Water diffusion parameters
    struct{
        bool             on                         =    false;
        Set::Scalar     diffusivity                =    1.0;
        Set::Scalar     refinement_threshold     =    0.01;
        std::string     ic_type;
        IC::IC            *ic;
        BC::BC<Set::Scalar>            *bc;
    } water;

    // Thermal diffusion parameters
    struct{
        bool            on                         =    false;
        Set::Scalar     diffusivity             =    1.0;
        Set::Scalar     refinement_threshold     =    0.01;
        std::string        ic_type;
        IC::IC            *ic;
        BC::BC<Set::Scalar>            *bc;
    } thermal;

    // Damage parameters
    struct{
        std::string                                    type;
        bool                                        anisotropy                     = false;
        int                                            number_of_eta                 = 1;
        amrex::Vector<Set::Scalar>                    d_final;
        amrex::Vector<int>                            number_of_terms;
        amrex::Vector<amrex::Vector<Set::Scalar> >     d_i;
        amrex::Vector<amrex::Vector<Set::Scalar> >     tau_i;
        amrex::Vector<amrex::Vector<Set::Scalar> >     t_start_i;
        Set::Scalar                                    refinement_threshold        = 0.01;
        std::string                                 ic_type;
        IC::IC                                        *ic;
        BC::BC<Set::Scalar>                                        *bc;
        BC::BC<Set::Scalar>                                        *bc_time;
    } damage;

    // Damage model: water relaxation
    struct{
        int                         number_of_terms    = 4;
        amrex::Vector<Set::Scalar>     d_i;
        Set::Scalar                 d_final         = 1.0;
        amrex::Vector<Set::Scalar>     tau_i;
        amrex::Vector<Set::Scalar>     t_start_i;
    } damage_w;

    // Damage model: temperature
    struct{
        Set::Scalar c0;
        Set::Scalar c1;
        Set::Scalar c2;
        Set::Scalar c3;
        Set::Scalar tau_T;
    } damage_T;

    // Elasticity parameters
    struct{
        Set::Scalar bottom_tol                 = 1.0E-6;
        bool        on                         = false;
        int         interval                 = 1;
        std::string type                    = "single";
        amrex::Vector<Set::Scalar> test_time;
        Set::Scalar test_rate                = 1.; 
        Set::Scalar test_duration            = 2.;
        Set::Scalar test_dt                    = 0.01;
        int            current_test            = 0;
        int         max_iter                 = 200;
        int         max_fmg_iter             = 0;
        int         bottom_max_iter            = 200;
        int         max_fixed_iter             = 5;
        int         verbose                 = 3;
        int         cgverbose                 = 3;
        Set::Scalar tol_rel                 = 1.0E-5;
        Set::Scalar tol_abs                 = 1.0E-5;
        Set::Scalar cg_tol_rel                 = 1.0E-5;
        Set::Scalar cg_tol_abs                 = 1.0E-5;
        Set::Scalar tstart                     = 0.0;
        Set::Scalar tend                     = 0.0;
        std::string bottom_solver            = "bicgstab";
        int         linop_maxorder              = 2;
        bool        use_fsmooth             = false;
        int         max_coarsening_level    = 0;
        bool         agglomeration               = true;
        bool         consolidation               = false;

        // Elastic BC
        std::array<BC::Operator::Elastic::Constant::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
        std::array<BC::Operator::Elastic::Constant::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

        amrex::Vector<Set::Scalar> body_force = {AMREX_D_DECL(0.0,0.0,0.0)};

        BC::Operator::Elastic::Constant bc;
        //Operator::Elastic<pd_model_type> op;
        std::map<std::string,BC::Operator::Elastic::Constant::Type >        bc_map;

    } elastic;

    struct{
        pd_model_type         modeltype;
        std::string     input_material = "isotropic";
        Set::Scalar        yieldstrength = 18.;
        Set::Field<pd_model_type> model;
    } material;

    std::vector<std::string> PlotFileNameNode (std::string plot_file_name, int lev) const;

    //Numeric::Interpolator::Linear<Set::Vector> AMREX_D_DECL(interpolate_left,interpolate_bottom,interpolate_back);
    //Numeric::Interpolator::Linear<Set::Vector> AMREX_D_DECL(interpolate_right,interpolate_top,interpolate_front);

};
}
#endif
