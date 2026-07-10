#ifndef INTEGRATOR_LOWMACH_H
#define INTEGRATOR_LOWMACH_H

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "BC/Nothing.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/IC.H"
#include "Model/Gas/Gas.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Set/Base.H"

namespace Integrator
{
class LowMach : virtual public Integrator
{
public:
    static constexpr const char* name = "lowmach";

    LowMach() {};
    LowMach(IO::ParmParse& pp);
    static void Parse(LowMach& value, IO::ParmParse& pp);

    virtual ~LowMach()
    {
        delete velocity_bc;
        delete temperature_bc;
        delete mass_fraction_bc;
        delete pressure_bc;
        delete eta_bc;
        delete xi_bc;

        delete velocity_ic;
        delete temperature_ic;
        delete mass_fraction_ic;
        delete pressure_ic;
        delete eta_ic;
    }

protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int iter) override;
    void PrintDiagnostics(Set::Scalar time, int iter);
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void SanitizeEta(amrex::MultiFab& eta_mf);
    void FillReferenceMapBoundary(int lev, amrex::MultiFab& xi_mf);
    void InitializeReferenceMap(int lev, amrex::MultiFab& xi_mf);
    void RebuildReferenceMapOutsideEta(int lev, const amrex::MultiFab& eta_mf, amrex::MultiFab& xi_mf, Set::Scalar time);
    void UpdateSolidStress(int lev,
                           const amrex::MultiFab& u_mf,
                           const amrex::MultiFab& T_mf,
                           const amrex::MultiFab& eta_mf,
                           const amrex::MultiFab& xi_mf);
    void UpdateDerived(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf, const amrex::MultiFab& Y_mf);
    void EnforceStateBounds(amrex::MultiFab& u_mf, amrex::MultiFab& T_mf, amrex::MultiFab& Y_mf);
    void ProjectVelocity(int lev, Set::Scalar time, Set::Scalar dt);
    void CompositeProjectVelocity(Set::Scalar time, Set::Scalar dt);
    void ImplicitElasticVelocitySolve(int lev, Set::Scalar time, Set::Scalar dt);
    void RHS(int lev, Set::Scalar time,
             amrex::MultiFab& u_rhs_mf,
             amrex::MultiFab& T_rhs_mf,
             amrex::MultiFab& Y_rhs_mf,
             amrex::MultiFab& eta_rhs_mf,
             amrex::MultiFab& xi_rhs_mf,
             const amrex::MultiFab& u_mf,
             const amrex::MultiFab& T_mf,
             const amrex::MultiFab& Y_mf,
             const amrex::MultiFab& eta_mf,
             const amrex::MultiFab& xi_mf);

    Set::Field<Set::Scalar> velocity_mf;
    Set::Field<Set::Scalar> velocity_old_mf;
    Set::Field<Set::Scalar> temperature_mf;
    Set::Field<Set::Scalar> temperature_old_mf;
    Set::Field<Set::Scalar> mass_fraction_mf;
    Set::Field<Set::Scalar> mass_fraction_old_mf;
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> xi_mf;
    Set::Field<Set::Scalar> xi_old_mf;

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> momentum_mf;
    Set::Field<Set::Scalar> pressure_mf;
    Set::Field<Set::Scalar> pressure_correction_mf;
    Set::Field<Set::Scalar> projection_rhs_mf;
    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> mole_fraction_mf;
    Set::Field<Set::Scalar> vorticity_mf;
    Set::Field<Set::Scalar> deformation_gradient_mf;
    Set::Field<Set::Scalar> piola_stress_mf;
    Set::Field<Set::Scalar> elastic_stress_mf;
    Set::Field<Set::Scalar> solid_stress_mf;
    Set::Field<Set::Matrix> cauchy_stress_mf;
    Set::Field<Set::Matrix> stress_rhs_mf;

    BC::BC<Set::Scalar>* velocity_bc = nullptr;
    BC::BC<Set::Scalar>* temperature_bc = nullptr;
    BC::BC<Set::Scalar>* mass_fraction_bc = nullptr;
    BC::BC<Set::Scalar>* pressure_bc = nullptr;
    BC::BC<Set::Scalar>* eta_bc = nullptr;
    BC::BC<Set::Scalar>* xi_bc = nullptr;
    BC::Nothing bc_nothing;

    IC::IC<Set::Scalar>* velocity_ic = nullptr;
    IC::IC<Set::Scalar>* temperature_ic = nullptr;
    IC::IC<Set::Scalar>* mass_fraction_ic = nullptr;
    IC::IC<Set::Scalar>* pressure_ic = nullptr;
    IC::IC<Set::Scalar>* eta_ic = nullptr;

    Model::Gas::Gas gas;
    Model::Solid::Finite::NeoHookean finite_solid_model;
    int nspecies = 0;

    Set::Vector g = Set::Vector::Zero();
    Set::Scalar cfl = NAN;
    Set::Scalar cfl_v = 1.0e100;
    Set::Scalar small = 1.0e-12;
    Set::Scalar density_floor = 1.0e-12;
    Set::Scalar temperature_floor = 1.0e-12;
    Set::Scalar pressure_floor = 1.0e-12;
    Set::Scalar pressure_scale = 1.0;
    Set::Scalar projection_tol_rel = 1.0e-11;
    Set::Scalar projection_tol_abs = 1.0e-12;
    Set::Scalar velocity_refinement_criterion = NAN;
    Set::Scalar pressure_refinement_criterion = NAN;
    Set::Scalar temperature_refinement_criterion = NAN;
    Set::Scalar eta_refinement_criterion = NAN;
    int projection_verbose = 0;
    int projection_node_reconstruction_sweeps = 1;
    int projection_nodal_iterations = 200;
    int diagnostics_interval = 0;
    bool projection_enabled = true;
    bool projection_amr_enabled = false;
    bool projection_update_pressure = false;
    bool include_viscosity = true;
    bool include_conduction = true;
    bool advect_temperature = true;
    Set::Scalar eta_initial_value = 0.0;
    Set::Scalar reference_map_eta_cutoff = 0.5;
    Set::Scalar reference_map_eta_core = 0.5;
    Set::Scalar reference_map_eta_extension = 1.0e-3;
    int reference_map_extrapolation_sweeps = 4;
    int reference_map_smoothing_sweeps = 2;
    bool finite_solid_enabled = false;
    Set::Scalar finite_solid_eta_threshold = 0.5;
    Set::Scalar finite_solid_J_floor = 1.0e-6;
    Set::Scalar finite_solid_viscosity = 0.0;
    Set::Scalar finite_solid_bulk_viscosity = 0.0;
    Set::Scalar finite_solid_stress_cap = 1.0e100;
    Set::Scalar finite_solid_rhs_sign = 0.0;
    bool finite_solid_implicit = true;
    Set::Scalar finite_solid_implicit_coeff_scale = 1.0e-4;
    bool finite_solid_implicit_preserve_rigid_modes = true;
    Set::Scalar finite_solid_implicit_tol_rel = 1.0e-8;
    Set::Scalar finite_solid_implicit_tol_abs = 1.0e-10;
    int finite_solid_implicit_verbose = 0;
};
}

#endif
