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
#include "Model/Chemistry/Chemistry.H"
#include "Model/Chemistry/FiniteRate.H"
#include "Model/Chemistry/Frozen.H"
#include "Model/Chemistry/Rocfire.H"
#include "Model/Chemistry/Rocfire_Aluminized.H"
#include "Model/Mechanism/Mechanism.H"
#include "Model/Mechanism/PhaseChange.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Set/Base.H"
#include "Numeric/Advect/Advect.H"
#include "Numeric/Advect/Centered.H"
#include "Numeric/Advect/MUSCL.H"
#include "Numeric/Advect/QUICK.H"
#include "Numeric/Advect/Upwind.H"
#include "Numeric/Advect/WENO5.H"
#include "Numeric/ReferenceMap/Reconstruction.H"
#include "Operator/Diffusion.H"
#include "Operator/PressurePoisson.H"

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
        delete component_density_bc;
        delete pressure_bc;
        delete xi_bc;

        delete velocity_ic;
        delete temperature_ic;
        for (auto* ic : component_density_ic)
            delete ic;
        delete pressure_ic;
        delete xi_ic;
    }

protected:
    void Initialize(int lev) override;
    void Regrid(int lev, Set::Scalar time) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int iter) override;
    void PreparePlotFile(Set::Scalar time, const amrex::Vector<int>& iter) override;
    void PrintDiagnostics(Set::Scalar time, int iter);
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void UpdateSolidStress(int lev,
                            const amrex::MultiFab& u_mf,
                            const amrex::MultiFab& eta_mf,
                            const amrex::MultiFab& xi_mf,
                            bool write_diagnostics = false);
    void UpdateComponentState(int lev, const amrex::MultiFab& component_density_mf);
    void UpdateDerivedDiagnostics(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf);
    void AdvanceChemistry(int lev, amrex::MultiFab& temperature_mf,
                            amrex::MultiFab& component_density_mf, Set::Scalar dt);
    void ApplyImplicitDiffusion(Set::Scalar time, Set::Scalar dt);
    void UpdateFreeRigidBodyState();

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::tuple<Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar>
    ComputeThermalState(
        Set::Patch<const Set::Scalar> component_density,
        Set::Scalar temperature, int i, int j, int k) const;

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::tuple<Model::Chemistry::SpeciesArray, Set::Scalar /*temperature*/, Set::Scalar /*dilatation*/>
    ComputeThermochemicalSource(
        Set::Patch<const Set::Scalar> component_density,
        Set::Patch<const Set::Scalar> temperature,
        int i, int j, int k, const Set::Scalar* dx, Set::Scalar dt,
        bool include_reaction = true) const;
    void ProjectVelocity(Set::Scalar time, Set::Scalar dt);
    void RHS(int lev, Set::Scalar time, Set::Scalar dt,
            amrex::MultiFab& u_rhs_mf,
            amrex::MultiFab& T_rhs_mf,
            amrex::MultiFab& component_density_rhs_mf,
            amrex::MultiFab* xi_rhs_mf,
            const amrex::MultiFab& u_mf,
            const amrex::MultiFab& T_mf,
            const amrex::MultiFab& component_density_mf,
            const amrex::MultiFab* xi_mf);

    Set::Field<Set::Scalar> velocity_mf;
    Set::Field<Set::Scalar> velocity_old_mf;
    Set::Field<Set::Scalar> temperature_mf;
    Set::Field<Set::Scalar> temperature_old_mf;
    Set::Field<Set::Scalar> component_density_mf;
    Set::Field<Set::Scalar> component_density_old_mf;
    Set::Field<Set::Scalar> mass_fraction_mf;
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> rigid_eta_mf;
    Set::Field<Set::Scalar> fixed_rigid_eta_mf;
    Set::Field<Set::Scalar> free_rigid_eta_mf;
    Set::Field<Set::Scalar> xi_mf;
    Set::Field<Set::Scalar> xi_old_mf;

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> momentum_mf;
    Set::Field<Set::Scalar> pressure_mf;
    Set::Field<Set::Scalar> pressure_correction_mf;
    Set::Field<Set::Scalar> chemistry_dilatation_mf;
    Set::Field<Set::Scalar> diffusion_dilatation_mf;
    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> mole_fraction_mf;
    Set::Field<Set::Scalar> vorticity_mf;
    Set::Field<Set::Scalar> viscosity_mf;
    Set::Field<Set::Scalar> thermal_conductivity_coeff_mf;
    Set::Field<Set::Scalar> diffusion_coeff_mf;
    Set::Field<Set::Scalar> wdot_mf;
    Set::Field<Set::Scalar> qdot_mf;
    Set::Field<Set::Matrix> deformation_gradient_mf;
    Set::Field<Set::Matrix> solid_deviatoric_stress_mf;

    BC::BC<Set::Scalar>* velocity_bc = nullptr;
    BC::BC<Set::Scalar>* temperature_bc = nullptr;
    BC::BC<Set::Scalar>* component_density_bc = nullptr;
    BC::BC<Set::Scalar>* pressure_bc = nullptr;
    BC::BC<Set::Scalar>* xi_bc = nullptr;
    BC::Nothing bc_nothing;

    IC::IC<Set::Scalar>* velocity_ic = nullptr;
    IC::IC<Set::Scalar>* temperature_ic = nullptr;
    std::vector<IC::IC<Set::Scalar>*> component_density_ic;
    IC::IC<Set::Scalar>* pressure_ic = nullptr;
    IC::IC<Set::Scalar>* xi_ic = nullptr;

    Model::Gas::Gas gas;
    using Chemistry = Model::Chemistry::Chemistry<
        Model::Chemistry::Frozen,
        Model::Chemistry::FiniteRate,
        Model::Chemistry::Rocfire,
        Model::Chemistry::Rocfire_Aluminized>;
    Chemistry chemistry;
    Model::Solid::Finite::NeoHookean finite_solid_model;
    using Mechanism = Model::Mechanism::Mechanism<Model::Mechanism::PhaseChange>;
    std::vector<Mechanism> mechanisms;
    std::vector<std::string> species_names;
    std::vector<Set::Scalar> reference_density;
    Model::Chemistry::SpeciesArray condensed_specific_heat;
    Model::Chemistry::SpeciesArray condensed_thermal_conductivity;
    Model::Chemistry::SpeciesArray condensed_inverse_reference_density;
    int ngas_species = 0;
    int nspecies = 0;
    int deformable_solid_species = -1;
    std::vector<int> rigid_solid_species;
    std::vector<int> fixed_rigid_solid_species;
    int free_rigid_solid_species = -1;

    Set::Vector g = Set::Vector::Zero();
    Set::Scalar cfl = NAN;
    Set::Scalar cfl_v = 1.0e100;
    Set::Scalar small = 1.0e-12;
    Set::Scalar density_floor = 1.0e-12;
    Set::Scalar pressure_floor = 1.0e-12;
    Set::Scalar pressure_scale = 1.0;
    Set::Scalar pressure_reference = NAN;
    Set::Scalar velocity_refinement_criterion = NAN;
    Set::Scalar pressure_refinement_criterion = NAN;
    Set::Scalar temperature_refinement_criterion = NAN;
    Set::Scalar eta_refinement_criterion = NAN;
    Set::Scalar reinitialize_condensed_composition_eta_min = 0.99;
    int diagnostics_interval = 0;
    bool diagnostics_extended_fields = false;
    bool reinitialize_condensed_composition = false;
    bool projection_enabled = true;
    bool projection_update_pressure = true;
    bool include_viscosity = true;
    bool include_conduction = true;
    bool implicit_momentum_diffusion = false;
    bool implicit_thermal_diffusion = false;
    bool implicit_species_diffusion = false;
    bool condensed_thermal_transport = false;
    bool advect_temperature = true;
    Set::Scalar finite_solid_reference_density = 1.0;
    Set::Scalar finite_solid_eta_threshold = 0.5;
    Set::Scalar finite_solid_J_floor = 1.0e-6;
    Set::Scalar finite_solid_viscosity = 0.0;
    Set::Scalar finite_solid_interface_viscosity = 0.0;
    Set::Scalar finite_solid_deviatoric_stress_divergence_sign = 0.0;
    Set::Scalar rigid_relaxation_time = 1.0;
    Set::Vector rigid_velocity = Set::Vector::Zero();
    int free_rigid_picard_iterations = 1;
    Set::Scalar free_rigid_picard_tolerance = 1.0e-6;
    int last_free_rigid_picard_iterations = 0;
    Set::Scalar last_free_rigid_picard_residual = 0.0;
    Set::Scalar previous_free_rigid_dt = NAN;

    struct FreeRigidBodyState
    {
        Set::Scalar mass = 0.0;
        Set::Vector center = Set::Vector::Zero();
        Set::Vector velocity = Set::Vector::Zero();
        Set::Vector angular_velocity = Set::Vector::Zero();
        Set::Matrix inertia = Set::Matrix::Zero();
        Set::Scalar radius_of_gyration = 0.0;
        bool valid = false;
    };
    FreeRigidBodyState free_rigid_body;
    FreeRigidBodyState previous_free_rigid_body;

    Numeric::Advect::Advect<Numeric::Advect::MUSCL,
                            Numeric::Advect::Upwind,
                            Numeric::Advect::Centered,
                            Numeric::Advect::QUICK,
                            Numeric::Advect::WENO5> advect;
    Numeric::ReferenceMap::Reconstruction reference_map_reconstruction;
    Operator::Diffusion diffusion;
    Operator::PressurePoisson pressure_poisson;
};
}

#endif
