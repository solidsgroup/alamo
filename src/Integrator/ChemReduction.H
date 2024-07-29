#ifndef INTEGRATOR_CHEMREDUCTION_H
#define INTEGRATOR_CHEMREDUCTION_H

#include <string>
#include <limits>
#include <memory>
#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"
#include "IC/Voronoi.H"
#include "IC/IC.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/PSRead.H" 
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Integrator
{

class ChemReduction : virtual public Integrator::Integrator
{
public:
    ChemReduction();
    ChemReduction(IO::ParmParse& pp);

    static void Parse(ChemReduction& value, IO::ParmParse& pp);

    virtual ~ChemReduction()
    {
        delete bc_phi;
        delete bc_chi;
        delete ic_chi;
        delete ic_phi;
    }

protected:
    void Initialize(int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
private:

    Set::Field<Set::Scalar> phi_f_mf; // Ferrite Field
    Set::Field<Set::Scalar> phi_w_mf; // Wustite Field
    Set::Field<Set::Scalar> phi_g_mf; // Water/Gas Field

    Set::Field<Set::Scalar> chi_o_mf; // Oxygen Field
    Set::Field<Set::Scalar> chi_h_mf; // Hydrogen Field
    Set::Field<Set::Scalar> chi_w_mf; // Water Field

    Set::Scalar phi_refinement_criterion;

    BC::BC<Set::Scalar>* bc_phi = nullptr;
    BC::BC<Set::Scalar>* bc_chi = nullptr;
    IC::IC* ic_phi = nullptr;
    IC::IC* ic_chi = nullptr;

    int ghost_count = 2; // Number of Ghost Cells

    Set::Scalar V; // Total Domain Molar Volume
    Set::Scalar R; // Universal Gas Constant
    Set::Scalar T; // Domain Temperature (Isothermal)

    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::Scalar L; // Diffusivity
        struct {
            Set::Scalar O; // Ferrite Oxygen Mobility
            Set::Scalar H; // Ferrite Hydrogen Mobility
            Set::Scalar W; // Ferrite Water Mobility
        } M; // Mobility
        struct {
            Set::Scalar O; // Oxygen Diffusivity
            Set::Scalar H; // Hydrogen Diffusivity
            Set::Scalar W; // Water Diffusivity
        } D;
        struct {
            Set::Scalar G; // Energy Gradient
            Set::Scalar B; // Energy Barrier
        } E;
    } ferrite;

    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::Scalar L; // Diffusivity
        struct {
            Set::Scalar O; // Wustite Oxygen Mobility
            Set::Scalar H; // Wustite Hydrogen Mobility
            Set::Scalar W; // Wustite Water Mobility
        } M; // Mobility
        struct {
            Set::Scalar O; // Oxygen Diffusivity
            Set::Scalar H; // Hydrogen Diffusivity
            Set::Scalar W; // Water Diffusivity
        } D;
        struct {
            Set::Scalar G; // Energy Gradient
            Set::Scalar B; // Energy Barrier
        } E;
    } wustite;

    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::Scalar L; // Diffusivity
        struct {
            Set::Scalar O; // Gas Oxygen Mobility
            Set::Scalar H; // Gas Hydrogen Mobility
            Set::Scalar W; // Gas Water Mobility
        } M; // Mobility
        struct {
            Set::Scalar O; // Oxygen Diffusivity
            Set::Scalar H; // Hydrogen Diffusivity
            Set::Scalar W; // Water Diffusivity
        } D; // Diffusivity
        struct {
            Set::Scalar G; // Energy Gradient
            Set::Scalar B; // Energy Barrier
        } E; // Energy
    } gas;
    // Integrated variables
    // Set::Scalar volume = 0.0;
    // Set::Scalar area = 0.0;
};
}
#endif
