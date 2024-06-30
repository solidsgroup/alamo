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

class ChemReduction : virtual public Integrator::Integrator,
{


public:
    ChemReduction();
    ChemReduction(IO::ParmParse& pp);

    static void Parse(ChemReduction& value, IO::ParmParse& pp);

    virtual ~ChemReduction()
    {
        delete bc_eta;
        delete ic_eta;
        delete ic_phi;
    }

protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar a_time, int a_iter) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override;
    void UpdateModel(int a_step, Set::Scalar a_time) override;
private:

    Set::Field<Set::Scalar> phi_o_mf;
    Set::Field<Set::Scalar> phi_h_mf;
    Set::Field<Set::Scalar> phi_fe_mf;

    Set::Field<Set::Scalar> eta_o_mf;
    Set::Field<Set::Scalar> eta_h_mf;
    Set::Field<Set::Scalar> eta_fe_mf;

    BC::BC<Set::Scalar>* bc_temp = nullptr;
    BC::BC<Set::Scalar>* bc_eta = nullptr;
    IC::IC* ic_phi = nullptr;
    IC::IC* ic_laser = nullptr;
    IC::IC* ic_eta = nullptr;

    int ghost_count = 2;

    Set::Scalar EnergyBarrier; 
    Set::Scalar GradientEnergy; 
    Set::Scalar ReactionConstant;
    Set::Scalar ChemicalGradientEnergy;
    Set::Scalar MolarVolume; 
    
    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::scalar D; // Diffusivity
        Set::Scalar L; // Mobility
        Set::Scalar M; // Mobility
    } fe;

    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::scalar D; // Diffusivity
        Set::Scalar L; // Mobility
        Set::Scalar M; // Mobility
    } h;

    struct {
        Set::Scalar G; // Shear Modulli 
        Set::Scalar K; // Bulk Modulli 
        Set::scalar D; // Diffusivity
        Set::Scalar L; // Mobility
        Set::Scalar M; // Mobility
    } o;
    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;
};
}

#endif
