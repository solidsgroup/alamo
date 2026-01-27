// This integrator implements a phase field model that simulates
// agglomeration of of a nanopowder additive in a solid composite
// propellant (SCP). The chemical potential :math:`W_C(\alpha)` is
// defined as
//
// .. math::
//    W_C(\alpha)=\frac{\sigma}{\epsilon}\alpha^2(1-\alpha)^2
//
// where :math:`\sigma` is the interface energy and
// :math:`\varepsilon` is the diffuse interface length parameter. To
// capture the surface tension of the agglomerate phase, we define
//
// .. math::
//    W_B(\alpha)=\frac{1}{2}\sigma\epsilon|\nabla\alpha|^2.
//
// The free energy is
//
// .. math::
//    F[\alpha]=\int_\Omega(W_c(\alpha)+W_B(\alpha))\mathrm{d}\boldsymbol{x}=\int_\Omega\left(\frac{\sigma}{\epsilon}\alpha^2(1-\alpha)^2+\frac{1}{2}\sigma\epsilon|\nabla\alpha|^2\right)\mathrm{d}\boldsymbol{x}.
//
// The kinetic coefficient ("mobility") :math:`L` governing the
// evolution of :math:`\alpha` is defined as
//
// .. math::
//    L=L_0(1-\eta)^n
//
// where :math:`L_0` is the free-flow parameter and :math:`n>0` is a
// model constant. This definition restricts the movement of
// agglomerate in the binder.
//
// The kinetics are modeled with Allen-Cahn kinetics with a Lagrangian
// constraint on the volume of agglomerate in the domain:
//
// .. math::
//    \frac{\partial\alpha}{\partial t}=-L\left[\frac{\partial F}{\partial\alpha}+\lambda\left(\frac{V}{V_0}-1\right)\right]
//
// where :math:`\lambda` is a Lagrangian multiplier that is suitably large to enforce the constraint up to tolerance, :math:`V` is the volume of agglomerate in the domain, and :math:`V_0` is the prescribed volume of agglomerate in the domain. :math:`V_0` is usually equal to the initial volume of agglomerate in the domain.
//
// For reference, the functional derivative of the free energy :math:`F` is
//
// .. math::
//    \frac{\delta F}{\delta\alpha}=2\frac{\sigma}{\epsilon}\alpha(1-\alpha)(1-2\alpha)-\sigma\epsilon\Delta\alpha
//
// where :math:`\Delta\alpha` is the Laplacian of :math:`\alpha`.
//
// See the following reference for more details on the model:
//
// TODO: Add reference to 2026 AIAA SciTech paper once published.

#ifndef INTEGRATOR_AGGLOMERATION_H
#define INTEGRATOR_AGGLOMERATION_H

#include "BC/BC.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Integrator/Flame.H"
#include "Set/Base.H"

#include "AMReX_TagBox.H"

namespace Integrator
{

class Agglomeration : public Flame
{
public:
    static constexpr const char *name = "agglomeration";

    Agglomeration() = default;
    ~Agglomeration() = default;
    Agglomeration(IO::ParmParse &pp) : Agglomeration() { Parse(*this, pp); };

    static void Parse(Agglomeration &value, IO::ParmParse &pp);

private:
    Set::Scalar CalculateInitialVolume(int target_resolution_level);
    void ScaleByComplement(MultiFab &dst, const MultiFab &src, int srccomp, int dstcomp, int numcomp, int nghost);
    void ScaleByPhiComplement(const Set::Field<Set::Scalar> &mf, int lev);
    void Initialize(int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;

    struct
    {
        Set::Field<Set::Scalar> alpha_0;
        Set::Field<Set::Scalar> alpha_old;
        Set::Field<Set::Scalar> alpha;
        Set::Field<Set::Scalar> alphadot_agglom;
        Set::Field<Set::Scalar> alphadot_reaction;
        IC::IC<Set::Scalar> *alpha_ic = nullptr;
        BC::BC<Set::Scalar> *alpha_bc = nullptr;
        Set::Scalar sigma = NAN;
        Set::Scalar epsilon = NAN;
        Set::Scalar lambda = NAN;
        Set::Scalar L_0_agglom = NAN;
        Set::Scalar L_0_reaction = NAN;
        Set::Scalar n = NAN;
        Set::Scalar V = NAN;
        Set::Scalar V_0 = NAN;
        bool calculate_initial_volume = true;
        Set::Scalar eta_refinement_threshold = NAN;
        Set::Scalar phi_refinement_threshold = NAN;
        Set::Scalar gradient_alpha_refinement_threshold = NAN;
    } agglom;
};
}
#endif
