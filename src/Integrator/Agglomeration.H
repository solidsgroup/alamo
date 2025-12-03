// TODO: Revise this to match the current model
//
// This integrator implements a phase field model for agglomeration of
// nanoparticles in a solid composite propellant (SCP). The chemical
// potential :math:`W_C(\alpha)` is defined as
//
// .. math::
//    W_C(\alpha)=\frac{\sigma}{\varepsilon}\alpha^2(1-\alpha)^2
//
// where :math:`\sigma` is the interface energy and
// :math:`\varepsilon` is the diffuse interface width. To capture the
// surface tension of the agglomerate phase, we define
//
// .. math::
//    W_B(\alpha)=\frac{1}{2}\sigma\varepsilon|\nabla\alpha|^2.
//
// The free energy is
//
// .. math::
//    F[\alpha]=\int_\Omega(W_c(\alpha)+W_B(\alpha))\mathrm{d}\boldsymbol{x}.
//
// The effective mobility :math:`L` of the agglomerate phase is
// defined as
//
// .. math::
//    L=L_0(1-\eta)^n
//
// where :math:`L_0` is the prescribed mobility parameter and
// :math:`n` is some constant greater than 0.
//
// The agglomerate phase can be evolved according to Cahn-Hilliard
// kinetics
//
// .. math::
//    \frac{\partial\alpha}{\partial t}=L\Delta\left(\frac{\delta F}{\delta\alpha}\right)
//
// or Allen-Cahn kinetics with a Lagrangian constraint
//
// .. math::
//    \frac{\partial\alpha}{\partial t}=-L\left(\frac{\partial F}{\partial\alpha}+\lambda\left(\frac{V}{V_0}-1\right)\right)
//
// where :math:`\lambda` is a suitably large Lagrangian multiplier to enforce the constraint up to tolerance, :math:`V` is the volume of agglomerate in the domain, and :math:`V_0` is the prescribed volume of agglomerate in the domain.
//
// For reference, the functional derivative of the free energy :math:`F` is
//
// .. math::
//    \frac{\delta F}{\delta\alpha}=2\frac{\sigma}{\varepsilon}\alpha(1-\alpha)(1-2\alpha)-\sigma\varepsilon\Delta\alpha.

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
    void Initialize(int lev) override;

    /**
     * @brief Scales the dst order parameter by the complement of src
     *
     * This function evaluates dst * (1 - src) and stores the result in dst.
     *
     * @param dst The order parameter to be scaled
     * @param src The scaling order parameter
     * @param srccomp The starting component in the src MultiFab
     * @param dstcomp The starting component in the dst MultiFab
     * @param numcomp The number of components to scale
     * @param nghost The number of ghost cells to scale
     *
     * @warning Assumes src is in the range [0.0,1.0].
     */
    void scaleByComplement(MultiFab &dst, const MultiFab &src, int srccomp, int dstcomp, int numcomp, int nghost);

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow) override;
    void Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mif, const amrex::Box &box) override;

    struct
    {
        Set::Field<Set::Scalar> alpha_old;
        Set::Field<Set::Scalar> alpha;
        IC::IC<Set::Scalar> *alpha_ic = nullptr;
        BC::BC<Set::Scalar> *alpha_bc = nullptr;
        Set::Scalar sigma = NAN;
        Set::Scalar epsilon = NAN;
        Set::Scalar lambda = NAN;
        Set::Scalar L_0 = NAN;
        Set::Scalar n = NAN;
        Set::Scalar V = NAN;
        bool calculate_initial_volume = true;
        Set::Scalar V_0 = NAN;
        Set::Scalar gradient_alpha_refinement_threshold = NAN;
        Set::Scalar eta_refinement_threshold = NAN;
    } agglom;
};
}
#endif
