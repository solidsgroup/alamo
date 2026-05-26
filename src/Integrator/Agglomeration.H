#ifndef INTEGRATOR_AGGLOMERATION_H
#define INTEGRATOR_AGGLOMERATION_H

#include "BC/BC.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Integrator/Flame.H"
#include "Set/Base.H"

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
    void ScaleByComplement(MultiFab &dst, const MultiFab &src, int srccomp, int dstcomp, int numcomp, int nghost);
    void ScaleByPhiComplement(const Set::Field<Set::Scalar> &mf, int lev);
    void Initialize(int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &box) override;

    struct
    {
        Set::Field<Set::Scalar> alpha_old;
        Set::Field<Set::Scalar> alpha;
        Set::Field<Set::Scalar> dalpha;
        Set::Field<Set::Scalar> driving_force;
        Set::Field<Set::Vector> evolution_driving_force;
        IC::IC<Set::Scalar> *alpha_ic = nullptr;
        BC::BC<Set::Scalar> *alpha_bc = nullptr;
        Set::Scalar p = NAN;
        Set::Scalar epsilon = NAN;
        Set::Scalar omega = NAN;
        Set::Scalar gamma = NAN;
        Set::Scalar kappa_g = NAN;
        Set::Scalar mu = NAN;
        Set::Scalar kappa_M = NAN;
        Set::Scalar q = NAN;
        Set::Scalar V = 0.0;
    } agglom;
};
}
#endif
