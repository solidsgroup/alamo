#ifndef NUMERIC_REFERENCEMAP_RECONSTRUCTION_H
#define NUMERIC_REFERENCEMAP_RECONSTRUCTION_H

#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "BC/BC.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"

namespace Numeric
{
namespace ReferenceMap
{
/// \brief Extrapolation and smoothing used to condition a reference map.
///
/// These operations repair xi where eta does not identify trusted solid data.
/// They are numerical processing and are deliberately kept separate from the
/// reference-map advection and solid constitutive response.
class Reconstruction
{
public:
    static void Parse(Reconstruction& value, IO::ParmParse& pp);

    void operator()(const amrex::Geometry& geom,
                    const amrex::MultiFab& eta,
                    amrex::MultiFab& xi,
                    BC::BC<Set::Scalar>& xi_bc,
                    Set::Scalar time) const;

    void SmoothForStress(const amrex::Geometry& geom,
                         const amrex::MultiFab& eta,
                         amrex::MultiFab& xi) const;

    bool SmoothsStressMap() const
    {
        return stress_smoothing_sweeps > 0 && stress_smoothing_alpha != 0.0;
    }

    Set::Scalar EtaExtension() const { return eta_extension; }

private:
    Set::Scalar eta_core = 0.5;
    Set::Scalar eta_extension = 1.0e-3;
    int extrapolation_sweeps = 4;
    Set::Scalar reconstruction_alpha = 1.0;
    Set::Scalar reconstruction_power = 1.0;
    Set::Scalar affine_tolerance = 1.0e-7;
    int smoothing_sweeps = 2;
    Set::Scalar smoothing_alpha = 1.0;
    Set::Scalar smoothing_power = 2.0;
    int stress_smoothing_sweeps = 0;
    Set::Scalar stress_smoothing_alpha = 0.1;
};
}
}

#endif
