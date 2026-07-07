#ifndef OPERATOR_DYNAMIC_ELASTICLOWMACH_H_
#define OPERATOR_DYNAMIC_ELASTICLOWMACH_H_

#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include "Model/Solid/Finite/NeoHookean.H"
#include "Set/Set.H"

namespace Operator
{
namespace Dynamic
{
class ElasticLowMach : public amrex::MLTensorOp
{
public:
    ElasticLowMach() = default;
    ElasticLowMach(const amrex::Vector<amrex::Geometry>& a_geom,
                   const amrex::Vector<amrex::BoxArray>& a_grids,
                   const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                   const amrex::LPInfo& a_info = amrex::LPInfo(),
                   const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory = {});

    void define(const amrex::Vector<amrex::Geometry>& a_geom,
                const amrex::Vector<amrex::BoxArray>& a_grids,
                const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                const amrex::LPInfo& a_info = amrex::LPInfo(),
                const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory = {});

    void SetCoefficients(int amrlev,
                         const amrex::MultiFab& eta,
                         const amrex::MultiFab& rho,
                         const Model::Solid::Finite::NeoHookean& model,
                         Set::Scalar eta_threshold,
                         Set::Scalar density_floor,
                         Set::Scalar coeff_scale);

private:
    void DefineCoefficientStorage();

    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> shear_face_mf;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> bulk_face_mf;
};
}
}

#endif
