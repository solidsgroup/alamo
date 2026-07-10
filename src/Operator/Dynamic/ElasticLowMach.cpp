#include "Operator/Dynamic/ElasticLowMach.H"

namespace Operator
{
namespace Dynamic
{
ElasticLowMach::ElasticLowMach(const amrex::Vector<amrex::Geometry>& a_geom,
                               const amrex::Vector<amrex::BoxArray>& a_grids,
                               const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                               const amrex::LPInfo& a_info,
                               const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
ElasticLowMach::define(const amrex::Vector<amrex::Geometry>& a_geom,
                       const amrex::Vector<amrex::BoxArray>& a_grids,
                       const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                       const amrex::LPInfo& a_info,
                       const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
    amrex::MLTensorOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
    DefineCoefficientStorage();
}

void
ElasticLowMach::DefineCoefficientStorage()
{
    shear_face_mf.resize(NAMRLevels());
    bulk_face_mf.resize(NAMRLevels());
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_ba = m_grids[amrlev][0];
            face_ba.surroundingNodes(d);
            shear_face_mf[amrlev][d].define(face_ba, m_dmap[amrlev][0], 1, 0);
            bulk_face_mf[amrlev][d].define(face_ba, m_dmap[amrlev][0], 1, 0);
            shear_face_mf[amrlev][d].setVal(0.0);
            bulk_face_mf[amrlev][d].setVal(0.0);
        }
    }
}

void
ElasticLowMach::SetCoefficients(int amrlev,
                                const amrex::MultiFab& eta,
                                const amrex::MultiFab& rho,
                                const Model::Solid::Finite::NeoHookean& model,
                                Set::Scalar eta_threshold,
                                Set::Scalar density_floor,
                                Set::Scalar coeff_scale)
{
    BL_PROFILE("Operator::Dynamic::ElasticLowMach::SetCoefficients");

    amrex::MultiFab shear_cc(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
    amrex::MultiFab bulk_cc (m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
    shear_cc.setVal(0.0);
    bulk_cc.setVal(0.0);

    eta_threshold = std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), eta_threshold));
    const Set::Scalar denom = std::max(Set::Scalar(1.0e-12), Set::Scalar(1.0) - eta_threshold);
    const Set::Scalar mu = coeff_scale * model.mu;
    const Set::Scalar kappa = coeff_scale * model.kappa;

    for (amrex::MFIter mfi(shear_cc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const Set::Scalar> const& eta_arr = eta.const_array(mfi);
        amrex::Array4<const Set::Scalar> const& rho_arr = rho.const_array(mfi);
        amrex::Array4<Set::Scalar> const& shear = shear_cc.array(mfi);
        amrex::Array4<Set::Scalar> const& bulk = bulk_cc.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta_val = eta_arr(i,j,k);
            eta_val = eta_val < 0.0 ? 0.0 : (eta_val > 1.0 ? 1.0 : eta_val);
            Set::Scalar solid_weight = 0.0;
            if (eta_val > eta_threshold)
            {
                Set::Scalar s = (eta_val - eta_threshold) / denom;
                s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
                solid_weight = s * s * s * (s * (s * 6.0 - 15.0) + 10.0);
            }
            Set::Scalar inv_rho = 1.0 / std::max(rho_arr(i,j,k), density_floor);
            shear(i,j,k) = solid_weight * mu * inv_rho;
            bulk(i,j,k) = solid_weight * kappa * inv_rho;
        });
    }
    shear_cc.FillBoundary(Geom(amrlev).periodicity());
    bulk_cc.FillBoundary(Geom(amrlev).periodicity());

    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> shear_ptr;
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> bulk_ptr;
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> shear_const_ptr;
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> bulk_const_ptr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        shear_ptr[d] = &shear_face_mf[amrlev][d];
        bulk_ptr[d] = &bulk_face_mf[amrlev][d];
        shear_const_ptr[d] = &shear_face_mf[amrlev][d];
        bulk_const_ptr[d] = &bulk_face_mf[amrlev][d];
    }

    amrex::average_cellcenter_to_face(shear_ptr, shear_cc, Geom(amrlev), 1, true, 0);
    amrex::average_cellcenter_to_face(bulk_ptr, bulk_cc, Geom(amrlev), 1, true, 0);

    setShearViscosity(amrlev, shear_const_ptr);
    setBulkViscosity(amrlev, bulk_const_ptr);
}
}
}
