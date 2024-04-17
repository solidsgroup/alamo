//
// This is a general purpose integrator that focuses on solving 
// elasticity/mechanics equations in the absence of other multiphysics
// simulations.
// It is enabled by :code:`alamo.program=mechanics` if used on its own, in
// which case there is no prefix.
// If it is being used by another integrator, see that integrator to determine
// the value of :code:`[prefix]` (often equal to :code:`elastic`).
//
// This integrator inherits from :ref:`Integrator::Base::Mechanics`; see 
// documentation for that integrator for additional parameters.
//
// :bdg-primary-line:`Model setup`
// There are two basic tools for setting up a mechanics problem in this
// integrator.
//
// 1. The :code:`eta` field: this is used to mix models of different types.
//    Use :code:`nmodels` to specify how many material models to use, and then
//    specify each model as :code:`model1`, :code:`model2`, etc.
//    The type of moel is specified using the :code:`alamo.program.mechanics.model`
//    input.
// 
//    Once you have done this, you must determine the spatial distribution of each
//    model. This is done by choosing an IC for the eta field with :code:`ic.type`.
//    The :code:`IC::Expression` is the most general and recommended.
//    The models are then mixed linearly, i.e.
//
//    .. math::
//
//       W_{\textrm{eff}} = \sum_{i=1}^N W_i\,\eta_i(\mathbf{x})
// 
//    See the :ref:`Eshelby` test for an example of model mixing.
// 
// 2. The :code:`psi` field: this is used specifically for cases where a "void" region
//    is desired. Its usage is similar to the :code:`eta` case, and is conceptually
//    similar in that it scales the model field to near-zero in order to mimic the
//    (lack of) mechanical behavior in an empty region.
//    It is important to use :code:`psi` here, for reasons that are discussed in detail
//    in
//    `this paper <https://doi.org/10.1007/s00466-023-02325-8>`_.
//    The initialization of :code:`psi` is similar to that for :code:`eta`.
//    
//    See the :ref:`PlateHole` and :ref:`RubberPlateHole` for canonical exmaples.
//    The :ref:`Integrator::Fracture` and :ref:`Integrator::TopOp` integrators are examples
//    of integrators that leverage the psi property.
//
// :bdg-primary-line:`Body forces` 
// currently have limited support due to the relatively low number of 
// times they are needed. See the :ref:`Integrator::Base::Mechanics` documentation for 
// detail. 
// See the :ref:`TrigTest` test for examples of current body force implementation.
//
// :bdg-primary-line:`Boundary conditions` 
// are implemented using the :ref:`BC::Operator::Elastic` classes.
// See the documentation on these classes for more detail.
// See any of the mechanics-based tests for examples of boundary condition application.
//

#ifndef INTEGRATOR_MECHANICS_H
#define INTEGRATOR_MECHANICS_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Base/Mechanics.H"


#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "IC/PSRead.H"
#include "IC/Constant.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"


namespace Integrator
{
template<class MODEL>
class Mechanics : virtual public Base::Mechanics<MODEL>
{
public:

    Mechanics() : Base::Mechanics<MODEL>() {
        Util::Message(INFO);
    }
    Mechanics(IO::ParmParse& pp) : Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    // Mechanics inputs. See also :ref:`Integrator::Base::Mechanics`
    static void Parse(Mechanics& value, IO::ParmParse& pp)
    {
        Base::Mechanics<MODEL>::Parse(value, pp);
        int nmodels = 1;
        pp.query("nmodels", nmodels); // Number of elastic model varieties
        for (int i = 0; i < nmodels; i++)
        {
            std::string name = "model" + std::to_string(i + 1);
            MODEL tmp_model;
            pp.queryclass(std::string(name.data()), tmp_model);
            value.models.push_back(tmp_model);
        }
        Util::Assert(INFO, TEST(value.models.size() > 0));
        value.RegisterNodalFab(value.eta_mf, value.models.size(), 2, "eta", true);
        // Refinement threshold for eta field
        pp.query("eta_ref_threshold", value.m_eta_ref_threshold);
        // Refinement threshold for strain gradient
        pp.query("ref_threshold", value.m_elastic_ref_threshold);

        pp.query("fields_verbose", value.plot_field);

        // Read in IC for eta
        if (nmodels>1 && pp.contains("ic.type"))
        {
            std::string type;
            pp.query("ic.type", type); // Read IC type for the eta field
            if (type == "ellipse")    value.ic_eta = new IC::Ellipse(value.geom, pp, "ic.ellipse");
            else if (type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "ic.constant");
            else if (type == "voronoi")    value.ic_eta = new IC::Voronoi(value.geom, pp, "ic.voronoi");
            else if (type == "bmp")        value.ic_eta = new IC::BMP(value.geom, pp, "ic.bmp");
            else if (type == "png")      value.ic_eta = new IC::PNG(value.geom, pp, "ic.png");
            else if (type == "expression") value.ic_eta = new IC::Expression(value.geom, pp, "ic.expression");
            else if (type == "psread") value.ic_eta = new IC::PSRead(value.geom, pp, "ic.psread");
            else Util::Abort(INFO,"Invalid value for ic.type: ",type);
            
            value.eta_reset_on_regrid = true;
            // Whether to re-initialize eta when re-gridding occurs.
            // Default is false unless eta ic is set, then default is.
            // true.
            pp.query("eta.reset_on_regrid",value.eta_reset_on_regrid);
        }

        // Read in IC for psi
        if (pp.contains("psi.ic.type"))
        {
            std::string type;
            pp.query("psi.ic.type", type); // Read IC type for the eta field
            if (type == "ellipse")         value.ic_psi = new IC::Ellipse(value.geom, pp, "psi.ic.ellipse");
            else if (type == "constant")   value.ic_psi = new IC::Constant(value.geom, pp, "psi.ic.constant");
            else if (type == "expression") value.ic_psi = new IC::Expression(value.geom, pp, "psi.ic.expression");
            else if (type == "psread")     value.ic_eta = new IC::PSRead(value.geom, pp, "psi.ic.psread");
            else Util::Abort(INFO, "Invalid value for psi.ic.type: ", type);

            value.bc_psi = new BC::Nothing();
            value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", value.plot_field);
            value.psi_on = true;
            
            value.psi_reset_on_regrid = true;
            // Whether to re-initialize psi when re-gridding occurs.
            // Default is false unless a psi ic is set, then default is
            // true.
            pp.query("psi.reset_on_regrid",value.psi_reset_on_regrid);
        }
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        eta_mf[lev]->setVal(0.0);
        if (models.size() > 1 && ic_eta) ic_eta->Initialize(lev, eta_mf);
        else eta_mf[lev]->setVal(1.0);

        if (psi_on) ic_psi->Initialize(lev, psi_mf);
    }

    virtual void UpdateModel(int a_step) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (a_step > 0) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            //rhs_mf[lev]->setVal(0.0);

            eta_mf[lev]->FillBoundary();

            //disp_mf[lev]->setVal(Set::Vector::Zero());

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();

                amrex::Array4<MODEL> const& model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& eta = eta_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    model(i, j, k) = MODEL::Zero();
                    for (unsigned int n = 0; n < models.size(); n++)
                        model(i, j, k) += eta(i, j, k, n) * models[n];
                });
            }

            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const& eta = eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data(), sten);
                    if (grad.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
            if (psi_on)
            {
                amrex::Array4<Set::Scalar> const& psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    {
                        Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(), sten);
                        if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }
        }
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        if (eta_reset_on_regrid && models.size() > 1 && ic_eta) ic_eta->Initialize(lev, eta_mf, time);
        if (psi_reset_on_regrid) ic_psi->Initialize(lev, psi_mf, time);
    }

protected:
    Set::Field<Set::Scalar> eta_mf;
    Set::Scalar m_eta_ref_threshold = 0.01;
    std::vector<MODEL> models;
    IC::IC* ic_eta = nullptr;
    IC::IC* ic_psi = nullptr;
    BC::BC<Set::Scalar>* bc_psi = nullptr;
    bool psi_reset_on_regrid = false;
    bool eta_reset_on_regrid = false;
    bool plot_field = false;

    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;
};










} // namespace Integrator
#endif
