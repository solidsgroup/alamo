#include "ChemReduction.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include <cmath>

namespace Integrator
{

ChemReduction::ChemReduction(IO::ParmParse& pp) : ChemReduction()
{
    pp.queryclass(*this);
}

void
ChemReduction::Parse(ChemReduction& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::ChemReduction::ChemReduction()");
    {
        value.bc_phi = new BC::Constant(1);
        pp.queryclass("pf.phi.bc", *static_cast<BC::Constant*>(value.bc_phi)); // See :ref:`BC::Constant`
        value.bc_chi = new BC::Constant(1);
        pp.queryclass("pf.chi.bc", *static_cast<BC::Constant*>(value.bc_chi)); // See :ref:`BC::Constant`

        { // Parse Properties
            pp.query_required("molar_volume", value.V); // Domain molar volume 
            pp.query_required("temperature", value.T); // Domain temperature (Isothermal)
            pp.query_default("gas_constant", value.R, 8.3145); // Universal gas constant

            pp.query_default("refinement_criterion", value.phi_refinement_criterion, 0.001); // Criterion for mesh refinement
            pp.query_default("ghost_cells_count", value.ghost_count, 2); // Number of ghost cells in each field

            pp.query_required("wustite.bulk", value.wustite.K); // Wustite bulk modulus 
            pp.query_required("wustite.shear", value.wustite.G); // Wustite shear modulus
            pp.query_required("wustite.mobility", value.wustite.L); // Wustite mobility coeficient
            pp.query_required("wustite.diffusivity.oxygen", value.wustite.D.O); // Wustite oxygen diffusion coeficient
            pp.query_required("wustite.diffusivity.hydrogen", value.wustite.D.H); // Wustite hydrogen diffusion coeficient
            pp.query_required("wustite.diffusivity.water", value.wustite.D.W); // Wustite water diffusion coeficient 
            pp.query_required("wustite.energy.gradient", value.wustite.E.G); // Wustite energy gradient
            pp.query_required("wustite.energy.barrier", value.wustite.E.B); // Wustite energy barrier

            pp.query_required("ferrite.bulk", value.ferrite.K); // Ferrite bulk modulus 
            pp.query_required("ferrite.shear", value.ferrite.G); // Ferrite shear modulus 
            pp.query_required("ferrite.mobility", value.ferrite.L); // Ferrite mobility coeficient
            pp.query_required("ferrite.diffusivity.oxygen", value.ferrite.D.O); // Ferrite oxygen diffusion coeficient
            pp.query_required("ferrite.diffusivity.hydrogen", value.ferrite.D.H);// Ferrite hydrogen diffusion coeficient
            pp.query_required("ferrite.diffusivity.water", value.ferrite.D.W);// Ferrite water diffusion coeficient
            pp.query_required("ferrite.energy.gradient", value.ferrite.E.G);// Ferrite energy gradient
            pp.query_required("ferrite.energy.barrier", value.ferrite.E.B);// Ferrite energy barrier

            pp.query_required("gas.bulk", value.gas.K); // Water/gas bulk modulus
            pp.query_required("gas.shear", value.gas.G); // Water/gas shear modulus
            pp.query_required("gas.mobility", value.gas.L); // Water/gas mobility coeficient 
            pp.query_required("gas.diffusivity.oxygen", value.gas.D.O); // Water/gas oxygen diffusion coeficient
            pp.query_required("gas.diffusivity.hydrogen", value.gas.D.H); // Water/gas hydrogen diffusion coeficient
            pp.query_required("gas.diffusivity.water", value.gas.D.W); // Water/gas water diffusion coeficient 
            pp.query_required("gas.energy.gradient", value.gas.E.G); // Water/gas energy gradient
            pp.query_required("gas.energy.barrier", value.gas.E.B); // Water/gas energy barrier

        } // End Parse Properties

        { // Allen-Cahn (Non-Conserved) fields
            value.RegisterNewFab(value.phi_g_mf, value.bc_phi, 1, value.ghost_count + 1, "phiGas", true);
            value.RegisterNewFab(value.phi_f_mf, value.bc_phi, 1, value.ghost_count + 1, "phiFerrite", true);
            value.RegisterNewFab(value.phi_w_mf, value.bc_phi, 1, value.ghost_count + 1, "phiWustite", true);

            std::string phi_bc_str = "constant";
            pp.query("pf.phi.ic.type", phi_bc_str); // phi boundary condition [constant, expression]
            if (phi_bc_str == "constant") value.ic_phi = new IC::Constant(value.geom, pp, "pf.phi.ic.constant");
            else if (phi_bc_str == "expression") value.ic_phi = new IC::Expression(value.geom, pp, "pf.phi.ic.expression");

            std::string phi_ic_type = "constant";
            pp.query("phi.ic.type", phi_ic_type); // phi initial condition [constant, laminate, expression, bmp]
            if (phi_ic_type == "constant") value.ic_phi = new IC::Constant(value.geom, pp, "phi.ic.constant");
            else if (phi_ic_type == "expression") value.ic_phi = new IC::Expression(value.geom, pp, "phi.ic.expression");
            else if (phi_ic_type == "bmp") value.ic_phi = new IC::BMP(value.geom, pp, "phi.ic.bmp");
            else if (phi_ic_type == "png") value.ic_phi = new IC::PNG(value.geom, pp, "phi.ic.png");
            else Util::Abort(INFO, "Invalid phi IC type", phi_ic_type);
        } // End Allen-Cahn Register Fab

        { // Cahn-Hillard (Conserved) fields
            value.RegisterNewFab(value.chi_o_mf, value.bc_chi, 1, value.ghost_count + 1, "chiOxygen", true);
            value.RegisterNewFab(value.chi_w_mf, value.bc_chi, 1, value.ghost_count + 1, "chiWater", true);
            value.RegisterNewFab(value.chi_h_mf, value.bc_chi, 1, value.ghost_count + 1, "chiHydrogen", true);

            std::string chi_bc_str = "constant";
            pp.query("pf.chi.ic.type", chi_bc_str); // chi boundary condition [constant, expression]
            if (chi_bc_str == "constant") value.ic_chi = new IC::Constant(value.geom, pp, "pf.chi.ic.constant");
            else if (chi_bc_str == "expression") value.ic_chi = new IC::Expression(value.geom, pp, "pf.chi.ic.expression");

            std::string chi_ic_type = "constant";
            pp.query("chi.ic.type", chi_ic_type); // chi initial condition [constant, laminate, expression, bmp]
            if (chi_ic_type == "constant") value.ic_chi = new IC::Constant(value.geom, pp, "chi.ic.constant");
            else if (chi_ic_type == "expression") value.ic_chi = new IC::Expression(value.geom, pp, "chi.ic.expression");
            else if (chi_ic_type == "bmp") value.ic_chi = new IC::BMP(value.geom, pp, "chi.ic.bmp");
            else if (chi_ic_type == "png") value.ic_chi = new IC::PNG(value.geom, pp, "chi.ic.png");
            else Util::Abort(INFO, "Invalid phi IC type", chi_ic_type);
        } // End Cahn-Hillard Register Fab       
    } // End BL_PROFILE
} // End Parse Function 

void ChemReduction::Initialize(int lev)
{
    BL_PROFILE("Integrator::ChemReduction::Initialize");
    {
        ferrite.M.O = (ferrite.D.O * V) / (R * T);
        ferrite.M.H = (ferrite.D.H * V) / (R * T);
        ferrite.M.W = (ferrite.D.W * V) / (R * T);

        wustite.M.O = (wustite.D.O * V) / (R * T);
        wustite.M.H = (wustite.D.H * V) / (R * T);
        wustite.M.W = (wustite.D.W * V) / (R * T);

        gas.M.O = (gas.D.O * V) / (R * T);
        gas.M.H = (gas.D.H * V) / (R * T);
        gas.M.W = (gas.D.W * V) / (R * T);

        ic_phi->Initialize(lev, phi_g_mf, phi_f_mf, phi_w_mf);
        //ic_phi->Initialize(lev, phi_g_old_mf, phi_f_old_mf, phi_w_old_mf);

        ic_chi->Initialize(lev, chi_h_mf, chi_o_mf, chi_w_mf);
        //ic_chi->Initialize(lev, chi_h_old_mf, chi_o_old_mf, chi_fe_old_mf);
    } // End BL_PROFILE
} // End Initialize Function

void ChemReduction::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::ChemReduction::Advance");
    {
        const Set::Scalar* DX = geom[lev].CellSize();

        // Free Energy Functions
        Numeric::Function::Polynomial<4> f_f(-40372, -36.031, 90.573, -182.51, 163.708);
        Numeric::Function::Polynomial<2> f_w(41.953, -213.46, 173.263);

        // Interpolation Functions
        Numeric::Function::Polynomial<5> p(0., 0., 0., 10., -15., 6.);
        Numeric::Function::Polynomial<4> g(0., 0., 1., -2., 1.);

        for (amrex::MFIter mfi(*phi_f_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();

            // Parallel Loop for Cahn-Hillard Fields
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                p(1.);
            });

            // Parallel Loop for Allen-Cahn Fields
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                g(1.);
            });
        }// End For Loop
    } // End BL_PROFILE
} //Function

void ChemReduction::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::ChemReduction::TagCellsForRefinement");
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Phi criterion for refinement
        for (amrex::MFIter mfi(*phi_w_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const& phiF = (*phi_f_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& phiW = (*phi_w_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradf = Numeric::Gradient(phiF, i, j, k, 0, DX);
                Set::Vector gradw = Numeric::Gradient(phiW, i, j, k, 0, DX);
                if (gradf.lpNorm<2>() * dr * 2 > phi_refinement_criterion || gradw.lpNorm<2>() * dr * 2 > phi_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    } // End BL_PROFILE
} // End TagCell

void ChemReduction::Regrid(int lev, Set::Scalar time)
{
    BL_PROFILE("Integrator::ChemReduction::Regrid");
    {
        if (lev < finest_level) return;
        Util::Message(INFO, "Regridding on level", lev);
    } // End BL_PROFILE
} // End Regrid 

} // namespace Integrator


