#ifndef SUTURECRACK_INTEGRATOR_H
#define SUTURECRACK_INTEGRATOR_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "BC/Constant.H"

#include "IC/IC.H"
#include "IC/Ellipsoid.H"
#include "IC/Notch.H"
#include "IC/BMP.H"

#include "Operator/Operator.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"
#include "Model/Solid/Affine/J2PlasticDegradable.H"
#include "Model/Solid/Affine/CrystalPlasticDegradable.H"

#include "Model/Interface/Crack/Crack.H"
#include "Model/Interface/Crack/Constant.H"
#include "Model/Interface/Crack/Sin.H"

#include "Numeric/Stencil.H"
#include <eigen3/Eigen/Dense>

namespace Integrator
{
using brittle_fracture_model_type = Model::Solid::Linear::IsotropicDegradable;

class SutureCrack : public Integrator
{

public:
    SutureCrack()
    {
        IO::ParmParse pp_crack("crack");
        pp_crack.query("modulus_scaling_max",crack.scaleModulusMax);
        pp_crack.query("refinement_threshold",crack.refinement_threshold);

        Model::Interface::Crack::Constant *tmpbdy = new Model::Interface::Crack::Constant();
        pp_crack.queryclass("constant1",*tmpbdy);
        crack.cracktype1 = tmpbdy;

        Model::Interface::Crack::Constant *tmpbdy2 = new Model::Interface::Crack::Constant();
        pp_crack.queryclass("constant2",*tmpbdy2);
        crack.cracktype2 = tmpbdy2;

        IO::ParmParse pp_crack_df("crack.df");
        pp_crack_df.query("mult_Gc", crack.mult_df_Gc);
        pp_crack_df.query("mult_Lap", crack.mult_df_lap);

        IO::ParmParse pp("crack.ic");
        IC::Notch *tmpic = new IC::Notch(geom);
        pp.queryclass("notch",*tmpic);
        crack.ic = tmpic;

        crack.bc = new BC::Constant(1);
        crack.bcdf = new BC::Constant(4);
        pp_crack.queryclass("bc",*static_cast<BC::Constant *>(crack.bc));
        pp_crack.queryclass("bc_df",*static_cast<BC::Constant *>(crack.bcdf));

        // RegisterNewFab(crack.field,     crack.bc, 1, number_of_ghost_cells, "c",		true);
        // RegisterNewFab(crack.field_old, crack.bc, 1, number_of_ghost_cells, "c_old",	true);
        // RegisterNewFab(crack.driving_force, crack.bcdf, 4, number_of_ghost_cells, "driving_force",true);

        RegisterNodalFab(crack.field,       1, number_of_ghost_nodes, "c",		true);
        RegisterNodalFab(crack.field_old,   1, number_of_ghost_nodes, "c_old",	true);
        RegisterNodalFab(crack.driving_force, 4, number_of_ghost_nodes, "driving_force",true);

        IO::ParmParse pp_material1("material1");
        pp_material1.queryclass("isotropic",material.brittlemodeltype1);
        IO::ParmParse pp_material2("material2");
        pp_material2.queryclass("isotropic",material.brittlemodeltype2);

        IO::ParmParse pp_material_ic("material.ic");
        std::string type;
        pp_material_ic.query("type",type);
        if (type=="ellipse")
        {
            material.ic = new IC::Ellipse(geom);
            pp_material_ic.queryclass("ellipse",static_cast<IC::Ellipse*>(material.ic));
        }
        else if (type=="bmp")
        {
            material.ic = new IC::BMP(geom);
            pp_material_ic.queryclass("bmp",static_cast<IC::BMP*>(material.ic));
        }

        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);
        pp_load.query("val", loading.val);

        IO::ParmParse pp_elastic("elastic");
        pp_elastic.query("df_mult", elastic.df_mult);
        pp_elastic.queryclass("bc",elastic.brittlebc);

        IO::ParmParse pp_solver("solver");
        pp_solver.query("int",				sol.interval);
        pp_solver.query("type",			sol.type);
        pp_solver.query("max_iter",		sol.max_iter);
        pp_solver.query("max_fmg_iter",	sol.max_fmg_iter);
        pp_solver.query("verbose",			sol.verbose);
        pp_solver.query("cgverbose",		sol.cgverbose);
        pp_solver.query("tol_rel",			sol.tol_rel);
        pp_solver.query("tol_abs",			sol.tol_abs);
        pp_solver.query("cg_tol_rel",		sol.cg_tol_rel);
        pp_solver.query("cg_tol_abs",		sol.cg_tol_abs);
        pp_solver.query("use_fsmooth",		sol.use_fsmooth);
        pp_solver.query("agglomeration", 	sol.agglomeration);
        pp_solver.query("consolidation", 	sol.consolidation);
        pp_solver.query("bottom_solver",       sol.bottom_solver);
        pp_solver.query("linop_maxorder",      sol.linop_maxorder);
        pp_solver.query("max_coarsening_level",sol.max_coarsening_level);
        pp_solver.query("verbose",             sol.verbose);
        pp_solver.query("cg_verbose",          sol.cgverbose);
        pp_solver.query("bottom_max_iter",     sol.bottom_max_iter);
        pp_solver.query("max_fixed_iter",      sol.max_fixed_iter);
        pp_solver.query("bottom_tol",          sol.bottom_tol);

        nlevels = maxLevel() + 1;
        RegisterNodalFab(elastic.disp,  AMREX_SPACEDIM, number_of_ghost_nodes, "disp", true);
        RegisterNodalFab(elastic.rhs,  AMREX_SPACEDIM, number_of_ghost_nodes, "rhs", true);
        RegisterNodalFab(elastic.residual,  AMREX_SPACEDIM, number_of_ghost_nodes, "res", true);
        RegisterNodalFab(elastic.strain,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strain", true);
        RegisterNodalFab(elastic.stress,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "stress", true);
        RegisterNodalFab(elastic.energy, 1, number_of_ghost_nodes, "energy", true);
        RegisterNodalFab(elastic.energy_pristine, 1, number_of_ghost_nodes, "energy_pristine", true);
        RegisterNodalFab(elastic.energy_pristine_old, 1, number_of_ghost_nodes, "energy_pristine_old", true);
        // RegisterNewFab(material.modulus_field, crack.bc, 1, number_of_ghost_cells, "modulus_field", true);
        RegisterNodalFab(material.modulus_field, 1, number_of_ghost_nodes, "modulus_field", true);
        RegisterGeneralFab(material.brittlemodel, 1, number_of_ghost_nodes);

    }

protected:
    void Initialize(int ilev) override
    {
        Util::Message(INFO, "ilev = ", ilev);
        crack.ic->Initialize(ilev,crack.field);
        crack.ic->Initialize(ilev,crack.field_old);
        
        Util::Message(INFO, "ilev = ", ilev);
        elastic.disp[ilev]->setVal(0.0);
        elastic.strain[ilev]->setVal(0.0);
        elastic.stress[ilev]->setVal(0.0);
        elastic.rhs[ilev]->setVal(0.0);
        elastic.energy[ilev]->setVal(0.0);
        elastic.residual[ilev]->setVal(0.0);
        elastic.energy_pristine[ilev] -> setVal(0.);
        elastic.energy_pristine_old[ilev] -> setVal(0.);

        Util::Message(INFO, "ilev = ", ilev);
        material.ic->Initialize(ilev,material.modulus_field);
        //material.brittlemodel[ilev]->setVal(material.brittlemodeltype1);
    }

    void TimeStepBegin(Set::Scalar /*time*/, int iter) override
    {
        if(iter%sol.interval) return;
        material.brittlemodeltype1.DegradeModulus(0.0);
        material.brittlemodeltype2.DegradeModulus(0.0);
        for (int ilev = 0; ilev < nlevels; ++ilev) std::swap(elastic.energy_pristine_old[ilev], elastic.energy_pristine[ilev]);

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            material.modulus_field[ilev]->FillBoundary();
            for (amrex::MFIter mfi(*elastic.disp[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.growntilebox(2);
                amrex::Array4<brittle_fracture_model_type> modelfab = (material.brittlemodel)[ilev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[ilev]).array(mfi);
                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    //Set::Scalar _temp = Numeric::Interpolate::CellToNodeAverage(modbox,i,j,k,0);
                    Set::Scalar _temp = modbox(i,j,k,0);
                    modelfab(i,j,k) = (1-_temp)*material.brittlemodeltype1 + (_temp)*material.brittlemodeltype2;
                });
            }
            Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }

        // Degrading the modulus
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            crack.field[ilev]->FillBoundary();
            material.modulus_field[ilev]->FillBoundary();
            for (amrex::MFIter mfi(*elastic.disp[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.validbox();//growntilebox(2);
                amrex::Array4<const Set::Scalar> const& c_new = (*crack.field[ilev]).array(mfi);
                amrex::Array4<brittle_fracture_model_type> modelfab = (material.brittlemodel)[ilev]->array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    // Set::Scalar _temp = 0;
                    Set::Scalar _temp = crack.cracktype1->g_phi(c_new(i,j,k,0),0.);
                    // Set::Scalar mul = AMREX_D_PICK(0.5,0.25,0.125);
                    // if(std::isinf(c_new(i,j,k,0))) Util::Abort(INFO, "infs in c at (",i,",",j,",",k,")");


                    // _temp =  mul*(AMREX_D_TERM(	
                    //             crack.cracktype1->g_phi(c_new(i,j,k,0),0.) + crack.cracktype1->g_phi(c_new(i-1,j,k,0),0.)
                    //             , 
                    //             + crack.cracktype1->g_phi(c_new(i,j-1,k,0),0.) + crack.cracktype1->g_phi(c_new(i-1,j-1,k,0),0.)
                    //             , 
                    //             + crack.cracktype1->g_phi(c_new(i,j,k-1,0),0.) + crack.cracktype1->g_phi(c_new(i-1,j,k-1,0),0.)
                    //             + crack.cracktype1->g_phi(c_new(i,j-1,k-1,0),0.) + crack.cracktype1->g_phi(c_new(i-1,j-1,k-1,0),0.))
                    //             );
                    if (std::isnan(_temp)) Util::Abort(INFO);
                    if(_temp < 0.0) _temp = 0.;
                    if(_temp > 1.0) _temp = 1.0;
                    modelfab(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-crack.scaleModulusMax));
                });
            }
            Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }

        // Body force
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            const Real* DX = geom[ilev].CellSize();
            Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

            AMREX_D_TERM(elastic.rhs[ilev]->setVal(loading.body_force(0)*volume,0,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(1)*volume,1,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(2)*volume,2,1););
        }

        // Mode I - displacement
        // elastic.brittlebc.Set(elastic.brittlebc.Face::YHI, elastic.brittlebc.Direction::Y, BC::Operator::Elastic<brittle_fracture_model_type>::Type::Displacement, loading.val);
        // elastic.brittlebc.Set(elastic.brittlebc.Face::XLO_YHI, elastic.brittlebc.Direction::Y, BC::Operator::Elastic<brittle_fracture_model_type>::Type::Displacement, loading.val);
        // elastic.brittlebc.Set(elastic.brittlebc.Face::XHI_YHI, elastic.brittlebc.Direction::Y, BC::Operator::Elastic<brittle_fracture_model_type>::Type::Displacement, loading.val);
        elastic.brittlebc.Init(elastic.rhs,geom);

        Operator::Elastic<brittle_fracture_model_type::sym> op_b;
        LPInfo info;
        info.setAgglomeration(sol.agglomeration);
        info.setConsolidation(sol.consolidation);
        info.setMaxCoarseningLevel(sol.max_coarsening_level);

        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp[ilev]->contains_nan()) Util::Warning(INFO);

        op_b.define(geom, grids, dmap, info);
        op_b.setMaxOrder(sol.linop_maxorder);
        op_b.SetUniform(false);
        op_b.SetBC(&elastic.brittlebc);

        Solver::Nonlocal::Newton<brittle_fracture_model_type>  solver(op_b);
        solver.setMaxIter(sol.max_iter);
        solver.setMaxFmgIter(sol.max_fmg_iter);
        solver.setFixedIter(sol.max_fixed_iter);
        solver.setVerbose(sol.verbose);
        solver.setCGVerbose(sol.cgverbose);
        solver.setBottomMaxIter(sol.bottom_max_iter);
        //solver.setBottomTolerance(sol.cg_tol_rel) ;
        //solver.setBottomToleranceAbs(sol.cg_tol_abs) ;
        if (sol.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
        else if (sol.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
        solver.solve(elastic.disp, elastic.rhs, material.brittlemodel, sol.tol_rel, sol.tol_abs);
        solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,material.brittlemodel);

        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            op_b.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
            op_b.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
            op_b.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);

            elastic.strain[ilev]->FillBoundary();
            elastic.energy_pristine[ilev]->setVal(0.0);
            elastic.energy_pristine_old[ilev]->FillBoundary();

            for (amrex::MFIter mfi(*elastic.strain[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.growntilebox(2);
                amrex::Array4<const Set::Scalar> const& strain_box 	= (*elastic.strain[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box 		= (*elastic.energy_pristine[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box_old 	= (*elastic.energy_pristine_old[ilev]).array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Matrix eps = Numeric::FieldToMatrix(strain_box,i,j,k);
                    Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(eps);
                    Set::Vector eValues = eigensolver.eigenvalues();
                    Set::Matrix eVectors = eigensolver.eigenvectors();

                    Set::Matrix epsp = Set::Matrix::Zero();
                    Set::Matrix epsn = Set::Matrix::Zero();

                    for (int n = 0; n < AMREX_SPACEDIM; n++)
                    {
                        if(eValues(n) > 0.0) epsp += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                        else epsn += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    }

                    Set::Scalar _temp = Numeric::Interpolate::CellToNodeAverage(modbox,i,j,k,0);
                    brittle_fracture_model_type brittlemodeltype = _temp*material.brittlemodeltype1 + (1.-_temp)*material.brittlemodeltype2;
                    
                    energy_box(i,j,k,0) = brittlemodeltype.W(epsp);
                    energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
                });
            }
        }
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        std::swap(crack.field_old[lev], crack.field[lev]);

        crack.field_old[lev]->FillBoundary();
        elastic.energy_pristine[lev]->FillBoundary();
        
        const Set::Scalar* DX = geom[lev].CellSize();

        for ( amrex::MFIter mfi(*crack.field[lev],true); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& c_old = (*crack.field_old[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& df = (*crack.driving_force[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& c_new = (*crack.field[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& energy_box = (*elastic.energy_pristine[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[lev]).array(mfi);
            
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
    			Set::Scalar rhs = 0.0;
                //Set::Scalar _temp = Numeric::Interpolate::NodeToCellAverage(modbox,i,j,k,0);
                Set::Scalar _temp = modbox(i,j,k,0);

                if(std::isnan(c_old(i,j,k,0))) Util::Abort(INFO, "nans in c at (",i,",",j,",",k,")");
                if(std::isinf(c_old(i,j,k,0))) Util::Abort(INFO, "infs in c at (",i,",",j,",",k,")");

                Set::Scalar Gc = _temp*crack.cracktype1->Gc(0.0) + (1-_temp)*crack.cracktype2->Gc(0.0);
                Set::Scalar Zeta = _temp*crack.cracktype1->Zeta(0.0) + (1-_temp)*crack.cracktype2->Zeta(0.0);
                Set::Scalar DrivingForceThreshold = _temp*crack.cracktype1->DrivingForceThreshold(c_old(i,j,k,0))
                                                    + (1-_temp)*crack.cracktype2->DrivingForceThreshold(c_old(i,j,k,0));
                Set::Scalar Mobility = _temp*crack.cracktype1->Mobility(c_old(i,j,k,0)) + (1-_temp)*crack.cracktype2->Mobility(c_old(i,j,k,0));

                //Set::Scalar _temp = Numeric::Interpolate::NodeToCellAverage(modbox,i,j,k,0);
                Set::Scalar en_cell = Numeric::Interpolate::NodeToCellAverage(energy_box,i,j,k,0);
                if (std::isnan(en_cell)) Util::Abort(INFO, "Nans detected in en_cell. energy_box(i,j,k,0) = ", energy_box(i,j,k,0));

                df(i,j,k,0) = crack.cracktype1->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;
                rhs += crack.cracktype1->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;

                Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
			    Set::Scalar laplacian = DDc.trace();

                df(i,j,k,1) = Gc*crack.cracktype1->Dw_phi(c_old(i,j,k,0),0.0)/(4.0*Zeta)*crack.mult_df_Gc;
				df(i,j,k,2) = 2.0*Zeta*Gc*laplacian*crack.mult_df_lap;

				rhs += Gc*crack.cracktype1->Dw_phi(c_old(i,j,k,0),0.)/(4.0*Zeta)*crack.mult_df_Gc;
				rhs -= 2.0*Zeta*Gc*laplacian*crack.mult_df_lap;

                //rhs *= crack.cracktype->g_phi(_temp,0.0);

                df(i,j,k,3) = std::max(0.,rhs - DrivingForceThreshold);
                if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.cracktype1->Dw_phi(c_old(i,j,k,0),0.0),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
			    c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0., rhs - DrivingForceThreshold)*Mobility;
            });
        }
        crack.field[lev]->FillBoundary();
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        const amrex::Real *DX = geom[lev].CellSize();
        const Set::Vector dx(DX);
        const Set::Scalar dxnorm = dx.lpNorm<2>();

        for (amrex::MFIter mfi(*material.modulus_field[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box 							&bx 	= mfi.tilebox();
            amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const 		&modbox 	= (*material.modulus_field[lev]).array(mfi);
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(modbox, i, j, k, 0, DX);
                if (dxnorm * grad.lpNorm<2>() > crack.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }

        for (amrex::MFIter mfi(*crack.field[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box 							&bx 	= mfi.tilebox();
            amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const 		&c_new 	= (*crack.field[lev]).array(mfi);
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX);
                if (dxnorm * grad.lpNorm<2>() > crack.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

private:

    int number_of_ghost_cells = 3;				///< Number of ghost cells
	int number_of_ghost_nodes = 2;				///< Number of ghost nodes
    int nlevels;

    struct{
        Set::Field<Set::Scalar> disp;             ///< displacement field
        Set::Field<Set::Scalar> strain;           ///< total strain field (gradient of displacement)
        Set::Field<Set::Scalar> stress;           ///< stress field
        Set::Field<Set::Scalar> rhs;              ///< rhs fab for elastic solution
        Set::Field<Set::Scalar> residual;         ///< residual field for solver
        Set::Field<Set::Scalar> energy;           ///< total elastic energy
        Set::Field<Set::Scalar> energy_pristine;      ///< energy of the prisitne material as if no crack is present
        Set::Field<Set::Scalar> energy_pristine_old;  ///< energy of the pristine material for previous time step.

        BC::Operator::Elastic 	brittlebc;  ///< elastic BC if using brittle fracture
		Set::Scalar df_mult = 1.0;              ///< mulitplier for elastic driving force.
    } elastic;

    struct{
        Set::Field<Set::Scalar> field;                ///< crack field at current time step
        Set::Field<Set::Scalar> field_old;            ///< crack field at previous time step
        Set::Field<Set::Scalar> driving_force;         ///< crack driving forces.

        Model::Interface::Crack::Crack *cracktype1;       ///< type of crack. See Crack/Constant or Crack/Sin
        Model::Interface::Crack::Crack *cracktype2;
		std::string ic_type;                            ///< crack IC type. See IC/Notch and IC/Ellipsoid
		IC::IC *ic;                                     ///< crack IC. See IC/Notch and IC/Ellipsoid
        BC::BC<Set::Scalar> *bc;                      ///< BC objects for crack fields
		BC::BC<Set::Scalar> *bcdf;                     ///< BC object for crack driving force

        Set::Scalar mobility;                       ///< crack mobility
        Set::Scalar scaleModulusMax = 0.02;         ///< material modulus ratio inside crack (default = 0.02).
        Set::Scalar refinement_threshold = 0.001;

		Set::Scalar mult_df_Gc = 1.0;				///< Multiplier for Gc/zeta term
		Set::Scalar mult_df_lap = 1.0;				///< Multiplier for the laplacian term
    } crack;

    struct{
        brittle_fracture_model_type brittlemodeltype1;
        brittle_fracture_model_type brittlemodeltype2;
        IC::IC *ic;
        Set::Field<Set::Scalar> modulus_field;
        Set::Field<brittle_fracture_model_type> brittlemodel;
        std::string input_material = "isotropic";
    } material;

    struct{
		Set::Scalar bottom_tol 				= 1.0E-6;
		int         interval	 			= 100;
		std::string type					= "single";
		int         max_iter 				= 1000;
		int         max_fmg_iter 			= 0;
		int 		bottom_max_iter			= 200;
		int 		max_fixed_iter 			= 500;
		int         verbose 				= 3;
		int         cgverbose 				= 3;
		Set::Scalar tol_rel 				= 1.0E-6;
		Set::Scalar tol_abs 				= 1.0E-6;
		Set::Scalar cg_tol_rel 				= 1.0E-6;
		Set::Scalar cg_tol_abs 				= 1.0E-6;
		Set::Scalar tstart 					= 0.0;
		Set::Scalar tend 					= 0.0;
		std::string bottom_solver			= "bicgstab";
		int 		linop_maxorder 	     	= 2;
		bool	    use_fsmooth 			= false;
		int 		max_coarsening_level	= 0;
		bool 		agglomeration 	  		= true;
		bool 		consolidation 	  		= false;
	} sol;

    struct{
        Set::Vector body_force 				= Set::Vector::Zero();
        Set::Scalar val 					= 0.;
    } loading;

};
}
#endif