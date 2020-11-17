#ifndef INTEGRATOR_FRACTURE_H
#define INTEGRATOR_FRACTURE_H

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
#include "IC/Laminate.H"

#include "Operator/Operator.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"
//#include "Model/Solid/Affine/J2PlasticDegradable.H"
//#include "Model/Solid/Affine/CrystalPlasticDegradable.H"

#include "Model/Interface/Crack/Crack.H"
#include "Model/Interface/Crack/Constant.H"
#include "Model/Interface/Crack/Sin.H"

#include "Numeric/Stencil.H"
#include <eigen3/Eigen/Dense>

namespace Integrator
{
//using brittle_fracture_model_type_test = Model::Solid::Linear::Isotropic;
using brittle_fracture_model_type_test = Model::Solid::Linear::IsotropicDegradable;

class Fracture : public Integrator
{

public:
    Fracture()
    {
        IO::ParmParse pp_crack("crack");
        pp_crack.query("modulus_scaling_max",crack.scaleModulusMax);
        pp_crack.query("refinement_threshold",crack.refinement_threshold);
        pp_crack.query("df_tol_rel",crack.driving_force_tolerance_rel);
        pp_crack.query("df_tol_abs",crack.driving_force_tolerance_abs);

        Model::Interface::Crack::Constant *tmpbdy = new Model::Interface::Crack::Constant();
        pp_crack.queryclass("constant",*tmpbdy);
        crack.cracktype = tmpbdy;

        IO::ParmParse pp_crack_df("crack.df");
        pp_crack_df.query("mult_Gc", crack.mult_df_Gc);
        pp_crack_df.query("mult_Lap", crack.mult_df_lap);

        {
            IO::ParmParse pp("pfcrack");
            pp.query("beta",crack.beta);
        }


        // IC for crack field
        IO::ParmParse pp("crack.ic");
        pp.query("type", crack.ic_type);
        if (crack.ic_type == "notch")
        {
            IC::Notch *tmpic = new IC::Notch(geom);
            pp.queryclass("notch",*tmpic);
            crack.ic = tmpic;
            crack.is_ic = true;
        }
        else
            crack.is_ic = false;

        RegisterNodalFab(crack.field, 1, number_of_ghost_nodes, "crack", true);
        RegisterNodalFab(crack.field_old, 1, number_of_ghost_nodes, "crack_old", true);
        RegisterNodalFab(crack.driving_force, 4, number_of_ghost_nodes, "driving_force", true);
        RegisterIntegratedVariable(&(crack.driving_force_norm),"driving_force_norm");

        IO::ParmParse pp_material("material");
        pp_material.queryclass("isotropic",material.brittlemodeltype);
        pp_material.query("refinement_threshold",material.refinement_threshold);

        IO::ParmParse pp_material_ic("material.ic");
        pp_material_ic.query("type",material.ic_type);
        if (material.ic_type == "laminate")
        {
            IC::Laminate *tmpic = new IC::Laminate(geom);
            pp_material_ic.queryclass("laminate",*tmpic);

            IO::ParmParse _temp("material.ic.laminate");
            _temp.query("number_of_inclusions",number_of_materials);
            number_of_materials++;

            material.ic = tmpic;
            material.is_ic = true;
            material.homogeneous = true;
            material.interface = true;
        }
        else if (material.ic_type == "ellipsoid")
        {
            IC::Ellipsoid *tmpic = new IC::Ellipsoid(geom);
            pp_material_ic.queryclass("ellipsoid",*tmpic);
            material.ic = tmpic;
            material.is_ic = true;
            material.homogeneous = false;
        }
        else if (material.ic_type == "notch")
        {
            IC::Notch *tmpic = new IC::Notch(geom);
            pp_material_ic.queryclass("notch",*tmpic);
            material.ic = tmpic;
            material.is_ic = true;
            material.homogeneous = false;
        }
        else
            material.is_ic = false;

        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);
        pp_load.query("val", loading.val);

        IO::ParmParse pp_elastic("elastic");
        pp_elastic.query("df_mult", elastic.df_mult);
        pp_elastic.queryclass("bc",elastic.brittlebc);
        pp_elastic.query("int",elastic.interval);


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
        RegisterNodalFab(material.modulus_field, number_of_materials, number_of_ghost_nodes, "modulus_field", true);
        RegisterGeneralFab(material.brittlemodel, 1, number_of_ghost_nodes);

    }

protected:
    void Initialize(int ilev) override
    {
        elastic.disp[ilev]->setVal(0.0);
        elastic.strain[ilev]->setVal(0.0);
        elastic.stress[ilev]->setVal(0.0);
        elastic.rhs[ilev]->setVal(0.0);
        elastic.energy[ilev]->setVal(0.0);
        elastic.residual[ilev]->setVal(0.0);
        elastic.energy_pristine[ilev] -> setVal(0.);
        elastic.energy_pristine_old[ilev] -> setVal(0.);

        //crack.ic->Initialize(ilev,material.modulus_field);
        if(material.is_ic) material.ic -> Initialize(ilev,material.modulus_field);
        else material.modulus_field[ilev]->setVal(1.0);

        if(crack.is_ic)
        {
            crack.ic->Initialize(ilev,crack.field);
            crack.ic->Initialize(ilev,crack.field_old);
        }
        else
        {
            crack.field[ilev]->setVal(1.0);
            crack.field_old[ilev]->setVal(1.0);
        }

        material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
    }

    void TimeStepBegin(Set::Scalar /*time*/, int iter) override
    {
        crack.driving_force_norminf = 0.0;
        for (int i = 0; i < crack.driving_force.size(); i++)
        {
            crack.driving_force_norminf = std::max(crack.driving_force_norminf, crack.driving_force[i]->norminf(3));
        }

        Util::Message(INFO,crack.driving_force_norminf, " ", crack.driving_force_norm," ",crack.driving_force_reference," ",crack.driving_force_tolerance_rel);
        if(elastic.interval > 0) 
        {  
            if (iter%elastic.interval) return;
        }
        else
        {
            if (crack.driving_force_norminf / crack.driving_force_reference < crack.driving_force_tolerance_rel)
                elastic.do_solve_now = true;
            if (crack.driving_force_norminf < crack.driving_force_tolerance_abs)
                elastic.do_solve_now = true;
            if (!elastic.do_solve_now) return;
            //if (crack.driving_force_norm / crack.driving_force_reference < crack.driving_force_tolerance_rel)
            //    elastic.do_solve_now = true;
            //if (crack.driving_force_norm < crack.driving_force_tolerance_abs)
            //    elastic.do_solve_now = true;
            //if (!elastic.do_solve_now) return;
        }
        
        material.brittlemodeltype.DegradeModulus(0.0);
        
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            std::swap(elastic.energy_pristine_old[ilev], elastic.energy_pristine[ilev]);
            material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
            Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }

        // Degrading the modulus
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            Util::RealFillBoundary(*crack.field[ilev],geom[ilev]);
            Util::RealFillBoundary(*crack.field_old[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.modulus_field[ilev],geom[ilev]);
            
            for (amrex::MFIter mfi(*elastic.disp[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.grownnodaltilebox();
                amrex::Array4<const Set::Scalar> const& c_new = (*crack.field[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[ilev]).array(mfi);
                amrex::Array4<brittle_fracture_model_type_test> modelfab = (material.brittlemodel)[ilev]->array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Scalar _temp = 0;
                    if (!material.homogeneous)
                    	_temp = std::min (crack.cracktype->g_phi(c_new(i,j,k,0),0), crack.cracktype->g_phi(modbox(i,j,k,0),0));
                    else
                    	_temp = crack.cracktype->g_phi(c_new(i,j,k,0),0);
                    if (std::isnan(_temp)) Util::Abort(INFO);
                    if(_temp < 0.0) Util::Abort(INFO);//_temp = 0.;
                    if(_temp > 1.0) Util::Abort(INFO);//_temp = 1.0;
					_temp = crack.scaleModulusMax + _temp * (1. - crack.scaleModulusMax);
                    //modelfab(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-crack.scaleModulusMax));
                    modelfab(i,j,k,0).DegradeModulus(1-_temp);
                });
            }
            Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }

        //
        // Update the body force to be consistent with grid size
        //
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            const Real* DX = geom[ilev].CellSize();
            Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

            AMREX_D_TERM(elastic.rhs[ilev]->setVal(loading.body_force(0)*volume,0,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(1)*volume,1,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(2)*volume,2,1););
        }

        elastic.brittlebc.Init(elastic.rhs,geom);

        //
        // Prepare the operator for solve
        //
        Operator::Elastic<brittle_fracture_model_type_test::sym> op_b;
        LPInfo info;
        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp[ilev]->contains_nan()) Util::Warning(INFO);
        Util::Message(INFO,geom.size());
        op_b.define(geom, grids, dmap, info);
        op_b.SetBC(&elastic.brittlebc);
        
        //
        // Actually do the elastic solve
        // 
        {
            IO::ParmParse pp("elastic");
            Solver::Nonlocal::Newton<brittle_fracture_model_type_test>  solver(op_b);
            pp.queryclass("solver",solver);
            //if (iter > 4585) solver.setFixedIter(10);
        
            Util::Message(INFO,elastic.disp.size());
            Util::Message(INFO,elastic.rhs.size());
    
            solver.solve(elastic.disp, elastic.rhs, material.brittlemodel);
            solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,material.brittlemodel);
            
            //if (iter > 4585) 
            //{
            //    WritePlotFile();
            //    Util::Abort(INFO);
            //}
        }
        
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            op_b.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
            op_b.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
            op_b.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);

            Util::RealFillBoundary(*elastic.strain[ilev],geom[ilev]);
            Util::RealFillBoundary(*elastic.energy_pristine_old[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.modulus_field[ilev],geom[ilev]);
            elastic.energy_pristine[ilev]->setVal(0.0);
            
            for (amrex::MFIter mfi(*elastic.strain[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.grownnodaltilebox();
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
                    if (!material.homogeneous)
                    {
                    	Set::Scalar _temp = crack.cracktype->g_phi(modbox(i,j,k,0),0.);
                    	_temp = crack.scaleModulusMax + _temp * (1. - crack.scaleModulusMax);
                    	material.brittlemodeltype.DegradeModulus(1.- _temp);
                    }
                    energy_box(i,j,k,0) = material.brittlemodeltype.W(epsp);
                    energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
                });
            }
            Util::RealFillBoundary(*elastic.energy_pristine[ilev],geom[ilev]);
        }

        integrate_variables_before_advance = false;
        integrate_variables_after_advance = true;
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        std::swap(crack.field_old[lev], crack.field[lev]);

        // crack.field_old[lev]->FillBoundary();
        // material.modulus_field[lev]->FillBoundary();
        Util::RealFillBoundary(*crack.field_old[lev],geom[lev]);
        Util::RealFillBoundary(*crack.field[lev],geom[lev]);
        Util::RealFillBoundary(*material.modulus_field[lev],geom[lev]);
        Util::RealFillBoundary(*elastic.energy_pristine[lev],geom[lev]);

        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain(geom[lev].Domain());
        // amrex::Box domain(geom.Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

        for ( amrex::MFIter mfi(*crack.field[lev],true); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& c_old = (*crack.field_old[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& df = (*crack.driving_force[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& c_new = (*crack.field[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& energy_box = (*elastic.energy_pristine[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[lev]).array(mfi);
            
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                
#if AMREX_SPACEDIM !=2
                Util::Abort(INFO, "This doesn't work with 1D or 3D yet");
#endif
                if (i == lo.x && j == lo.y) c_new(i,j,k,0) = c_new(i+1,j+1,k,0);
                else if (i == lo.x && j == hi.y) c_new(i,j,k,0) = c_new(i+1,j-1,k,0);
                else if (i == hi.x && j == lo.y) c_new(i,j,k,0) = c_new(i-1,j+1,k,0);
                else if (i == hi.x && j == hi.y) c_new(i,j,k,0) = c_new(i-1,j-1,k,0);
                else if (i == lo.x) c_new(i,j,k,0) = c_new(i+1,j,k,0);
                else if (j == lo.y) c_new(i,j,k,0) = c_new(i,j+1,k,0);
                else if (i == hi.x) c_new(i,j,k,0) = c_new(i-1,j,k,0);
                else if (j == hi.y) c_new(i,j,k,0) = c_new(i,j-1,k,0);
                else
                {
                    Set::Scalar rhs = 0.0;
                    //Set::Vector Dc = Numeric::Gradient(c_old, i, j, k, 0, DX);

                    Set::Scalar bilap = 
                              Numeric::Stencil<Set::Scalar,4,0,0>::D(c_old,i,j,k,0,DX) +
                        2.0 * Numeric::Stencil<Set::Scalar,2,2,0>::D(c_old,i,j,k,0,DX) + 
                              Numeric::Stencil<Set::Scalar,0,4,0>::D(c_old,i,j,k,0,DX);
                    
                    Set::Scalar _temp = modbox(i,j,k,0);
                    Set::Scalar en_cell = energy_box(i,j,k,0);
                    
                    if (std::isnan(en_cell)) Util::Abort(INFO, "Nans detected in en_cell. energy_box(i,j,k,0) = ", energy_box(i,j,k,0));
                    if (std::isinf(c_old(i,j,k,0))) Util::Abort(INFO, "Infs detected in c_old");
                    //if (c_old(i,j,k,0) > 1) Util::Abort(INFO, "Very large values of c_old at (",i,",",j,",",k,") = ", c_old(i,j,k,0));

                    df(i,j,k,0) = crack.cracktype->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;
                    rhs += crack.cracktype->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;

                    //if(rhs > 1.e10) Util::Abort(INFO, "Very large values of rhs at (",i,",",j,",",k,") = ", rhs);

                    Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
                    Set::Scalar laplacian = DDc.trace();

                    Set::Scalar Gc = crack.cracktype->Gc(0.0);
                    if (material.interface)
                    {
                    	Set::Scalar _temp_product = 1.0;
                    	for (int m = 0; m < number_of_materials; m++)
                        	_temp_product *= modbox(i,j,k,m);
                    	Gc *= (1.0 - _temp_product);
                    }

                    df(i,j,k,1) = Gc*crack.cracktype->Dw_phi(c_old(i,j,k,0),0.0)/(4.0*crack.cracktype->Zeta(0.0))*crack.mult_df_Gc;
                    df(i,j,k,2) = 2.0*crack.cracktype->Zeta(0.0)*Gc*laplacian*crack.mult_df_lap;

                    rhs += Gc*crack.cracktype->Dw_phi(c_old(i,j,k,0),0.)/(4.0*crack.cracktype->Zeta(0.0))*crack.mult_df_Gc;

                    rhs -= 2.0*crack.cracktype->Zeta(0.0)*Gc*laplacian*crack.mult_df_lap;
                    
                    if (!material.homogeneous) 
                    	rhs *= crack.cracktype->g_phi(_temp,0.0);

                    rhs += crack.beta*bilap;

                    df(i,j,k,3) = std::max(0.,rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)));
                    if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.cracktype->Dw_phi(c_old(i,j,k,0),0.0),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
                    c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0., rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)))*crack.cracktype->Mobility(c_old(i,j,k,0));

                    if (c_new (i,j,k,0) < 0.0) c_new(i,j,k,0) = 0.0;
                    if (c_new (i,j,k,0) > 1.0) c_new(i,j,k,0) = 1.0;
                }
            });
        }
        Util::RealFillBoundary(*crack.field[lev],geom[lev]);
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        const amrex::Real *DX = geom[lev].CellSize();
        const Set::Vector dx(DX);
        const Set::Scalar dxnorm = dx.lpNorm<2>();
        amrex::Box domain(geom[lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());

        for (amrex::MFIter mfi(*crack.field[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box 							bx 	= mfi.tilebox() & domain.grow(-1);
            amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const 		&c_new 	= (*crack.field[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const 		&modbox 	= (*material.modulus_field[lev]).array(mfi);
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX);
                Set::Vector grad2 = Numeric::Gradient(modbox, i, j, k, 0, DX);
                if (dxnorm * grad.lpNorm<2>() >= crack.refinement_threshold || dxnorm * grad2.lpNorm<2>() >= material.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
    {
        const amrex::Real* DX = geom[amrlev].CellSize();
        const Set::Scalar DV = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

        amrex::Array4<const Set::Scalar> const &df = (*crack.driving_force[amrlev]).array(mfi);
        
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
        {
            crack.driving_force_norm += df(i,j,k,3) * DV;
        });
    }

    void TimeStepComplete(amrex::Real /*time*/,int /*iter*/)
    {
        if (elastic.do_solve_now)
            //crack.driving_force_reference = crack.driving_force_norm;
            crack.driving_force_reference = crack.driving_force_norminf;

        elastic.do_solve_now = false;
    }

private:

    int number_of_ghost_nodes = 2;				///< Number of ghost nodes
    int number_of_materials = 1;
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

        BC::Operator::Elastic::Constant 	brittlebc;  ///< elastic BC if using brittle fracture
		Set::Scalar df_mult = 1.0;              ///< mulitplier for elastic driving force.
        bool do_solve_now = false;
        int interval = 0;
    } elastic;

    struct{
        Set::Field<Set::Scalar> field;                ///< crack field at current time step
        Set::Field<Set::Scalar> field_old;            ///< crack field at previous time step
        Set::Field<Set::Scalar> driving_force;         ///< crack driving forces.
        Set::Scalar driving_force_reference = 1.0;
        Set::Scalar driving_force_norm = 0.0;
        Set::Scalar driving_force_norminf = 0.0;
		Set::Scalar driving_force_tolerance_rel = 1E-4;
		Set::Scalar driving_force_tolerance_abs = 0.0;

        Model::Interface::Crack::Crack *cracktype;       ///< type of crack. See Crack/Constant or Crack/Sin
		std::string ic_type;                            ///< crack IC type. See IC/Notch and IC/Ellipsoid
		IC::IC *ic;                                     ///< crack IC. See IC/Notch and IC/Ellipsoid
        bool is_ic = false;

        Set::Scalar scaleModulusMax = 0.02;         ///< material modulus ratio inside crack (default = 0.02).
        Set::Scalar refinement_threshold = 0.001;

		Set::Scalar mult_df_Gc = 1.0;				///< Multiplier for Gc/zeta term
		Set::Scalar mult_df_lap = 1.0;				///< Multiplier for the laplacian term
        Set::Scalar beta = 1.0;
    } crack;

    struct{
        brittle_fracture_model_type_test brittlemodeltype;
        Set::Field<Set::Scalar> modulus_field;
        Set::Field<brittle_fracture_model_type_test> brittlemodel;
        std::string input_material = "isotropic";
        Set::Scalar refinement_threshold = 0.1;
        
        IC::IC *ic;
        std::string ic_type;
        bool is_ic = false;
        bool homogeneous = true;
	bool interface = false;
    } material;

    struct{
        Set::Vector body_force 				= Set::Vector::Zero();
        Set::Scalar val 					= 0.;
    } loading;

};
}
#endif
