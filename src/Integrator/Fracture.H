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
#include "IC/PerturbedInterface.H"

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
        

        pp_crack.queryclass("constant",crack.cracktype);

        IO::ParmParse pp_crack_df("crack.df");
        pp_crack_df.query("mult_Gc", crack.mult_df_Gc);
        pp_crack_df.query("mult_Lap", crack.mult_df_lap);
        pp_crack_df.query("beta", crack.beta);
        pp_crack_df.query("tol_rel",crack.driving_force_tolerance_rel);
        pp_crack_df.query("tol_abs",crack.driving_force_tolerance_abs);

        // IC for crack field
        IO::ParmParse pp("crack.ic");
        pp.queryarr("type", crack.ic_type);
        
        if (crack.ic_type.size() == 0) crack.is_ic = false;
        else
        {
            crack.is_ic = true;
            for (int i = 0; i < crack.ic_type.size(); i++)
            {
                if(crack.ic_type[i] == "notch")
                {
                    IC::Notch *tmpic = new IC::Notch(geom);
                    pp.queryclass("notch",*tmpic);
                    crack.ic.push_back(tmpic);
                }
                else if (crack.ic_type[i] == "ellipsoid")
                {
                    IC::Ellipsoid *tmpic = new IC::Ellipsoid(geom);
                    pp.queryclass("ellipsoid",*tmpic);
                    crack.ic.push_back(tmpic);
                }
                else
                    Util::Abort(INFO, "This IC hasn't been implemented yet");
            }
        }

        RegisterNodalFab(crack.c_mf, 1, number_of_ghost_nodes, "crack", true);
        RegisterNodalFab(crack.c_old_mf, 1, number_of_ghost_nodes, "crack_old", false);
        RegisterNodalFab(crack.driving_force_mf, 5, number_of_ghost_nodes, "driving_force", false);

        RegisterIntegratedVariable(&(crack.driving_force_norm),"driving_force_norm");
        RegisterIntegratedVariable(&(crack.int_crack),"int_c");
        RegisterIntegratedVariable(&(elastic.int_energy),"int_W");

        IO::ParmParse pp_material("material");
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
        }
        else if (material.ic_type == "ellipse")
        {
            IC::Ellipse *tmpic = new IC::Ellipse(geom);
            pp_material_ic.queryclass("ellipse",*tmpic);
            number_of_materials = 2;
            
            material.ic = tmpic;
            material.is_ic = true;
        }
        else if (material.ic_type == "perturbed_interface")
        {
            IC::PerturbedInterface *tmpic = new IC::PerturbedInterface(geom);
            pp_material_ic.queryclass("perturbed_interface", *tmpic);
            number_of_materials = 2;

            material.ic = tmpic;
            material.is_ic = true;
        }
        else
            material.is_ic = false;

        for (int p = 0; p < number_of_materials; p++)
        {
            std::string pp_string = "material.model"+std::to_string(p+1);
            Util::Message(INFO, pp_string);
            IO::ParmParse pp_model(pp_string);

            if(pp_model.contains("type"))
            {
                std::string model_type_str;
                pp_model.query("type", model_type_str);
                if(model_type_str == "isotropic")
                {
                    brittle_fracture_model_type_test *tmpmodel = new brittle_fracture_model_type_test();
                    pp_model.queryclass("isotropic",*tmpmodel);
                    material.brittlemodeltype.push_back(*tmpmodel);
                }
            }
            else
            {
                if (p == 0) Util::Abort(INFO, "Need material information for at least one material");
                material.brittlemodeltype.push_back(material.brittlemodeltype[p-1]);
            }
        }
        
        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);
        pp_load.query("val", loading.val);

        IO::ParmParse pp_elastic("elastic");
        pp_elastic.query("df_mult", elastic.df_mult);
        pp_elastic.queryclass("bc",elastic.brittlebc);
        pp_elastic.query("int",elastic.interval);
        pp_elastic.query("omega",elastic.omega);


        nlevels = maxLevel() + 1;
        RegisterNodalFab(elastic.disp_mf,  AMREX_SPACEDIM, number_of_ghost_nodes, "disp", true);
        RegisterNodalFab(elastic.rhs_mf,  AMREX_SPACEDIM, number_of_ghost_nodes, "rhs", false);
        RegisterNodalFab(elastic.residual_mf,  AMREX_SPACEDIM, number_of_ghost_nodes, "res", false);
        RegisterNodalFab(elastic.strain_mf,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strain", false);
        RegisterNodalFab(elastic.stress_mf,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "stress", true);
        RegisterNodalFab(elastic.energy_mf, 1, number_of_ghost_nodes, "energy", false);
        RegisterNodalFab(elastic.energy_pristine_mf, 1, number_of_ghost_nodes, "energy_pristine", false);
        RegisterNodalFab(elastic.energy_pristine_old_mf, 1, number_of_ghost_nodes, "energy_pristine_old", false);
        RegisterNodalFab(material.material_mf, number_of_materials, number_of_ghost_nodes, "materials", true);
        RegisterGeneralFab(material.model_mf, 1, number_of_ghost_nodes);

    }

protected:
    void Initialize(int ilev) override
    {
        elastic.disp_mf[ilev]->setVal(0.0);
        elastic.strain_mf[ilev]->setVal(0.0);
        elastic.stress_mf[ilev]->setVal(0.0);
        elastic.rhs_mf[ilev]->setVal(0.0);
        elastic.energy_mf[ilev]->setVal(0.0);
        elastic.residual_mf[ilev]->setVal(0.0);
        elastic.energy_pristine_mf[ilev] -> setVal(0.);
        elastic.energy_pristine_old_mf[ilev] -> setVal(0.);

        if(crack.is_ic)
        {
            crack.c_mf[ilev]->setVal(1.0);
            crack.c_old_mf[ilev]->setVal(1.0);

            for (int i = 0; i < crack.ic.size(); i++)
            {
                crack.ic[i]->Add(ilev,crack.c_mf);
                crack.ic[i]->Add(ilev,crack.c_old_mf);
            }
        }
        else
        {
            crack.c_mf[ilev]->setVal(1.0);
            crack.c_old_mf[ilev]->setVal(1.0);
        }

        if(material.is_ic) material.ic->Initialize(ilev,material.material_mf);
        else material.material_mf[ilev]->setVal(1.0);
    }

    void TimeStepBegin(Set::Scalar /*time*/, int /*iter*/) override
    {
        Util::Message(INFO,crack.driving_force_norm," ",crack.driving_force_reference," ",crack.driving_force_tolerance_rel);
        if (crack.driving_force_norm / crack.driving_force_reference < crack.driving_force_tolerance_rel)
            elastic.do_solve_now = true;
        if (crack.driving_force_norm < crack.driving_force_tolerance_abs)
            elastic.do_solve_now = true;

        if (!elastic.do_solve_now) return;
        // if(iter%elastic.interval != 0) return;
        
        Util::Message(INFO);
        for (int p = 0; p < number_of_materials; p++)
            material.brittlemodeltype[p].DegradeModulus(0.0);
        
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            //material.model_mf[ilev]->setVal(material.brittlemodeltype[0]);
            Util::RealFillBoundary(*material.material_mf[ilev],geom[ilev]);
            std::swap(elastic.energy_pristine_old_mf[ilev], elastic.energy_pristine_mf[ilev]);
            
            for (MFIter mfi(*material.material_mf[ilev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<brittle_fracture_model_type_test> const &model = material.model_mf[ilev]->array(mfi);
                amrex::Array4<const Set::Scalar>                const &mat = material.material_mf[ilev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                            Set::Scalar matsum = 0.0;
                                            brittle_fracture_model_type_test modelsum(0.0,0.0);
                                            for (int p = 0; p < number_of_materials; p++)
                                            {
                                                modelsum += mat(i,j,k,p)*material.brittlemodeltype[p];
                                                matsum += mat(i,j,k,p);

                                            }
                                            model(i,j,k) = /*(1./matsum) * */modelsum ;
                                        });
            }
            Util::RealFillBoundary(*material.model_mf[ilev],geom[ilev]);
        }

        // Degrading the modulus
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            Util::RealFillBoundary(*crack.c_mf[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.material_mf[ilev],geom[ilev]);
            
            for (amrex::MFIter mfi(*elastic.disp_mf[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.grownnodaltilebox();
                amrex::Array4<const Set::Scalar> const& c = (*crack.c_mf[ilev]).array(mfi);
                amrex::Array4<brittle_fracture_model_type_test> model = (material.model_mf)[ilev]->array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                                            Set::Scalar _temp = 0.0;
                                            // _temp = crack.cracktype->g_phi(c(i,j,k,0),0);
                                            _temp = crack.cracktype.g_phi(c(i,j,k,0),0) + crack.scaleModulusMax;
                    
                                            if (std::isnan(_temp)) Util::Abort(INFO);
                                            if(_temp < 0.0) _temp = 0.;
                                            if(_temp > 1.0) _temp = 1.0;
                    
                                            // _temp = std::max(crack.scaleModulusMax, _temp);
                                            // _temp = std::max(crack.scaleModulusMax, crack.scaleModulusMax + _temp * (1. - crack.scaleModulusMax));
                                            model(i,j,k).DegradeModulus(1-_temp);
                                        });
            }
            Util::RealFillBoundary(*material.model_mf[ilev],geom[ilev]);
        }

        //
        // Update the body force to be consistent with grid size
        //
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            const Real* DX = geom[ilev].CellSize();
            Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

            AMREX_D_TERM(elastic.rhs_mf[ilev]->setVal(loading.body_force(0)*volume,0,1);,
                        elastic.rhs_mf[ilev]->setVal(loading.body_force(1)*volume,1,1);,
                        elastic.rhs_mf[ilev]->setVal(loading.body_force(2)*volume,2,1););
        }

        elastic.brittlebc.Init(elastic.rhs_mf,geom);

        //
        // Prepare the operator for solve
        //
        Operator::Elastic<brittle_fracture_model_type_test::sym> op_b;
        LPInfo info;
        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp_mf[ilev]->contains_nan()) Util::Warning(INFO);
        
        op_b.define(geom, grids, dmap, info);
        op_b.SetBC(&elastic.brittlebc);
        op_b.SetAverageDownCoeffs(true); // <<< Important! Do not change.
        op_b.SetOmega(elastic.omega);
        //
        // Actually do the elastic solve
        // 
        //for (int ilev = 0; ilev < nlevels; ilev++)
        //{
        //    elastic.disp_mf[ilev]->setVal(0.0);
        //}

        {
            IO::ParmParse pp("elastic");
            Solver::Nonlocal::Newton<brittle_fracture_model_type_test>  solver(op_b);
            //Solver::Nonlocal::Linear<brittle_fracture_model_type_test>  solver(op_b);
            pp.queryclass("solver",solver);
            solver.solve(elastic.disp_mf, elastic.rhs_mf, material.model_mf);
            solver.compResidual(elastic.residual_mf,elastic.disp_mf,elastic.rhs_mf,material.model_mf);
        }
        
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            op_b.Strain(ilev,*elastic.strain_mf[ilev],*elastic.disp_mf[ilev]);
            op_b.Stress(ilev,*elastic.stress_mf[ilev],*elastic.disp_mf[ilev]);
            op_b.Energy(ilev,*elastic.energy_mf[ilev],*elastic.disp_mf[ilev]);

            Util::RealFillBoundary(*elastic.strain_mf[ilev],geom[ilev]);
            Util::RealFillBoundary(*elastic.energy_pristine_old_mf[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.material_mf[ilev],geom[ilev]);
            elastic.energy_pristine_mf[ilev]->setVal(0.0);
            
            for (amrex::MFIter mfi(*elastic.strain_mf[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.grownnodaltilebox();
                amrex::Array4<const Set::Scalar> const& strain  = (*elastic.strain_mf[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& mat = (*material.material_mf[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy        = (*elastic.energy_pristine_mf[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_old    = (*elastic.energy_pristine_old_mf[ilev]).array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    
                                            Set::Matrix eps = Numeric::FieldToMatrix(strain,i,j,k);
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
                                            // brittle_fracture_model_type_test _tmpmodel(0.0,0.0);
                                            // for (int p = 0; p < number_of_materials; p++)
                                            //     _tmpmodel += mat(i,j,k,p)*material.brittlemodeltype[p];


                                            for (int p = 0; p < number_of_materials; p++)
                                                energy(i,j,k) += mat(i,j,k,p)*material.brittlemodeltype[p].W(epsp);

                                            if (energy(i,j,k) < energy_old(i,j,k)) energy(i,j,k) = energy_old(i,j,k);
                                        });
            }
            Util::RealFillBoundary(*elastic.energy_pristine_mf[ilev],geom[ilev]);
        }

        integrate_variables_before_advance = false;
        integrate_variables_after_advance = true;
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        std::swap(crack.c_old_mf[lev], crack.c_mf[lev]);

        //Util::RealFillBoundary(*crack.c_old_mf[lev],geom[lev]);
        //Util::RealFillBoundary(*crack.c_mf[lev],geom[lev]);
        //Util::RealFillBoundary(*material.material_mf[lev],geom[lev]);
        //Util::RealFillBoundary(*elastic.energy_pristine_mf[lev],geom[lev]);

        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain(geom[lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

        for ( amrex::MFIter mfi(*crack.c_mf[lev],true); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.nodaltilebox();
            amrex::Array4<const Set::Scalar> const& c_old = (*crack.c_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& df = (*crack.driving_force_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& c = (*crack.c_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& energy = (*elastic.energy_pristine_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& mat = (*material.material_mf[lev]).array(mfi);
            
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){

#if AMREX_SPACEDIM !=2
                                        Util::Abort(INFO, "This doesn't work with 1D or 3D yet");
#endif
                                        if      (i == lo.x && j == lo.y) c(i,j,k,0) = c(i+1,j+1,k,0);
                                        else if (i == lo.x && j == hi.y) c(i,j,k,0) = c(i+1,j-1,k,0);
                                        else if (i == hi.x && j == lo.y) c(i,j,k,0) = c(i-1,j+1,k,0);
                                        else if (i == hi.x && j == hi.y) c(i,j,k,0) = c(i-1,j-1,k,0);
                                        else if (i == lo.x) c(i,j,k,0) = c(i+1,j,k,0);
                                        else if (j == lo.y) c(i,j,k,0) = c(i,j+1,k,0);
                                        else if (i == hi.x) c(i,j,k,0) = c(i-1,j,k,0);
                                        else if (j == hi.y) c(i,j,k,0) = c(i,j-1,k,0);
                                        else
                                        {
                                            Set::Scalar rhs = 0.0;
                                            //Set::Vector Dc = Numeric::Gradient(c_old, i, j, k, 0, DX);

                                            Set::Scalar bilap = 
                                                Numeric::Stencil<Set::Scalar, 4, 0, 0>::D(c_old,i,j,k,0,DX) + 
                                                2.0 * Numeric::Stencil<Set::Scalar, 2, 2, 0>::D(c_old,i,j,k,0,DX) +
                                                Numeric::Stencil<Set::Scalar, 0, 4, 0>::D(c_old,i,j,k,0,DX); 
                    
                                            Set::Scalar en_cell = energy(i,j,k,0);
                    
                                            //if (std::isnan(en_cell)) Util::Abort(INFO, "Nans detected in en_cell. energy_box(i,j,k,0) = ", energy(i,j,k,0));
                                            //if (std::isinf(c_old(i,j,k,0))) Util::Abort(INFO, "Infs detected in c_old");
                                            //if (c_old(i,j,k,0) > 1) Util::Abort(INFO, "Very large values of c_old at (",i,",",j,",",k,") = ", c_old(i,j,k,0));

                                            df(i,j,k,0) = crack.cracktype.Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;
                                            rhs += crack.cracktype.Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;

                                            //if(rhs > 1.e10) Util::Abort(INFO, "Very large values of rhs at (",i,",",j,",",k,") = ", rhs);

                                            //Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
                                            Set::Scalar laplacian = Numeric::Laplacian(c_old,i,j,k,0,DX);//DDc.trace();

                                            Set::Scalar Gc = crack.cracktype.Gc(0.0);
                                            Set::Scalar Zeta = crack.cracktype.Zeta(0.0);
                                            Set::Scalar _temp_product = 1.0;
                                            for (int m = 0; m < number_of_materials; m++)
                                                _temp_product *= mat(i,j,k,m);
                    
                                            if (number_of_materials > 1) Gc *= (1.0 - _temp_product*(1.-crack.mult_df_Gc)/0.25);

                                            // Util::Message(INFO, "Gc = ", Gc);
                    
                                            df(i,j,k,1) = Gc*crack.cracktype.Dw_phi(c_old(i,j,k,0),0.0)/(4.0*Zeta);
                                            rhs += Gc*crack.cracktype.Dw_phi(c_old(i,j,k,0),0.)/(4.0*Zeta);

                                            df(i,j,k,2) = 2.0*Zeta*Gc*laplacian*crack.mult_df_lap;
                                            rhs -= 2.0*Zeta*Gc*laplacian*crack.mult_df_lap;
                    
                                            df(i,j,k,3) = crack.beta*bilap;
                                            rhs += crack.beta * bilap;
                    
                                            //if (!material.homogeneous) 
                                            //rhs *= crack.cracktype->g_phi(_temp,0.0);

                                            //if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.cracktype->Dw_phi(c_old(i,j,k,0),0.0),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
                                            df(i,j,k,4) =                    std::max(0.,rhs - crack.cracktype.DrivingForceThreshold(c_old(i,j,k,0)));
                                            c(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0.,rhs - crack.cracktype.DrivingForceThreshold(c_old(i,j,k,0)))*crack.cracktype.Mobility(c_old(i,j,k,0));

                                            if (c (i,j,k,0) < 0.0) c(i,j,k,0) = 0.0;
                                            if (c (i,j,k,0) > 1.0) c(i,j,k,0) = 1.0;
                                        }
                                    });
        }
        Util::RealFillBoundary(*crack.c_mf[lev],geom[lev]);
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        const amrex::Real *DX = geom[lev].CellSize();
        const Set::Vector dx(DX);
        const Set::Scalar dxnorm = dx.lpNorm<2>();
        amrex::Box domain(geom[lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());

        for (amrex::MFIter mfi(*crack.c_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box                            bx  = mfi.tilebox() & domain.grow(-1);
            amrex::Array4<char> const                   &tags   = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const      &c  = (*crack.c_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const      &mat    = (*material.material_mf[lev]).array(mfi);
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        Set::Vector grad = Numeric::Gradient(c, i, j, k, 0, DX);
                                        Set::Vector grad2 = Numeric::Gradient(mat, i, j, k, 0, DX);
                                        if (dxnorm * grad.lpNorm<2>() >= crack.refinement_threshold || dxnorm * grad2.lpNorm<2>() >= material.refinement_threshold)
                                            tags(i, j, k) = amrex::TagBox::SET;
                                    });
        }
    }

    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
    {
        const amrex::Real* DX = geom[amrlev].CellSize();
        const Set::Scalar DV = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

        amrex::Array4<const Set::Scalar> const &df = (*crack.driving_force_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Scalar> const &c_new = (*crack.c_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Scalar> const &energy = (*elastic.energy_mf[amrlev]).array(mfi);
        
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                {
                                    crack.driving_force_norm += df(i,j,k,4) * DV;
                                    crack.int_crack += c_new(i,j,k) * DV;
                                    elastic.int_energy += energy(i,j,k) * DV;
                                });
    }

    void TimeStepComplete(amrex::Real /*time*/,int /*iter*/)
    {
        if (elastic.do_solve_now)
            crack.driving_force_reference = crack.driving_force_norm;
        //crack.driving_force_reference = crack.driving_force_norminf;

        elastic.do_solve_now = false;
    }

private:

    int number_of_ghost_nodes = 2;              ///< Number of ghost nodes
    int number_of_materials = 1;
    int nlevels;

    struct{
        Set::Field<Set::Scalar> disp_mf;             ///< displacement field
        Set::Field<Set::Scalar> strain_mf;           ///< total strain field (gradient of displacement)
        Set::Field<Set::Scalar> stress_mf;           ///< stress field
        Set::Field<Set::Scalar> rhs_mf;              ///< rhs fab for elastic solution
        Set::Field<Set::Scalar> residual_mf;         ///< residual field for solver
        Set::Field<Set::Scalar> energy_mf;           ///< total elastic energy
        Set::Field<Set::Scalar> energy_pristine_mf;      ///< energy of the prisitne material as if no crack is present
        Set::Field<Set::Scalar> energy_pristine_old_mf;  ///< energy of the pristine material for previous time step.
        Set::Scalar int_energy = 0.0;               ///< integrated energy over the entire domain.

        BC::Operator::Elastic::Constant     brittlebc;  ///< elastic BC if using brittle fracture
        Set::Scalar df_mult = 1.0;              ///< mulitplier for elastic driving force.
        bool do_solve_now = false;
        int interval = 0;
        Set::Scalar omega = 2./3.;
    } elastic;

    

    struct{
        Set::Field<Set::Scalar> c_mf;                ///< crack field at current time step
        Set::Field<Set::Scalar> c_old_mf;            ///< crack field at previous time step
        Set::Field<Set::Scalar> driving_force_mf;         ///< crack driving forces.
        Set::Scalar int_crack = 0.0;                    ///< integrated crack field over the entire domain.
        Set::Scalar driving_force_reference = 1.0;
        Set::Scalar driving_force_norm = 0.0;
        Set::Scalar driving_force_norminf = 0.0;
        Set::Scalar driving_force_tolerance_rel = 1E-4;
        Set::Scalar driving_force_tolerance_abs = 0.0;

        Model::Interface::Crack::Constant cracktype;                      ///< type of crack. See Crack/Constant or Crack/Sin
        amrex::Vector<std::string> ic_type;                            ///< crack IC type. See IC/Notch and IC/Ellipsoid
        amrex::Vector<IC::IC*> ic;                                     ///< crack IC. See IC/Notch and IC/Ellipsoid
        bool is_ic = false;

        Set::Scalar scaleModulusMax = 0.02;         ///< material modulus ratio inside crack (default = 0.02).
        Set::Scalar refinement_threshold = 0.001;

        Set::Scalar mult_df_Gc = 1.0;               ///< Multiplier for Gc/zeta term
        Set::Scalar mult_df_lap = 1.0;              ///< Multiplier for the laplacian term
        Set::Scalar beta = 0.0;
    } crack;

    struct{
        amrex::Vector<brittle_fracture_model_type_test> brittlemodeltype;
        Set::Field<Set::Scalar> material_mf;
        Set::Field<brittle_fracture_model_type_test> model_mf;
        std::string input_material = "isotropic";
        Set::Scalar refinement_threshold = 0.1;
        
        IC::IC *ic;
        std::string ic_type;
        bool is_ic = false;
    } material;

    struct{
        Set::Vector body_force              = Set::Vector::Zero();
        Set::Scalar val                     = 0.;
    } loading;

};
}
#endif
