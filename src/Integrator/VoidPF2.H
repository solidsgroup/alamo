#ifndef INTEGRATOR_VOIDPF2_H
#define INTEGRATOR_VOIDPF2_H
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
#include "IC/PSRead.H"
#include "IC/PS.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "IC/BMP.H"
#include "IC/EulerAngles.H"
#include "IC/PNG4.H"
#include "IC/PNG3.H"
#include "IC/Quaternions.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Linear/Cubic.H"
#include "Operator/Operator.H"
#include "IC/DynamicNucleation.H"
#include "/usr/include/eigen3/Eigen/Dense"


namespace Integrator
{
template<class MODEL>
class VoidPF2: virtual public Base::Mechanics<MODEL>
{
public:

    VoidPF2(): Base::Mechanics<MODEL>() {}
    VoidPF2(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(VoidPF2& value, IO::ParmParse& pp)
    {
        
        //Base::Mechanics<MODEL>::Parse(value, pp);
        pp_queryclass("mechanics",static_cast<Base::Mechanics<MODEL>&>(value));

        //pp_queryclass("model", value.model);


    

        value.bc_d = new BC::Constant(1,pp,"d.bc.constant");
        value.bc_p = new BC::Nothing();
        value.bc_pr = new BC::Nothing();

        value.RegisterNewFab(value.d_mf, value.bc_d, 1, 3, "d",true);
        value.RegisterNewFab(value.d_old_mf, value.bc_d,  1, 3,"dold",false);
        value.RegisterNewFab(value.p_mf, value.bc_d, 1, 3, "p",true);
        value.RegisterNewFab(value.pr_mf, value.bc_d, 1, 3, "pr",true);
                
        
        

        value.d_on = true;
        
        //value.spall.ic = new IC::PSRead(value.geom, pp, "spall.ic.psread.filename");
              
        
        if (pp.contains("euler.ic.type"))
        {
        
            value.RegisterNodalFab(value.euler_mf, 4, 3, "euler", true);   //newline for euler
            std::string type;
            pp_query("euler.ic.type", type); // Read IC type for the euler angles field
            if (type == "eulerangles")  
            {
                 value.ic_euler = new IC::PNG4(value.geom, pp, "euler.ic.eulerangles");   
                 value.euler_ic_type_flag = 1;
            }                               //newly added line
            else Util::Abort(INFO,"Invalid value for ic.type: ",type);
            
            value.euler_reset_on_regrid = true;
            // Whether to re-initialize eta when re-gridding occurs.
            // Default is false unless eta ic is set, then default is.
            // true.
            pp_query("euler.reset_on_regrid",value.euler_reset_on_regrid);
        }
        
         if (pp.contains("heatmap.ic.type")) {
        value.RegisterNodalFab(value.heatmap_mf, 3, 3, "heatmap", true);
        std::string type;
        pp_query("heatmap.ic.type", type); // Read IC type for heatmap field
        if (type == "hm") {
            value.ic_heatmap = new IC::PNG3(value.geom, pp, "heatmap.ic.hm"); // Specify the file
            value.heatmap_ic_type_flag = 1;
        } else {
            //Util::Abort(INFO, "Invalid value for heatmap.ic.type: ", type);
        }
        }
        
         //if (pp.contains("quat.ic.type"))
        //{
        
            //value.RegisterNodalFab(value.quat_mf, 4, 3, "quat", true);   //newline for quaernions
            //std::string type;
            //pp_query("quat.ic.type", type); // Read IC type for the quaternions angles field
            //if (type == "quaternions")  
            //{
                 //value.ic_quat = new IC::Quaternions(value.geom, pp, "quat.ic.quaternions");   
                 //value.quat_ic_type_flag = 1;
            //}                               //newly added line
            //else Util::Abort(INFO,"Invalid value for ic.type: ",type);
            
            //value.quat_reset_on_regrid = true;
            // Whether to re-initialize eta when re-gridding 
            // Default is false unless eta ic is set, then default is.
            // true.
            //pp_query("quat.reset_on_regrid",value.quat_reset_on_regrid);
        //}
        
        
        // Initialize PSRead with parameters from the input file
        //std::string filename;
        std::string ictype;
        pp_query_validate("spall.ic.type",ictype,{"psread"});
        if (ictype == "psread") value.spall.ic = new IC::PSRead(value.geom, pp,"spall.ic.psread");

        pp_query_required("spall.M",value.spall.M);
        pp_query_required("spall.kappa",value.spall.kappa);
        pp_query_required("spall.cp",value.spall.cp);
        pp_query_required("spall.po",value.spall.po);
        pp_query_required("spall.elastic_mult",value.spall.elastic_mult);
        pp_query_required("spall.small",value.spall.small); // regularization to set minimum modulus multiplier

        pp_query_required("eta_ref_threshold", value.m_eta_ref_threshold); // Refinement threshold based on eta

        
        // Read material properties
        pp.query("C11", value.C11);
        pp.query("C12", value.C12);
        pp.query("C44", value.C44);
        
        
        pp_query_required("nucleation.kb",value.nucleation.kb);
        pp_query_required("nucleation.T",value.nucleation.T);
        pp_query_required("nucleation.threshold",value.nucleation.threshold);
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        //if (d_on)  spall.ic->Initialize(lev,d_mf);
        //if (d_on)  spall.ic->Initialize(lev,d_old_mf);
        
        
       if (d_on)  d_mf[lev]->setVal(0.0);
       if (d_on)  d_old_mf[lev]->setVal(0.0);
        

        p_mf[lev]->setVal(0.0);
        pr_mf[lev]->setVal(0.0);
        
        euler_mf[lev]->setVal(0.0); 
        ic_euler->Initialize(lev, euler_mf);        //newline for euler angles
        
        
        ic_heatmap->Initialize(lev, heatmap_mf);
        
        //quat_mf[lev]->setVal(0.0); 
        //ic_quat->Initialize(lev, quat_mf);        //newline for quaternions
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        BL_PROFILE("VoidPF2::TimeStepBegin");
        Base::Mechanics<MODEL>::TimeStepBegin(a_time, a_step);

        // Constants and parameters
        const Set::Scalar k_B = 1.0;
        const Set::Scalar T = 300.0;            // Temperature
        const Set::Scalar epsilon_e = 1e-6;     // Regularization term to avoid division by zero
        const Set::Scalar nucleation_threshold = 0.8;  // Threshold for nucleation
        const Set::Scalar void_radius = 0.025;   // Radius of the nucleated void
        const Set::Scalar diffuse_thickness = 0.01;
        const Set::Scalar n = 2;

        Set::Scalar E_max = 0.0;  // Variable to store maximum strain energy



        // Correctly passing geometry as argument to DynamicNucleation
        IC::DynamicNucleation dynamic_ic(this->geom);

        // Track if nucleation happens
        bool nucleation_occurred = false;  // Flag to track nucleation

        // Store nucleation points temporarily
        std::vector<std::pair<Set::Scalar, Set::Scalar>> nucleation_points;  // Store (x, y) coordinates


        if (a_step >0) return;
        
        if (a_step == 0)  // Only perform nucleation during the first time step
        {
            // Compute maximum strain energy (E_max)
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar* DX = geom[lev].CellSize();

                for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.tilebox();
                    //amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
                    //amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);

                    amrex::ParallelFor(bx, [=, &E_max] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        // Initialize stress and strain matrices
                        //Eigen::Matrix2d sigma, epsilon;
                        //sigma(0, 0) = sig(i, j, k)(0, 0); sigma(0, 1) = sig(i, j, k)(0, 1);
                        //sigma(1, 0) = sig(i, j, k)(1, 0); sigma(1, 1) = sig(i, j, k)(1, 1);
                        //epsilon(0, 0) = eps(i, j, k)(0, 0); epsilon(0, 1) = eps(i, j, k)(0, 1);
                        //epsilon(1, 0) = eps(i, j, k)(1, 0); epsilon(1, 1) = eps(i, j, k)(1, 1);

                        // Compute strain energy density in 2D
                        //Set::Scalar E_strain = 0.5 * (sigma(0, 0) * epsilon(0, 0) +
                                                      //sigma(1, 1) * epsilon(1, 1) +
                                                      //2.0 * sigma(1, 0) * epsilon(1, 0));

                        // Update maximum strain energy
                        //amrex::Gpu::Atomic::Max(&E_max, E_strain);
                    });
                }
            }

            // Check for nucleation points and create voids
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar* DX = geom[lev].CellSize();

                for (amrex::MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.tilebox();
                    //amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
                    //amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& pr = (*pr_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& d = (*d_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& dold = (*d_old_mf[lev]).array(mfi);
                    amrex::Array4<const Set::Scalar> const& heatmap = (*heatmap_mf[lev]).array(mfi);

                    // Temporarily collect nucleation points from parallel for
                    //amrex::ParallelFor(bx, [=, &dynamic_ic, &nucleation_occurred, E_max, &nucleation_points] AMREX_GPU_DEVICE(int i, int j, int k)
                    amrex::ParallelFor(bx, [=, &dynamic_ic, &nucleation_occurred, &nucleation_points] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        // Initialize stress and strain matrices
                        //Eigen::Matrix2d sigma, epsilon;
                        //sigma(0, 0) = sig(i, j, k)(0, 0); sigma(0, 1) = sig(i, j, k)(0, 1);
                        //sigma(1, 0) = sig(i, j, k)(1, 0); sigma(1, 1) = sig(i, j, k)(1, 1);
                        //epsilon(0, 0) = eps(i, j, k)(0, 0); epsilon(0, 1) = eps(i, j, k)(0, 1);
                        //epsilon(1, 0) = eps(i, j, k)(1, 0); epsilon(1, 1) = eps(i, j, k)(1, 1);

                        // Compute strain energy density
                        //Set::Scalar E_strain = 0.5 * (sigma(0, 0) * epsilon(0, 0) +
                                                      //sigma(1, 1) * epsilon(1, 1) +
                                                      //2.0 * sigma(1, 0) * epsilon(1, 0));

                        // Normalize strain energy by E_max
                        //Set::Scalar E_strain_norm = E_strain / E_max;

                        // Check principal stresses (positive indicates tension)
                        //Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(sigma);
                        //Set::Scalar max_principal_stress = solver.eigenvalues().maxCoeff();

                        //if (max_principal_stress > 0)  // Tension exists, proceed with nucleation
                        {
                            // Compute nucleation energy
                            Set::Scalar R = heatmap(i, j, k, 0);
                            Set::Scalar G = heatmap(i, j, k, 1);
                            Set::Scalar B = heatmap(i, j, k, 2);
                            Set::Scalar H = (R + G + B) / 3.0;

                            Set::Scalar energy = (1.0 - std::max(H, epsilon_e)); // - E_strain_norm)/2;
                            //Set::Scalar p_nucleation = std::exp(-energy / (k_B * T));
                            //p(i, j, k) = std::clamp(p_nucleation, 0.0, 1.0);
     
     
                            pr(i,j,k) = 1.0 - B;

                            // If nucleation probability exceeds the threshold, create a void
                            if (pr(i, j, k) > nucleation_threshold)
                            {
                                // Get the physical coordinates of the nucleation point
                                const Set::Scalar* prob_lo = geom[lev].ProbLo();
                                Set::Scalar center_x = prob_lo[0] + (i + 0.5) * DX[0];
                                Set::Scalar center_y = prob_lo[1] + (j + 0.5) * DX[1];

                                // Add nucleation points to a vector to print later
                                nucleation_points.emplace_back(center_x, center_y);

                                // Add this nucleation point to the dynamic nucleation class
                                Set::Vector nucleation_point(2);
                                nucleation_point(0) = center_x;
                                nucleation_point(1) = center_y;
                                dynamic_ic.AddNucleationPoint(nucleation_point, void_radius);

                                // Mark that nucleation occurred
                                nucleation_occurred = true;
                            }
                        }
                    });
                }

                // Apply the nucleation points to the domain (after the parallel region)
                dynamic_ic.Add(lev, d_mf, 1E-6);
                //std::swap(d_old_mf[lev], d_mf[lev]);
                Util::RealFillBoundary(*d_mf[lev], geom[lev]);
                Util::RealFillBoundary(*d_old_mf[lev], geom[lev]);            
                
            }
        }
        
        

        // Print nucleation points after processing
        if (nucleation_occurred)
        {
            for (const auto& point : nucleation_points)
            {
                std::cout << "Nucleation occurred at physical coordinates (x, y): (" 
                          << point.first << ", " << point.second << ")\n";
            }
        }
        else
        {
            std::cout << "No nucleation occurred at any point during the time step.\n";
        }
    }





    virtual void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        //Base::Mechanics<MODEL>::UpdateModel(a_step, a_time);
        
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        
        //if (a_step > 0) return;

        if (!d_on) return;
        for (int lev = 0; lev <= finest_level; ++lev)
            {
                euler_mf[lev]->FillBoundary();    //euler angles
                
                //quat_mf[lev]->FillBoundary();    //quaternions
                
                heatmap_mf[lev]->FillBoundary(); 
                
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());

                //this->rhs_mf[lev]->setVal(Set::Vector::Zero(0.0));
                d_mf[lev]->FillBoundary();
                //d_old_mf[lev]->FillBoundary();

                // model.Define(C11,C12,C44);
                // model_mf[lev]->setVal(model);
                Set::Vector DX(geom[lev].CellSize());

                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                    {
                        amrex::Array4<MODEL> const& model = model_mf[lev]->array(mfi);
                        amrex::Box bx = mfi.grownnodaltilebox();
                        bx = bx & domain;
                        //Set::Patch<MODEL> model_sum = model_mf.Patch(lev,mfi);
                        Set::Patch<const Set::Scalar> d = d_mf.Patch(lev,mfi);
                        amrex::Array4<const Set::Scalar> const& euler = euler_mf[lev]->array(mfi);
                        
                        //amrex::Array4<const Set::Scalar> const& quat = quat_mf[lev]->array(mfi);
                
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            //auto sten = Numeric::GetStencil(i, j, k, bx);
                            Set::Scalar d_inv_avg = 1.0 - Numeric::Interpolate::CellToNodeAverage(d,i,j,k,0);
                            
                            //Set::Scalar d_inv_avg = 1.0 - d(i,j,k);
                            
                            if (d_inv_avg > 1.0)   d_inv_avg = 1.0;
                            if (d_inv_avg < spall.small) d_inv_avg = spall.small;
                            model(i, j, k) = MODEL::Zero();
                            
                            // Extract quaternion components
                            Set::Scalar w = euler(i, j, k, 0);
                            Set::Scalar x = euler(i, j, k, 1);
                            Set::Scalar y = euler(i, j, k, 2);
                            Set::Scalar z = euler(i, j, k, 3);
                            
                              // Extract quaternion components
                            //Set::Scalar w = quat(i, j, k, 0);
                            //Set::Scalar x = quat(i, j, k, 1);
                            //Set::Scalar y = quat(i, j, k, 2);
                            //Set::Scalar z = quat(i, j, k, 3);
                
                            Eigen::Quaterniond q(w, x, y, z);
                            
                          
                            q.normalize();
                            
                           
                            
                            
                            // ZYX -> yaw (phi1), pitch (Phi), roll (phi2)
                            
                            //Eigen::Matrix3d rotationMatrix = q.toRotationMatrix();
                            
                            //Eigen::Vector3d eulerAngles = rotationMatrix.eulerAngles(2, 1, 0); // ZYX -> yaw (phi1), pitch (Phi), roll (phi2)
                                                       
                            //double phi1 = eulerAngles[0];  // Yaw (Z)
                            //double Phi  = eulerAngles[1];  // Pitch (Y)
                            //double phi2 = eulerAngles[2];  // Roll (X)
                                                     
                            // Step: Use the Euler angles in the model definition
                            //model(i, j, k).Define(C11, C12, C44, phi1, Phi, phi2); // Pass the Euler angles to the model
                        
                            
                            
                            //rotation about ZYZ
                            
                            Eigen::Matrix3d R = q.toRotationMatrix();
                            
                            // Beta (rotation about Y-axis)
                            double beta = std::acos(R(2, 2));  // R33
                            
                            // Alpha (first rotation about Z-axis)
                            double alpha = std::atan2(R(1, 2), R(0, 2));  // R23, R13
                            
                            // Gamma (second rotation about Z-axis)
                            double gamma = std::atan2(R(2, 1), -R(2, 0)); // R32, -R31
                            // Step 6: Define the model with the Euler angles
                            model(i, j, k).Define(C11, C12, C44, alpha, beta, gamma);                            


                            
                            
                            

                            // Calculate Phi
                            //double Phi = acos(2.0 * (w * w + x * x) - 1.0);

                            // Calculate phi1
                            //double phi1 = atan2(z * x - w * y, w * x + y * z);

                            // Calculate phi2
                            //double phi2 = atan2(w * z + x * y, w * y - x * z);

                            // Define model with converted Euler angles
                            //model(i, j, k).Define(C11, C12, C44, phi1, Phi, phi2);
			                
			                
			                
			                //model(i,j,k).Define(C11, C12,C44 ,euler(i,j,k,0)*180.0*(3.14/180.), euler(i,j,k,1)*360.0*(3.14/180.), euler(i,j,k,2)*360.0*(3.14/180.));
			                
			                //model(i,j,k).Define(C11, C12,C44 ,euler(i,j,k,0)*180.0*(3.14/180.), euler(i,j,k,1)*360.0*(3.14/180.), euler(i,j,k,2)*360.0*(3.14/180.));
			                
			                //model(i,j,k).Define(C11, C12,C44,0.0,0.0,0.0);
                            model(i,j,k).ddw *= (d_inv_avg);
                        });
                    }
                Util::RealFillBoundary(*model_mf[lev], this->geom[lev]);
                //Util::RealFillBoundary(*model_mf[lev], geom[lev]);
                Util::RealFillBoundary(*d_mf[lev], geom[lev]);
                Util::RealFillBoundary(*d_old_mf[lev], geom[lev]);

            }

    }

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        BL_PROFILE("VoidPF2::Advance");
        Base::Mechanics<MODEL>::Advance(lev, time, dt);
        
        std::swap(d_old_mf[lev], d_mf[lev]);
        
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();



        //Util::RealFillBoundary(*stress_mf[lev], geom[lev]);
        //Util::RealFillBoundary(*strain_mf[lev], geom[lev]);

        for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(1);
            bx = bx & domain;
            
            amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& dold = (*d_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& d = (*d_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& p = (*p_mf[lev]).array(mfi);
            
            
            //Set::Patcht::Vector grad <const Set::Matrix> sig  = stress_mf.Patch(lev,mfi);
            //Set::Patch<const Set::Matrix> eps  = strain_mf.Patch(lev,mfi);
            //Set::Patch<const Set::Scalar> dold = d_old_mf.Patch(lev,mfi);
            //Set::Patch<Set::Scalar>       d    = d_mf.Patch(lev,mfi);
            //Set::Patch<Set::Scalar>       p    = p_mf.Patch(lev,mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                //auto sten = Numeric::GetStencil(i, j, k, bx);
                Set::Scalar df = 0.0;

                df += spall.cp * 2.0 * dold(i,j,k)*(1 - 2.0*dold(i,j,k))*(1.0 - dold(i,j,k));
                df += spall.kappa * Numeric::Laplacian(dold,i,j,k,0,DX);
                

                Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k, 0);
                Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);
  
                p(i,j,k) = 0.5 * sig_avg.trace();
      
                if (p(i,j,k) > spall.po)
                    df += spall.elastic_mult*(p(i,j,k) - spall.po);

                d(i,j,k) = dold(i,j,k) - spall.M * dt * df;
                
                
                if (df < 0)
                    df = 0.0;
                            
                
                if (d(i,j,k) > 1.0)   d(i,j,k) = 1.0;
                if (d(i,j,k) < 0.0)   d(i,j,k) = 0.0;
                //d(i,j,k) = (eps_avg.transpose()*eps_avg).trace();

            });
        
        }

    }


    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);

           if (d_on)
           {    
            amrex::Array4<Set::Scalar> const& d = d_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                {
                    Set::Vector gradpsi = Numeric::Gradient(d, i, j, k, 0, DX.data(), sten);
                    if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                    //if (gradpsi.lpNorm<2>() * DXnorm > 0.0001)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });

           }
           
        }
        
         {
             for (amrex::MFIter mfi(*euler_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
                amrex::Box bx = mfi.nodaltilebox();
                amrex::Array4<char> const& tags = a_tags.array(mfi);
                amrex::Array4<Set::Scalar> const& euler = euler_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    {
                       Set::Vector grad = Numeric::Gradient(euler, i, j, k, 0, DX.data(), sten);
                        if (grad.lpNorm<2>() * DXnorm > m_euler_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }
            
            
                //for (amrex::MFIter mfi(*quat_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
              //{
                    //amrex::Box bx = mfi.nodaltilebox();
                    //amrex::Array4<char> const& tags = a_tags.array(mfi);
                    //amrex::Array4<Set::Scalar> const& quat = quat_mf[lev]->array(mfi);

                    //amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    //{
                        //auto sten = Numeric::GetStencil(i, j, k, bx);
                        //{
                           //Set::Vector grad = Numeric::Gradient(quat, i, j, k, 0, DX.data(), sten);
                            //if (grad.lpNorm<2>() * DXnorm > m_quat_ref_threshold)
                                //tags(i, j, k) = amrex::TagBox::SET;
                        //}
                    //});
            //}
        }
    }
    
    //void Regrid(int lev, Set::Scalar time) override
    //{
        //ic_euler->Initialize(lev, euler_mf, time);                 //newline added for euler angles
        
        //ic_quat->Initialize(lev, quat_mf, time);                 //newline added for euler angles
    //}



protected:
    MODEL model;
    IC::IC* ic;
    BC::BC<Set::Scalar>*bc_d, *bc_p, *bc_pr;
    Set::Scalar m_eta_ref_threshold = 1E-2;
    Set::Field<Set::Scalar> d_old_mf;
    Set::Field<Set::Scalar> d_mf;
    Set::Field<Set::Scalar> p_mf;
    Set::Field<Set::Scalar> pr_mf;

    bool d_on = true;

    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;

    Set::Field<Set::Scalar> euler_mf;          //new line added
    Set::Scalar m_euler_ref_threshold = 0.01;   //new line added
    IC::IC* ic_euler = nullptr;     //newline added
    int euler_ic_type_flag; // new member variable for the flag
    bool euler_reset_on_regrid = false;
    
    
    //Set::Field<Set::Scalar> quat_mf;          //new line added
    //Set::Scalar m_quat_ref_threshold = 0.01;   //new line added
    //IC::IC* ic_quat = nullptr;     //newline added
    //int quat_ic_type_flag; // new member variable for the flag
    //bool quat_reset_on_regrid = false;
    
    Set::Field<Set::Scalar> heatmap_mf;          //new line added
    Set::Scalar m_heatmap_ref_threshold = 0.01;   //new line added
    IC::IC* ic_heatmap = nullptr;     //newline added
    int heatmap_ic_type_flag; // new member variable for the flag
    bool heatmap_reset_on_regrid = false;
    
    std::string filename; 
    


    struct {
        IC::IC *ic;
        Set::Scalar M = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar po = NAN;
        Set::Scalar cp = NAN;
        Set::Scalar elastic_mult = NAN;
        Set::Scalar small = NAN;
    } spall;
    
    // Material properties
    Set::Scalar C11 = NAN;
    Set::Scalar C12 = NAN;
    Set::Scalar C44 = NAN;
    
    
    struct {    
        Set::Scalar kb = NAN;
        Set::Scalar T  = NAN;
        Set::Scalar threshold = NAN;
        } nucleation;
        
};










} // namespace Integrator
#endif
