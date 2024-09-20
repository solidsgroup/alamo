#ifndef INTEGRATOR_TENSION_H
#define INTEGRATOR_TENSION_H
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


namespace Integrator
{
template<class MODEL>
class Tension: virtual public Base::Mechanics<MODEL>
{
public:

    Tension(): Base::Mechanics<MODEL>() {}
    Tension(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(Tension& value, IO::ParmParse& pp)
    {
        
        //Base::Mechanics<MODEL>::Parse(value, pp);
        pp_queryclass("mechanics",static_cast<Base::Mechanics<MODEL>&>(value));


              
        
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



        pp_query_required("eta_ref_threshold", value.m_eta_ref_threshold); // Refinement threshold based on eta

        
        // Read material properties
        pp.query("C11", value.C11);
        pp.query("C12", value.C12);
        pp.query("C44", value.C44);
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);

        
        euler_mf[lev]->setVal(0.0); 
        ic_euler->Initialize(lev, euler_mf);        //newline for euler angles

    }




    //void TimeStepBegin(Set::Scalar a_time, int a_step) override
    //{ 
    //    Base::Mechanics<MODEL>::TimeStepBegin(a_time, a_step);
    //}

    virtual void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        //Base::Mechanics<MODEL>::UpdateModel(a_step, a_time);
        
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        
        if (a_step > 0) return;

        for (int lev = 0; lev <= finest_level; ++lev)
            {
                euler_mf[lev]->FillBoundary();    //euler angles
                
                //quat_mf[lev]->FillBoundary();    //quaternions
                
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());

                //this->rhs_mf[lev]->setVal(Set::Vector::Zero(0.0));
  
                // model.Define(C11,C12,C44);
                // model_mf[lev]->setVal(model);
                Set::Vector DX(geom[lev].CellSize());

                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                    {
                        amrex::Array4<MODEL> const& model = model_mf[lev]->array(mfi);
                        amrex::Box bx = mfi.grownnodaltilebox();
                        bx = bx & domain;
                        //Set::Patch<MODEL> model_sum = model_mf.Patch(lev,mfi);
                        amrex::Array4<const Set::Scalar> const& euler = euler_mf[lev]->array(mfi);
                        
                        //amrex::Array4<const Set::Scalar> const& quat = quat_mf[lev]->array(mfi);
                
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            //auto sten = Numeric::GetStencil(i, j, k, bx);


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

                        });
                    }
                Util::RealFillBoundary(*model_mf[lev], this->geom[lev]);

            }

    }



    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        
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
            

        }
    }
    
    void Regrid(int lev, Set::Scalar time) override
    {
        ic_euler->Initialize(lev, euler_mf, time);                 //newline added for euler angles
        
        //ic_quat->Initialize(lev, quat_mf, time);                 //newline added for euler angles
    }



protected:
    MODEL model;
    IC::IC* ic;
    Set::Scalar m_eta_ref_threshold = 1E-2;
    Set::Field<Set::Scalar> p_mf;



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

    
    // Material properties
    Set::Scalar C11 = NAN;
    Set::Scalar C12 = NAN;
    Set::Scalar C44 = NAN;
};










} // namespace Integrator
#endif
