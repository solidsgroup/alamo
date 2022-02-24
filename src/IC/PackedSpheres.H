#ifndef IC_PACKEDSPHERES_H_
#define IC_PACKEDSPHERES_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "mpi.h"
#include "IO/ParmParse.H"

namespace IC
{
/// \brief Create a random spheres with different radius
///
///
/// Let `N` be the number of patches in the tesselation,
/// and let \f$\alpha_1,\ldots,\alpha_n\f$ be a set of scalar values.
/// There are two methods for applying the tesselation IC:
/// Here partition method is used for applying the Packed spheres IC:

/// This requires an `N`-component multifab.
/// (There is currently no check for this.)
///  Each component of the fab corresponds to one of the patches
///  Then fab component i is equal to \f$\alpha_i\f$ within the non ovelapped sphere  and
///  zero elsewhere. (This is mostly useful for multiphase field
///   microstructure)

///
/// Note on periodicity: This IC uses the `amrex::Geometry` object to determine
/// if it is periodic, and if so, it will generate a periodic spheres.
///
class PackedSpheres : public IC
{
public:
    enum Type
    {
        Partition,
        Values
    };

    PackedSpheres(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}

    void Define()
    {
        std::vector<double> sph_vol;
        std::vector<Set::Scalar> sph_nos;
        Set::Scalar sph_tot_vol = 0;
        R.resize(interval);
        R_n.resize(interval);
        sph_vol.resize(interval);
        sph_nos.resize(interval);
        R_new.resize(interval);
        pack_size = 0;
        sph_tot_vol = 0;
        std::vector<double> ins(1, 1);
        int pack_size1=0;
        int pack_size2=0;
        int pack_size3=0;
        std::vector<Set::Scalar> V_per{V_per1,V_per2,V_per3,V_per4};
        std::vector<Set::Scalar> R_per{R_per1,R_per2,R_per3,R_per4};
            

        int p = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Set::Vector size;
        //size of the domain in 3-D
            
        AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
                     size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,       
                     size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)
            int total = 0;
  
        // converting volume fraction in to no of spheres in 3-D 
#if AMREX_SPACEDIM > 2
        if (geom[0].isPeriodic(2))
        {
            double dom_vol = size(0) * size(1) * size(2);
            for (int n = 0; n < interval; n++)
            {
                Util::Message(INFO, "3 dimensions");
                if (n == 0)
                {
                    // R[n]=R_min+((R_max-R_min)/interval);
                    R[n] = ((R_mean) /interval);
                    sph_vol[n] = (1.33) * 3.14 * verbosity * (R[n] * R[n] * R[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size = pack_size + sph_nos[n];
                }

                else if (n > 0)
                {
                    R[n] = ((R_mean) / interval);
                    sph_vol[n] = (1.33) * 3.14 * verbosity * (R[n] * R[n] * R[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size = pack_size + sph_nos[n];
                }
            }
        }
//#endif

// converting volume fraction in to no of spheres in 2-D 
#elif AMREX_SPACEDIM > 1
        if (geom[0].isPeriodic(1))
        {
            double dom_vol = size(0) * size(1);
            Util::Message(INFO, "dom_vol ",dom_vol);
            // converting volume fraction in to no of spheres
            for (int n = 0; n < interval; n++)  
            {
                Util::Message(INFO, "2 dimensions");
                if (n == 0)
                {
                    // R[n]=R_min+((R_max-R_min)/interval);
                    R[n] = ((R_per[n]*R_mean));
                    sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                    // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil(((V_per[n]*dom_vol) / (sph_vol[n] * 1)) * vf);
                    //Util::Message(INFO, "sph_nos ",sph_nos[n]);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size1= sph_nos[n];
                }

                else if (n == 1)
                {
                    R[n] = ((R_per[n]*R_mean));
                    sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                    // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil(((V_per[n]*dom_vol) / (sph_vol[n] * 1)) * vf);
                    // Util::Message(INFO, "sph_nos ",sph_nos[n]);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size2 =  pack_size1+sph_nos[n];
                }
                    
                else if (n == 2)
                {
                    R[n] = ((R_per[n]*R_mean));
                    sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                    // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil(((V_per[n]*dom_vol) / (sph_vol[n] * 1)) * vf);
                    // Util::Message(INFO, "sph_nos ",sph_nos[n]);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size3 =  pack_size2+ sph_nos[n];
                }
                    
                else if (n == 3)
                {
                    R[n] = ((R_per[n]*R_mean));
                    sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                    // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                    // sph_nos[n]=ceil(ceil(dom_vol/interval/(sph_vol[n]*1))*1);
                    sph_nos[n] = ceil(((V_per[n]*dom_vol) / (sph_vol[n] * 1)) * vf);
                    // Util::Message(INFO, "sph_nos ",sph_nos[n]);
                    sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                    pack_size =  pack_size3+sph_nos[n];
                }
                    
                   
            }
        }
#endif

             
        R_new.resize(pack_size);
            
        if (rank == 0)
        {
                
            int n=0;
            while(n<pack_size)
            {
                if(n<pack_size1)
                {
                    R_new[n]=R_per1*R_mean;     
                    
                }
                
                if((n>=pack_size1)&&(n<pack_size2))
                {
                    R_new[n]=R_per2*R_mean;
                    
                }
                
                if((n>=pack_size2)&&(n<pack_size3))
                {
                    R_new[n]=R_per3*R_mean;
                    
                }
                if((n>=pack_size3)&&(n<pack_size))
                {
                    R_new[n]=R_per4 *R_mean;
                    
                }
                n++;
            }

            int oLoop=0 ; // max no of times this counter checks for the new sphere centers which satisfies non overlap condition 

          
            dd.resize(pack_size - 1);
            dd1.resize(pack_size - 1);
            R_n.resize(pack_size - 1);
            packedspheres.resize(pack_size);
            packedspheres1.resize(pack_size);
            R_new.resize(pack_size);
            count.resize(pack_size);

            for (int n = 0; n < pack_size; n++)
            {

                std::fill(count.begin(), count.end(), 0);
                //first spheres center point 
                if (n == 0)
                {
                    AMREX_D_TERM(packedspheres1[0](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0) - geom[0].ProbLo(0)) * Util::Random();,
                                 packedspheres1[0](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1) - geom[0].ProbLo(1)) * Util::Random();,
                                 packedspheres1[0](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2) - geom[0].ProbLo(2)) * Util::Random(););
                    packedspheres[n] = packedspheres1[n];
                }
                // center points for the second and more spheres
                else if (n > 0)
                {

                    AMREX_D_TERM(packedspheres1[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0) - geom[0].ProbLo(0)) * Util::Random();,
                                 packedspheres1[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1) - geom[0].ProbLo(1)) * Util::Random();,
                                 packedspheres1[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2) - geom[0].ProbLo(2)) * Util::Random(););

                    p = 0;

                    while (p < n)
                    {

                        dd[p] = (packedspheres1[n] - packedspheres1[p]).lpNorm<2>(); //distance between spheres centers

                        //caculating the distance between sphere centers for the periodic boundary condtion
                        if (geom[0].isPeriodic(0))
                        {
                            dd1[p] = std::min((packedspheres1[n] - packedspheres1[p] + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                                              (packedspheres1[n] - packedspheres1[p] - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                        }
#if AMREX_SPACEDIM > 1
                        if (geom[0].isPeriodic(1))
                        {
                            dd1[p] = std::min((packedspheres1[n] - packedspheres1[p] + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                                              (packedspheres1[n] - packedspheres1[p] - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                        }
#endif
#if AMREX_SPACEDIM > 2
                        if (geom[0].isPeriodic(2))
                        {
                            dd1[p] = std::min((packedspheres1[n] - packedspheres1[p] + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                                              (packedspheres1[n] - packedspheres1[p] - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                        }
#endif

                        dd[p] = std::min(dd[p], dd1[p]);   

                        R_n[p] = R_new[n] + R_new[p];  //adding two radius
                        if (dd[p] > R_n[p]) // non overlap checking. distance between spheres > the sum of radius 
                        {
                            count[p] = 0; //if condition is satisfied it stores the value as 0
                        }
                        else
                        {
                            count[p] = 1;   //if condition is not satisfied it stores the value as 1
                        }
                        p++;
                    }
                    int count1 = 0;
                    count1 = std::accumulate(count.begin(), count.end(), 0);

                    if (count1 > 0) // if the non overlapping condtion is not satisfied 
                    {


                        std::fill(count.begin(), count.end(), 0);

                        oLoop = oLoop + 1;
                        total = 0;     
                        n=n-1;
                    }
                    else if (count1 == 0) //if non overalpping condition is satisfied then accept the point
                    {
                        oLoop = 0; 
                        packedspheres[n] = packedspheres1[n];
                        total = total + 1;
                    }
                }
                                       
                if (oLoop >= 100000)
                {
                    Util::Abort(INFO, " spheres will not be generated with the given parameters. Try adjusting the V_per and R_per values ");
                }
                    
                Util::Message(INFO, "count  ", n);
                    
            }
            // Util::Message(INFO, "count  ", pack_size);
              

            std::ofstream file("/home/bkanagar/Desktop/alamo/tests/Flame/sphere_info.txt");
              
              

              
            file << "location"<< "  "<<"radius"<<std::endl;
            for(int i=0; i<pack_size;i++)
            {
                file << packedspheres[i]<< "  "<<R_new[i]<<std::endl;
                file << " " <<std::endl;
            }
              
            file.close();
             
                   
        }
        //Util::ParallelMessage(INFO);
        R_new.resize(pack_size);
        // Util::ParallelMessage(INFO);
        packedspheres.resize(pack_size);
        // Util::ParallelMessage(INFO);
        MPI_Bcast(R_new.data(), R_new.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_C_LONG_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        //MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);i
        //Util::ParallelMessage(INFO);
        MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_2DOUBLE_PRECISION, 0, MPI_COMM_WORLD);        
        //Util::ParallelMessage(INFO);    
       
    };

    // void Add(const int &lev, Set::FlameSpeed_mf<Set::Scalar>  &a_FlameSpeed)
    void Add(const int &lev, Set::Field<Set::Scalar> &a_phi)
    // void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
    {
        //Util::ParallelMessage(INFO);
        Set::Vector size;
        AMREX_D_TERM(size(0) = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];,
                     size(1) = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];,
                     size(2) = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2];)

            for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                // for (amrex::MFIter mfi(*a_FlameSpeed[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
        
                amrex::Box bx = mfi.tilebox();
                bx.grow(a_phi[lev]->nGrow());
                // bx.grow(a_FlameSpeed[lev]->nGrow());
                amrex::Array4<Set::Scalar> const &phi = a_phi[lev]->array(mfi);
                // amrex::Array4<Set::Scalar> const& FlameSpeed_mf= a_FlameSpeed[lev]->array(mfi);
                packedspheres.resize(pack_size);
                R_new.resize(pack_size);
                //Util::ParallelMessage(INFO);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                       {
                   
                                           Set::Vector x;
                                           Set::Scalar min_sph=0;
                                       
                                           AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                                        x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                                        x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

                                       

                                           for (int n = 0; n < pack_size; n++)
                                           {
                                               Set::Scalar d;
                                               Set::Scalar d1;

                                               d = (x - packedspheres[n]).lpNorm<2>();

                                               if (geom[0].isPeriodic(0))
                                               {
                                                   d1 = std::min((x - packedspheres[n] + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                                                                 (x - packedspheres[n] - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                                               }
#if AMREX_SPACEDIM > 1
                                               if (geom[0].isPeriodic(1))
                                               {
                                                   d1 = std::min((x - packedspheres[n] + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                                                                 (x - packedspheres[n] - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                                               }
#endif
#if AMREX_SPACEDIM > 2
                                               if (geom[0].isPeriodic(2))
                                               {
                                                   d1 = std::min((x - packedspheres[n] + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                                                                 (x - packedspheres[n] - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                                               }
#endif

                                               d = std::min(d, d1);

                                               if (d <= (tau*R_new[n]+ fact)) //the points d within the radius are assigned the 0 (HTPB) or 1 (AP)
                                           
                                               {                 
                                                   Set::Scalar m = 0.5*(1+erf((-d+R_new[n])/(eps))); //error function is used to avoid the sharp interface
                                                   min_sph = min_sph + m * (1. - min_sph);
                                               } 
                                            
                                               /*if (d< R_new[n])
                                                 {
                                                 min_sph=1;
                                                 }*/
                                            
                                           }

                                   
                                           if (type == Type::Values)
                                               phi(i, j, k) = min_sph;
                                           else if (type == Type::Partition) 
                                               phi(i, j, k, min_sph) = min_sph;
                                       
                                       });
            }
        //Util::ParallelMessage(INFO);
    }


    static void Parse(PackedSpheres &value, IO::ParmParse &pp)
    {
        pp.query("interval", value.interval);  // How many the radius between R_min and R_max. This is not used when R_mean used
        pp.query("R_min", value.R_min);  // Max radius-Introduced for varying the radius between R_min and R_max along with variable interval. This is not used when R_mean used
        pp.query("R_max", value.R_max); // Min radius- Introduced for varying the radius between R_min and R_max along with variable interval. This is not used when R_mean used
        pp.query("R_mean", value.R_mean); //Mean radius of sphere
        pp.query("vf", value.vf);  //volume fraction
        pp.query("verbosity", value.verbosity); // This factor is used when one radius is uded spheres for all spheres. If the no of spheres are not able to fit in the domain due to non overlapping condition checking, this will reduce the no of sphere and fit the spheres.
        pp.query("eps", value.eps);   //Interface thickness between AP and HTPB
        pp.query("fact", value.fact);     //This is used to control the interface thickness direction (towards HTPB or towards AP)
        pp.query("tau", value.tau);   // This controls the variation of the interface thickness in terms of the spheres radius in percentage
        pp.query("V_per1",value.V_per1);
        pp.query("V_per2",value.V_per2);
        pp.query("V_per3",value.V_per3);
        pp.query("V_per4",value.V_per4);
        pp.query("R_per1",value.R_per1);
        pp.query("R_per2",value.R_per2);
        pp.query("R_per3",value.R_per3);
        pp.query("R_per4",value.R_per4);    
        std::string type = "values";
        pp.query("type", type);
        if (type == "values")
            value.type = Type::Values;
        else if (type == "partition")
            value.type = Type::Partition;
        else
            Util::Abort(INFO);
        value.Define();
    }

private:
        
    int pack_size;
    std::vector<Set::Vector> packedspheres; //centre points of spheres that is accepted after checking non overlapping condition 
    std::vector<Set::Vector> packedspheres1; //temporary centre points of spheres produced randomly before checking non overlapping condition 
    std::vector<Set::Scalar> R_new; //
    std::vector<Set::Scalar> R;        //If R_min and R_max are used R is calculated using the formula. If R_mean is used then R=R_mean
    std::vector<Set::Scalar> dd; //distance between the spheres centres 
    std::vector<Set::Scalar> dd1; // distance between the sphere centers that includes periodic conditions
    std::vector<Set::Scalar> R_n; //adding radius of two spheres 
    std::vector<Set::Scalar> count; //if the spheres center is accepted after checking non overlapping condition , this variable is assigned 0. otherwise 1
    std::vector<double>::iterator it;
    Set::Scalar interval;
    Set::Scalar vf ; // Volume fraction
    Set::Scalar R_min; // Min radius of the spheres
    Set::Scalar R_max;  // Max radius of the spheres
    Set::Scalar R_mean;  // R_mean radius of the spheres
    Set::Scalar verbosity;   
    Set::Scalar eps;
    Set::Scalar tau;
    Set::Scalar fact;
    Set::Scalar V_per1,V_per2,V_per3,V_per4;
    Set::Scalar R_per1,R_per2,R_per3,R_per4;
    int rank;
    Type type = Type::Values;
    amrex::Vector<long double> packedspheres_x;
    amrex::Vector<long double> packedspheres_y;
#if BL_SPACEDIM == 3
    amrex::Vector<long double> packedspheres_z;
#endif
};
}
#endif
