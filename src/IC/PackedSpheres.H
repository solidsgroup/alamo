#ifndef IC_PACKEDSPHERES_H_
#define IC_PACKEDSPHERES_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
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
            std::vector<amrex::Real> sph_nos;
            Set::Scalar sph_tot_vol = 0;
            R.resize(number_of_grains);
            R_n.resize(number_of_grains);
            sph_vol.resize(number_of_grains);
            sph_nos.resize(number_of_grains);
            R_new.resize(number_of_grains);
            pack_size = 0;
            sph_tot_vol = 0;
            int inc = 0;
            int jj = 0;
            std::vector<double> ins(1, 1);
            packedspheres.resize(number_of_grains);

            int p = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            Set::Vector size;
            AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
                         size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,       
                         size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)
            int total = 0;
  
#if AMREX_SPACEDIM > 2
            if (geom[0].isPeriodic(2))
            {
                double dom_vol = size(0) * size(1) * size(2);
                // converting volume fraction in to no of spheres
                for (int n = 0; n < number_of_grains; n++)
                {
                    Util::Message(INFO, "3 dimensions");
                    if (n == 0)
                    {
                        // R[n]=R_min+((R_max-R_min)/number_of_grains);
                        R[n] = ((R_mean) / number_of_grains);
                        sph_vol[n] = (1.33) * 3.14 * verbosity * (R[n] * R[n] * R[n]);
                        // sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*1);
                        sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                        sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                        pack_size = pack_size + sph_nos[n];
                    }

                    else if (n > 0)
                    {
                        R[n] = ((R_mean) / number_of_grains);
                        sph_vol[n] = (1.33) * 3.14 * verbosity * (R[n] * R[n] * R[n]);
                        // sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*1);
                        sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                        sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                        pack_size = pack_size + sph_nos[n];
                    }
                }
            }
//#endif

Util::Message(INFO, "R_mean ", R_mean);

#elif AMREX_SPACEDIM > 1
            if (geom[0].isPeriodic(1))
            {
                double dom_vol = size(0) * size(1);
                Util::Message(INFO, "dom_vol ",dom_vol);
                // converting volume fraction in to no of spheres
                for (int n = 0; n < number_of_grains; n++)
                {
                    Util::Message(INFO, "2 dimensions");
                    if (n == 0)
                    {
                        // R[n]=R_min+((R_max-R_min)/number_of_grains);
                        R[n] = ((R_mean) / number_of_grains);
                        sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                        // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                        // sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*1);
                        sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                        //Util::Message(INFO, "sph_nos ",sph_nos[n]);
                        sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                        pack_size = pack_size + sph_nos[n];
                    }

                    else if (n > 0)
                    {
                        R[n] = ((R_mean) / number_of_grains);
                        sph_vol[n] = 3.14 * verbosity * (R[n] * R[n]);
                       // Util::Message(INFO, "sph_vol ",sph_vol[n]);
                        // sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*1);
                        sph_nos[n] = ceil((dom_vol / (sph_vol[n] * 1)) * vf);
                        // Util::Message(INFO, "sph_nos ",sph_nos[n]);
                        sph_tot_vol = sph_tot_vol + sph_vol[n] * sph_nos[n];
                        pack_size = pack_size + sph_nos[n];
                    }
                }
            }
#endif
            
            if (rank == 0)
            {
                R_new.assign(R.begin(), R.end());
                
                // matching the total no of elements with the total no spheres
                for (int i = 0; i < number_of_grains; i++)
                {

                    if (sph_nos[i] > 1)
                    {

                        inc = sph_nos[i] - 1;

                        if (i == 0)
                        {
                            it = R_new.begin();
                            ins.clear();
                            std::vector<double> ins(inc, R[i]);
                            R_new.insert(it + 1, ins.begin(), ins.end());
                            jj = jj + inc + 1;
                        }

                        else
                        {

                            it = R_new.begin();
                            ins.clear();
                            std::vector<double> ins(inc, R[i]);
                            R_new.insert(it + jj, ins.begin(), ins.end());
                            jj = jj + inc + 1;
                        }
                    }
                }

                std::reverse(R_new.begin(), R_new.end());

                int oLoop = 0;

         label1:

                dd.resize(pack_size - 1);
                dd1.resize(pack_size - 1);
                R_n.resize(pack_size - 1);
                packedspheres.resize(pack_size);
                packedspheres1.resize(pack_size);
                R_new.resize(pack_size);
                count.resize(pack_size);
                int iLoop = 0;

                for (int n = 0; n < pack_size; n++)
                {

                    std::fill(count.begin(), count.end(), 0);

                    if (n == 0)
                    {
                        AMREX_D_TERM(packedspheres1[0](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0) - geom[0].ProbLo(0)) * Util::Random();,
                                     packedspheres1[0](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1) - geom[0].ProbLo(1)) * Util::Random();,
                                     packedspheres1[0](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2) - geom[0].ProbLo(2)) * Util::Random(););
                        packedspheres[n] = packedspheres1[n];
                    }

                    else if (n > 0)
                    {

                    label:

                        AMREX_D_TERM(packedspheres1[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0) - geom[0].ProbLo(0)) * Util::Random();,
                                     packedspheres1[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1) - geom[0].ProbLo(1)) * Util::Random();,
                                     packedspheres1[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2) - geom[0].ProbLo(2)) * Util::Random(););

                        p = 0;

                        while (p < n)
                        {

                            dd[p] = (packedspheres1[n] - packedspheres1[p]).lpNorm<2>();

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

                            R_n[p] = R_new[n] + R_new[p];
                            if (dd[p] > R_n[p])
                            {
                                count[p] = 0;
                            }
                            else
                            {
                                count[p] = 1;
                            }
                            p++;
                        }
                        int count1 = 0;
                        count1 = std::accumulate(count.begin(), count.end(), 0);

                        if (count1 > 0)
                        {

                            if (iLoop < 100000)
                            {

                                iLoop = iLoop + 1;
                                std::fill(count.begin(), count.end(), 0);
                                goto label;
                            }
                            else
                            {
                                oLoop = oLoop + 1;
                                total = 0;
                                double dom_vol = size(0) * size(1);

                                double sph_vol1 = 3.14 * verbosity * (R_mean * R_mean);
                                double sph_nos1 = ceil(((dom_vol * 0.6) / (sph_vol1 * 1)) * vf);
                                //double sph_tot_vol1 = sph_vol1 * sph_nos1;
                                double pack_size1 = sph_nos1;
                               // Util::Message(INFO, "sph_tot_vol1 ", sph_tot_vol1);

                                

                                double sph_vol2 = 3.14 * verbosity * (R_mean * R_mean * 0.5 * 0.5);
                                double sph_nos2 = ceil(((dom_vol * 0.2) / (sph_vol2 * 1)) * vf);
                               // double sph_tot_vol2 = sph_vol2 * sph_nos2;
                                double pack_size2 = pack_size1 + sph_nos2;
                               //Util::Message(INFO, "sph_tot_vol2 ", sph_tot_vol2);

                                

                                double sph_vol3 = 3.14 * verbosity * (R_mean * R_mean * 0.3 * 0.3);
                                double sph_nos3 = ceil(((dom_vol * 0.1) / (sph_vol3 * 1)) * vf);
                               // double sph_tot_vol3 = sph_vol3 * sph_nos3;
                                double pack_size3 = pack_size2 + sph_nos3;
                                //Util::Message(INFO, "sph_tot_vol3 ", sph_tot_vol3);

                                double sph_vol4 = 3.14 * verbosity * (R_mean * R_mean * 0.2 * 0.2);
                                double sph_nos4 = ceil(((dom_vol * 0.1) / (sph_vol4 * 1)) * vf);
                                //double sph_tot_vol4 = sph_vol4 * sph_nos4;
                                double pack_size4 = pack_size3 + sph_nos4;
                                //Util::Message(INFO, "sph_tot_vol4 ", sph_tot_vol4);

                                
                                pack_size = pack_size4;
                                R_new.resize(pack_size);

                                for (int i = 0; i < (pack_size1); i++)
                                {
                                    R_new[i] = R_mean;
                                }

                                for (int i = pack_size1; i < (pack_size2); i++)
                                {
                                    R_new[i] = R_mean * 0.5;
                                }

                                for (int i = pack_size2; i < (pack_size3); i++)
                                {
                                    R_new[i] = R_mean * 0.3;
                                }

                                for (int i = pack_size3; i < (pack_size); i++)
                                {
                                    R_new[i] = R_mean * 0.2;
                                }
                                goto label1;
                            }
                        }
                        else if (count1 == 0)
                        {
                            iLoop = 0;
                            packedspheres[n] = packedspheres1[n];
                            total = total + 1;
                        }
                    }

                    if (oLoop > 25)
                    {
                        Util::Abort(INFO, " spheres will not be generated with the given parameters. Try adjusting the verbosity to produce spheres ");
                    }

                    Util::Message(INFO, "count  ", total);
                }
            }

            R_new.resize(pack_size);
            packedspheres.resize(pack_size);
            MPI_Bcast(R_new.data(), R_new.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            // MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_C_LONG_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
        };

        // void Add(const int &lev, Set::FlameSpeed_mf<Set::Scalar>  &a_FlameSpeed)
        void Add(const int &lev, Set::Field<Set::Scalar> &a_phi)
        // void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
        {

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
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {
                                       Set::Vector x;
                                       double min_grain_id=0;
                                       
                                       AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                                    x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                                    x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

                                       packedspheres.resize(pack_size);
                                       R_new.resize(pack_size);

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

                                           R_new.resize(pack_size);
                                           
                                           
                                           
                                          /* if (d < R_new[n])
                                           {
                                               min_grain_id = 1;
                                           }*/

                                           /* if (x(0)>-1.5 && x(0)<1.5 )
                                               {
                                                   min_grain_id = 1;
                                               }*/

                                         /*  if (x(1)>-0.05 && x(1)<0.05 )
                                               {
                                                   min_grain_id = 0;
                                               }*/
                                               
                                           //if (x(1)>=-0.05 && x(1)<=0.05 )
                                          // {    

                                           //if (d <= R_new[n])
                                            {
												//min_grain_id = (1/2)*(1+erf(x(0)/sqrt(2)));
												//min_grain_id=(erf((x(0)*40)));
												//min_grain_id=(1/2)*(1+erf(x(0)/sqrt(2)));
												//min_grain_id=(1/2)*(1+erf((x(0))/sqrt(2)));
												
												//min_grain_id=(erf((-d)/0.05));
												//min_grain_id=(erf((x(0)-1)/0.05));
												min_grain_id=(erf((x(0)-1)/0.05))-(erf((x(0))/0.05));
												
												//min_grain_id=(1/2)*(0.5+erf((x(0))/0.05));
												//a=a+0.001;
												 //Util::Message(INFO, "min_grain_id   ", min_grain_id );
												
											} 
										//}
                                               
                                       }

                                       /*if (type == Type::Values) FlameSpeed(i, j, k) = alpha[min_grain_id];
                                       else if (type == Type::Partition) FlameSpeed(i, j, k,min_grain_id) = alpha[min_grain_id];*/

                                       if (type == Type::Values)
                                           phi(i, j, k) = min_grain_id;
                                       else if (type == Type::Partition)
                                           phi(i, j, k, min_grain_id) = min_grain_id;
                                   });
            }
        }


        static void Parse(PackedSpheres &value, IO::ParmParse &pp)
        {
            pp.query("number_of_grains", value.number_of_grains);
            pp.query("R_min", value.R_min);
            pp.query("R_max", value.R_max);
            pp.query("R_mean", value.R_mean);
            pp.query("vf", value.vf);
            pp.query("verbosity", value.verbosity);
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
        std::vector<Set::Vector> packedspheres;
        std::vector<Set::Vector> packedspheres1;
        std::vector<Set::Scalar> R_new;
        std::vector<Set::Scalar> R;
        std::vector<Set::Scalar> dd;
        std::vector<Set::Scalar> dd1;
        std::vector<Set::Scalar> R_n;
        std::vector<Set::Scalar> count;
        std::vector<double>::iterator it;
        Set::Scalar number_of_grains;
        Set::Scalar vf ; // Volume fraction
        Set::Scalar R_min;
        Set::Scalar R_max;
        Set::Scalar R_mean;
        Set::Scalar verbosity;
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
