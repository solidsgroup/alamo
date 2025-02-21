//
// This simulates the creation of "disconnection pairs" at a GB
// by perturbing the order parameter :math:`\eta` with a gaussian.
//
// Disconnections can be nucleated randomly, or they can be added
// in a controlled way through the "fixed" option.
//
// See the following reference for further details
// 
// .. bibliography::
//    :list: none
//    :filter: False
// 
//    gokuli2021multiphase
// 
//

#ifndef MODEL_DEFECT_DISCONNECTION
#define MODEL_DEFECT_DISCONNECTION

#include <random>

#include "AMReX_SPACE.H"
#include "IO/ParmParse.H"
#include "Set/Base.H"
#include "Util/MPI.H"

namespace Model
{
namespace Defect
{
class Disconnection
{
public:

    Disconnection(){};
    ~Disconnection() {};

    static void
    Parse(Disconnection &value, IO::ParmParse &pp)
    {
        Util::Assert(INFO,TEST(AMREX_SPACEDIM==2),"2D only");
        // time to start applying disconnections
        pp_query_default ("tstart", value.tstart,0.0); 
        // nucleation energy
        pp_query_default ("nucleation_energy", value.nucleation_energy,0.0); 
        // characteristic time
        pp_query_default("tau_vol", value.tau_vol, 1.0);
        // temperature
        pp_query_default("temp", value.temp, 0.0);
        // characteristic size
        pp_query_default("box_size", value.box_size, 0.0);
        // interval between generation events
        pp_query_required("interval", value.interval);
        // regularization epsilon
        pp_query_default ("epsilon",value.epsilon,1E-20); 

        // whether to manually specify disconnection nucleation points
        pp.query("disconnection.fixed.on",value.fixed.on); 
        if (value.fixed.on)
        {
            // array of x locations
            pp.queryarr("fixed.sitex",value.fixed.sitex); 
            // array of y locations
            pp.queryarr("fixed.sitey",value.fixed.sitey); 
            // array of order parameter number
            pp.queryarr("fixed.phases",value.fixed.phases); 
            // time to appear
            pp.queryarr("fixed.time",value.fixed.time); 
            Util::Assert(INFO,TEST(value.fixed.sitex.size() == value.fixed.sitey.size()));
            Util::Assert(INFO,TEST(value.fixed.sitex.size() == value.fixed.phases.size()));
            Util::Assert(INFO,TEST(value.fixed.sitex.size() == value.fixed.time.size()));
            value.fixed.done.resize(value.fixed.sitex.size(),false);
        }
        else
        {
            value.unif_dist = std::uniform_real_distribution<double>(0.0,1.0);
            value.int_dist = std::uniform_int_distribution<int>(0,1);
            value.rand_num_gen.seed(amrex::ParallelDescriptor::MyProc());
        }

        // verbosity
        pp_query_default("verbose",value.verbose,false); 
    }


    /// This operates on an entire field, and manages all of the MPI
    /// communication necessary for consistent nucleation.
    void
    Nucleate(   Set::Field<Set::Scalar> &eta_mf,
                std::vector<amrex::Geometry> &geom,
                Set::Scalar timestep,
                Set::Scalar time,
                int iter
        )
    {
        Util::Assert(INFO,TEST(eta_mf[0]->nComp() == 2), "This only works for 2 component phase fields");

        if (time < tstart) // wait until it's time to go
            return;
        if (iter % interval) // skip every [interval] timesteps
            return;

        sitex.clear();
        sitey.clear();
        phases.clear();

        int max_lev = eta_mf.finest_level;

        const Set::Scalar *DX = geom[max_lev].CellSize();
        Set::Scalar exponent = DX[0] * DX[0] * (timestep * interval) / tau_vol;

        // Determine the nucleation sites in the finest grid only
        for (amrex::MFIter mfi(*eta_mf[max_lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<amrex::Real> const &eta = (*eta_mf[max_lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar E0 = 2.0*nucleation_energy;
                E0 /= epsilon + 256.0*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1);
                Set::Scalar p = std::exp(-E0/(K_b*temp));
                Set::Scalar P = 1.0 - std::pow(1.0 - p,exponent);
                if (eta(i,j,k,0) < 0 || eta(i,j,k,0) > 1.0 || eta(i,j,k,1) < 0 || eta(i,j,k,1) > 1.0) P = 0.0;
                Set::Scalar q = 0.0;
                q = unif_dist(rand_num_gen);
                if (q < P)
                {
                    sitex.push_back(geom[max_lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0]);
                    sitey.push_back(geom[max_lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1]);
                    int phase = int_dist(rand_num_gen);
                    phases.push_back(phase);
                } });
        }

        // Sync up all the nucleation sites among processors
        Util::MPI::Allgather(sitex);
        Util::MPI::Allgather(sitey);
        Util::MPI::Allgather(phases);

        if (verbose)
        {
            Util::Message(INFO, "Nucleating ", phases.size(), " disconnections");
        }

        if (sitex.size() > 0)
        {
            // Now that we all know the nucleation locations, perform the nucleation
            for (int lev = 0; lev <= max_lev; lev++)
            {
                amrex::Box domain = geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                const amrex::Real *DX = geom[lev].CellSize();
                for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box bx = mfi.grownnodaltilebox() & domain;
                    amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Vector x;
                        AMREX_D_TERM(
                            x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0];,
                            x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1];,
                            x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * DX[2];);
                        for (unsigned int m = 0; m < phases.size(); m++)
                        {
                            amrex::Real r_squared = 0;
                            Set::Vector nucleation_site(AMREX_D_DECL(sitex[m], sitey[m], 0.0));
                            for (int n = 0; n < AMREX_SPACEDIM; n++)
                            {
                                amrex::Real dist = nucleation_site(n) - x(n);
                                r_squared += dist * dist;
                            }
                            amrex::Real bump = exp(-r_squared / box_size);
                            eta(i, j, k, phases[m]) = bump * (1 - eta(i, j, k, phases[m])) + eta(i, j, k, phases[m]);
                            eta(i, j, k, 1 - phases[m]) = (1. - bump) * eta(i, j, k, 1 - phases[m]);
                        }
                    });
                }
            }
        }
    }
    


private:

    bool verbose = false;

    Set::Scalar tstart = NAN;
    Set::Scalar nucleation_energy = NAN;
    Set::Scalar tau_vol = NAN;
    Set::Scalar temp = NAN;
    Set::Scalar box_size = NAN;
    Set::Scalar epsilon = NAN;

    int interval = -1;

    std::uniform_real_distribution<double> unif_dist;  /// random number distribution for spatial location
    std::uniform_int_distribution<int> int_dist;       /// random number generator for phase
    std::default_random_engine rand_num_gen;           /// generator object

    struct {
        int on = 0;
        std::vector<Set::Scalar> sitex;
        std::vector<Set::Scalar> sitey;
        std::vector<int> phases;
        std::vector<Set::Scalar> time;
        std::vector<bool> done;
    } fixed;


    std::vector<Set::Scalar> sitex; /// list of nucleation site x coordinates
    std::vector<Set::Scalar> sitey; /// list of nucleation stie y coordinates
    std::vector<int> phases;        /// list of nucleation site phases (up or down)

    const Set::Scalar K_b = 8.617333262145e-5; // eV/K
};
}
}

#endif
