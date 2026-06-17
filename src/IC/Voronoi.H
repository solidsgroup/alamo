// Initialize an N component field with a Voronoi tessellation.
// (Mostly used by :code:`Integrator::PhaseFieldMicrostructure`).

#ifndef IC_VORONOI_H_
#define IC_VORONOI_H_

#include "AMReX_Arena.H"
#include "AMReX_GpuContainers.H"
#include "AMReX_GpuDevice.H"
#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"

namespace IC
{
class Voronoi : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "voronoi";

    enum Type {Partition, Values};

    Voronoi (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
    Voronoi (amrex::Vector<amrex::Geometry> &_geom,IO::ParmParse &pp, std::string name) : IC(_geom) 
    {pp_queryclass(name,*this);}
    Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains, Set::Scalar a_alpha) : IC(_geom) 
    {
        Define(_number_of_grains,a_alpha);
    }
    Voronoi (amrex::Vector<amrex::Geometry> &a_geom, int a_number_of_grains) : IC(a_geom) 
    {
        Define(a_number_of_grains,1.0);
    }

    void Define (int a_number_of_grains, Set::Scalar a_alpha)
    {
        Define(a_number_of_grains, std::vector<Set::Scalar>(a_number_of_grains,a_alpha),Type::Partition);
    }

    void Define (int a_number_of_grains,
                std::vector<Set::Scalar> a_alpha,
                Type a_type = Type::Values)
    {
        number_of_grains = a_number_of_grains;
        alpha = a_alpha;
        type = a_type;

        voronoi.resize(number_of_grains);

        for (int n = 0; n<number_of_grains; n++)
        {
            AMREX_D_TERM(voronoi[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random(true);,
                        voronoi[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random(true);,
                        voronoi[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(true););
        }

        alpha_device.resize(number_of_grains);
        amrex::Gpu::copyAsync(  amrex::Gpu::hostToDevice,
                                alpha.begin(), alpha.end(), alpha_device.begin());
        alpha_ptr = alpha_device.data();
        
#ifdef ALAMO_GPU
        voronoi_ptr = static_cast<Set::Vector*>(amrex::The_Device_Arena()->alloc(number_of_grains * sizeof(Set::Vector)));
        amrex::Gpu::htod_memcpy(voronoi_ptr, voronoi.data(), number_of_grains*sizeof(Set::Vector));
#else
        voronoi_ptr = voronoi.data();
#endif
    };
    
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {
        Set::Vector size;
        AMREX_D_TERM(   size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
                        size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
                        size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];);

        const Set::Vector problo(geom[lev].ProbLo());
        const Set::Vector cellsize(geom[lev].CellSize());
        const auto periodic = geom[lev].isPeriodic();
        const auto number_of_grains = this->number_of_grains;
        const auto type = this->type;
        const auto alpha_ptr = this->alpha_ptr;
        const auto voronoi_ptr = this->voronoi_ptr;

        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType ixtype = a_field[lev]->ixType();
            if (ixtype == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (ixtype == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();

            int ncomp = a_field[lev]->nComp();
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                Set::Vector x;
                AMREX_D_TERM(   x(0) = problo(0) + ((Set::Scalar)(i) + 0.5) * cellsize(0);,
                                x(1) = problo(1) + ((Set::Scalar)(j) + 0.5) * cellsize(1);,
                                x(2) = problo(2) + ((Set::Scalar)(k) + 0.5) * cellsize(2););
                            
                Set::Scalar min_distance = 1E100; //std::numeric_limits<amrex::Real>::infinity();
                int min_grain_id = -1;

                for (int n = 0; n<number_of_grains; n++)
                {
                    Set::Scalar d = (x - voronoi_ptr[n]).lpNorm<2>();

                    if (periodic[0])
                        {
                            d = std::min(   d,
                                            std::min(   (x-voronoi_ptr[n] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
                                                        (x-voronoi_ptr[n] - size(0)*Set::Vector::Unit(0)).lpNorm<2>()));
                        }
#if AMREX_SPACEDIM>1
                    if (periodic[1])
                        {
                            d = std::min(d,
                                        std::min( (x-voronoi_ptr[n] + size(0)*Set::Vector::Unit(1)).lpNorm<2>(),
                                                    (x-voronoi_ptr[n] - size(0)*Set::Vector::Unit(1)).lpNorm<2>()));
                        }
#endif
#if AMREX_SPACEDIM>2
                    if (periodic[2])
                        {
                            d = std::min(d,
                                        std::min( (x-voronoi_ptr[n] + size(0)*Set::Vector::Unit(2)).lpNorm<2>(),
                                                    (x-voronoi_ptr[n] - size(0)*Set::Vector::Unit(2)).lpNorm<2>()));
                        }
#endif
                    if (d<min_distance)
                        {
                            min_distance = d;
                            min_grain_id = n;
                        }
                }

                if (type == Type::Values) field(i,j,k) = alpha_ptr[min_grain_id];
                else if (type == Type::Partition) field(i,j,k,min_grain_id % ncomp) = alpha_ptr[min_grain_id];
            });
        }
    }

    static void Parse(Voronoi &value, IO::ParmParse &pp)
    {
        pp_query("number_of_grains",value.number_of_grains); // Number of grains
        if (pp.contains("alpha"))
        {
            pp_queryarr("alpha",value.alpha); // Value to take in the region [1.0]
            value.Define(value.number_of_grains,value.alpha);    
        }
        else
        {
            value.Define(value.number_of_grains,1.0);
        }
        pp_query("seed",value.seed); // Random seed to use
    }
    
private:
    int number_of_grains;
    int seed = 1;
    std::vector<Set::Scalar> alpha;
    std::vector<Set::Vector> voronoi;

    amrex::Gpu::DeviceVector<Set::Scalar> alpha_device;
    
    Set::Scalar *alpha_ptr;
    Set::Vector* voronoi_ptr;

    Type type;
};
}
#endif
