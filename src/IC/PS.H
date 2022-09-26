#ifndef IC_PS_H_
#define IC_PS_H_

#include "Set/Set.H"
#include "IC/IC.H"

namespace IC
{
class PS : public IC
{
public:
    enum Type {Partition, Values};

    PS (amrex::Vector<amrex::Geometry> &a_geom) : IC(a_geom) {}
    
    PS (amrex::Vector<amrex::Geometry> &a_geom, int a_nspheres, Set::Scalar a_matrix, Set::Scalar a_inclusion) 
    : IC(a_geom)
    {
        Define(a_nspheres,a_matrix,a_inclusion);
    }
    
    void Define(int a_nspheres, Set::Scalar a_matrix, Set::Scalar a_inclusion)
    {
        nspheres = a_nspheres;
        matrix = a_matrix;
        inclusion = a_inclusion;
        points.resize(nspheres);
        radii.resize(nspheres);
        for (int n = 0; n<nspheres; n++)
        {
            AMREX_D_TERM(points[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random();,
                        points[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random();,
                        points[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(););
            radii[n] = 0.25*Util::Random();
        }
    }

    //void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
    void Add(const int &lev,Set::Field<Set::Scalar> &a_field)
    {
        Set::Vector size;
        AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
                    size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
                    size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)
        amrex::IndexType type = a_field[lev]->ixType();

        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            //bx.grow(a_field[lev]->nGrow());
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                Set::Vector x = Set::Position(i,j,k,geom[lev],type);
                            
                bool inside = false;
                for (int n = 0; n<nspheres; n++)
                {
                    Set::Scalar d = (x - points[n]).lpNorm<2>();

                    if (geom[lev].isPeriodic(0))
                        {
                            d = std::min(d,std::min((x-points[n] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
                                                    (x-points[n] - size(0)*Set::Vector::Unit(0)).lpNorm<2>()));
                        }
#if AMREX_SPACEDIM>1
                    if (geom[lev].isPeriodic(1))
                        {
                            d = std::min(d,std::min((x-points[n] + size(1)*Set::Vector::Unit(1)).lpNorm<2>(),
                                                    (x-points[n] - size(1)*Set::Vector::Unit(1)).lpNorm<2>()));
                        }
#endif
#if AMREX_SPACEDIM>2
                    if (geom[lev].isPeriodic(2))
                        {
                            d = std::min(d,std::min((x-points[n] + size(2)*Set::Vector::Unit(2)).lpNorm<2>(),
                                                    (x-points[n] - size(2)*Set::Vector::Unit(2)).lpNorm<2>()));
                        }
#endif
                    if (d<radii[n]) inside = true;
                }

                if (inside) field(i,j,k) += inclusion;
                else field(i,j,k) += matrix;
            });
        }
    }
    
private:
    int nspheres;
    Set::Scalar matrix=0.0, inclusion=1.0;
    std::vector<Set::Scalar> radii;
    std::vector<Set::Vector> points;

public:
    static void Parse(PS & value, IO::ParmParse & pp)
    {
        int tmp_nspheres;
        Set::Scalar tmp_matrix = 0.0, tmp_inclusion = 1.0;
        pp.query("nspheres",tmp_nspheres);
        pp.query("matrix",tmp_matrix);
        pp.query("inclusion",tmp_inclusion);
        value.Define(tmp_nspheres,tmp_matrix,tmp_inclusion);
    }
};
}
#endif
