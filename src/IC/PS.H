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

	PS (amrex::Vector<amrex::Geometry> &a_geom, int a_nspheres, Set::Scalar a_matrix, Set::Scalar a_inclusion) 
    : IC(a_geom), nspheres(a_nspheres), matrix(a_matrix), inclusion(a_inclusion)
	{
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

	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
	{
		Set::Vector size;
		AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
					 size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
					 size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			//bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

				Set::Vector x;
				AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
							 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
							 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
							
                bool inside = false;
				for (int n = 0; n<nspheres; n++)
				{
					Set::Scalar d = (x - points[n]).lpNorm<2>();

					if (geom[0].isPeriodic(0))
						{
							d = std::min(d,std::min((x-points[n] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
										            (x-points[n] - size(0)*Set::Vector::Unit(0)).lpNorm<2>()));
						}
#if AMREX_SPACEDIM>1
					if (geom[0].isPeriodic(1))
						{
							d = std::min(d,std::min((x-points[n] + size(1)*Set::Vector::Unit(1)).lpNorm<2>(),
										            (x-points[n] - size(1)*Set::Vector::Unit(1)).lpNorm<2>()));
						}
#endif
#if AMREX_SPACEDIM>2
					if (geom[0].isPeriodic(2))
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
};
}
#endif
