#include <AMReX_ParallelDescriptor.H>
#include "Voronoi.H"

namespace IC
{

  
void Voronoi::Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
{
	Set::Vector size;
	AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
		     	 size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
		     	 size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)

	for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		bx.grow(a_field[lev]->nGrow());
		amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

			Set::Vector x;
			AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
						 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
						 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
						 
			amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
			int min_grain_id = -1;

			for (int n = 0; n<number_of_grains; n++)
			{
				Set::Scalar d = (x - voronoi[n]).lpNorm<2>();

				if (geom[0].isPeriodic(0))
					{
						d = std::min( (x-voronoi[n] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
									  (x-voronoi[n] - size(0)*Set::Vector::Unit(0)).lpNorm<2>());
					}
#if AMREX_SPACEDIM>1
				if (geom[0].isPeriodic(1))
					{
						d = std::min( (x-voronoi[n] + size(1)*Set::Vector::Unit(1)).lpNorm<2>(),
									  (x-voronoi[n] - size(1)*Set::Vector::Unit(1)).lpNorm<2>());
					}
#endif
#if AMREX_SPACEDIM>2
				if (geom[0].isPeriodic(2))
					{
						d = std::min( (x-voronoi[n] + size(2)*Set::Vector::Unit(2)).lpNorm<2>(),
									  (x-voronoi[n] - size(2)*Set::Vector::Unit(2)).lpNorm<2>());
					}
#endif
				if (d<min_distance)
					{
						min_distance = d;
						min_grain_id = n;
					}
			}

			if (type == Type::Values) field(i,j,k) = alpha[min_grain_id];
			else if (type == Type::Partition) field(i,j,k,min_grain_id) = alpha[min_grain_id];
		});
	}
}
}
