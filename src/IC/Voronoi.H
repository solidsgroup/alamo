#ifndef IC_VORONOI_H_
#define IC_VORONOI_H_

#include "Set/Set.H"
#include "IC/IC.H"

namespace IC
{
/// \brief Create a series of patches using the Voronoi Tesselation
///
/// For a thorough discussion of the Voronoi tesselation algorithm, see
/// [this wikipedia article](https://en.wikipedia.org/wiki/Voronoi_diagram)
///
/// Let `N` be the number of patches in the tesselation, 
/// and let \f$\alpha_1,\ldots,\alpha_n\f$ be a set of scalar values.
/// There are two methods for applying the tesselation IC:
///
///   * `Partition`: This requires an `N`-component multifab. 
///     (There is currently no check for this.)
///     Each component of the fab corresponds to one of the patches
///     Then fab component i is equal to \f$\alpha_i\f$ within the patch and
///     zero elsewhere. (This is mostly useful for multiphase field 
///     microstructure)
///   * `Values`: This uses a 1-component fab. It sets the value of the ith fab 
///     equal to \f$\alpha_i\f$.
///
/// Note on periodicity: This IC uses the `amrex::Geometry` object to determine
/// if it is periodic, and if so, it will generate a periodic tesselation.
///
class Voronoi : public IC
{
public:
	enum Type {Partition, Values};

	Voronoi (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
	Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains, Set::Scalar a_alpha = 1.0) : IC(_geom) 
	{
    	Define(_number_of_grains,a_alpha);
  	}

  	void Define (int a_number_of_grains, Set::Scalar a_alpha = 1.0)
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
			AMREX_D_TERM(voronoi[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random();,
						 voronoi[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random();,
						 voronoi[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(););
		}
	};
	
	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
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
	
private:
	int number_of_grains;
	std::vector<Set::Scalar> alpha;
	std::vector<Set::Vector> voronoi;
	Type type;
	amrex::Vector<amrex::Real> voronoi_x;
	amrex::Vector<amrex::Real> voronoi_y;
#if BL_SPACEDIM==3
	amrex::Vector<amrex::Real> voronoi_z;
#endif
};
}
#endif
