#ifndef IC_PACKEDSPHERES_H_
#define IC_PACKEDSPHERES_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>

namespace IC
{
/// \brief Create a series of patches using the packedspheres Tesselation
///
/// For a thorough discussion of the packedspheres tesselation algorithm, see
/// [this wikipedia article](https://en.wikipedia.org/wiki/packedspheres_diagram)
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
class PackedSpheres : public IC
{
public:
	enum Type {Partition, Values};

	PackedSpheres (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
	//PackedSpheres (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains, amrex::Real _mean, amrex::Real _std_deviation,Set::Scalar a_alpha = 1.0) : IC(_geom)
	PackedSpheres (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains, Set::Scalar a_alpha = 1.0) : IC(_geom) 
 
	{
    	//Define(_number_of_grains,_mean,_std_deviation,a_alpha);
		Define(_number_of_grains,a_alpha);
  	}

  	//void Define (int a_number_of_grains,amrex::Real a_mean, amrex::Real a_std_deviation, Set::Scalar a_alpha = 1.0)
	void Define (int a_number_of_grains, Set::Scalar a_alpha = 1.0)  
	{
    	//Define(a_number_of_grains,a_mean, a_std_deviation, std::vector<Set::Scalar>(a_number_of_grains,a_mean, a_std_deviation,a_alpha),Type::Partition);
		Define(a_number_of_grains, std::vector<Set::Scalar>(a_number_of_grains,a_alpha),Type::Partition);
  	}

	void Define (int a_number_of_grains,
                 //amrex::Real a_mean,
                 //amrex::Real a_std_deviation,
				std::vector<Set::Scalar> a_alpha,
				Type a_type = Type::Values)
				//Type b_type = Type::Value2,
				//Type c_type = Type::Value3)
	{
		number_of_grains = a_number_of_grains;
       // mean= a_mean;
       // std_deviation=a_std_deviation;
		alpha = a_alpha;
		type = a_type;
		//type_b = b_type;
		//type_c = c_type;



	    std::vector <amrex::Real> dd(number_of_grains-1);
		std::vector <amrex::Real> R_n(number_of_grains-1);
		std::vector <amrex::Real> count(number_of_grains);
       

         
		
        int p=0;
		packedspheres.resize(number_of_grains);
		
		Set::Vector size;
		AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
					 size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
					 size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)	
					 int total=0;
		do
		{			 

		for (int n = 0; n<number_of_grains; n++)
		{
			AMREX_D_TERM(packedspheres[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random();,
						 packedspheres[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random();,
						 packedspheres[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(););
			R_new.resize(number_of_grains);	
			//amrex::Real R_max=0.15;
			//amrex::Real R_min=0.1;
			//amrex::Real R=R_max-R_min;
			//R_new[n]=(((((amrex::Real)rand())/((amrex::Real)RAND_MAX))*R+R_min)*100.0)/100.0;
           
			amrex::Real R_old=Util::Gaussian (mean=0.1,std_deviation=0);
            Util::Message(INFO,"R_old	    ", R_old);
            if (R_old>0)
	     {
	        R_new[n]= R_old;
	        Util::Message(INFO,"R_new	    ", R_new[n]);
	     }
		


			std::vector<bool> df(number_of_grains-1);
            std::fill(df.begin(), df.end(), true);	

	if (n>0)								
    {
 		p=0;

 while (p<n)
	 {		

                         dd[p] = (packedspheres[n] - packedspheres[p]).lpNorm<2>();

	   if (geom[0].isPeriodic(0))
						{
						 dd[p] = std::min( (packedspheres[n]-packedspheres[p] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
										   (packedspheres[n]-packedspheres[p] - size(0)*Set::Vector::Unit(0)).lpNorm<2>());
						}
		#if AMREX_SPACEDIM>1
					if (geom[0].isPeriodic(1))
						{
						 dd[p] = std::min( (packedspheres[n]-packedspheres[p] + size(1)*Set::Vector::Unit(1)).lpNorm<2>(),
										   (packedspheres[n]-packedspheres[p] - size(1)*Set::Vector::Unit(1)).lpNorm<2>());
						}
		#endif
		#if AMREX_SPACEDIM>2
					if (geom[0].isPeriodic(2))
						{
						 dd[p] = std::min( (packedspheres[n]-packedspheres[p] + size(2)*Set::Vector::Unit(2)).lpNorm<2>(),
										   (packedspheres[n]-packedspheres[p] - size(2)*Set::Vector::Unit(2)).lpNorm<2>());
						}
		#endif

       R_n[p]=R_new[n]+R_new[p];
           if (dd[p]>(R_n[p]))
			{
				df[p]=true;
			}
			else
			{
				df[p]=false;
			}          
	   
	   		p++;
	 }

	 		if (std::all_of(df.begin(),df.end(),[](bool i) {return i==true;}))
										
			{
												
				count[n]=1;
						
			}
		
	}
		total = accumulate(count.begin(), count.end(), 0, plus<int>());
		total=total+1;
		Util::Message(INFO,"count:	    ", total);
	    
    }
	    
	}
	while (total<number_of_grains);
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
							
				//amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
				//amrex::Real min_distance = 0.2;
				int min_grain_id = -1;

				for (int n = 0; n<number_of_grains; n++)
				{
					Set::Scalar d = (x - packedspheres[n]).lpNorm<2>();
					 

					if (geom[0].isPeriodic(0))
						{
							d = std::min( (x-packedspheres[n] + size(0)*Set::Vector::Unit(0)).lpNorm<2>(),
										  (x-packedspheres[n] - size(0)*Set::Vector::Unit(0)).lpNorm<2>());
						}
#if AMREX_SPACEDIM>1
					if (geom[0].isPeriodic(1))
						{
							d = std::min( (x-packedspheres[n] + size(1)*Set::Vector::Unit(1)).lpNorm<2>(),
										  (x-packedspheres[n] - size(1)*Set::Vector::Unit(1)).lpNorm<2>());
						}
#endif
#if AMREX_SPACEDIM>2
					if (geom[0].isPeriodic(2))
						{
							d = std::min( (x-packedspheres[n] + size(2)*Set::Vector::Unit(2)).lpNorm<2>(),
										  (x-packedspheres[n] - size(2)*Set::Vector::Unit(2)).lpNorm<2>());
						}
#endif
					if (d<R_new[n])
						{
							//min_distance = d;
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
	std::vector<Set::Vector> packedspheres;
	std::vector <amrex::Real> R_new;
    amrex::Real mean;
    amrex::Real std_deviation;
	Type type;
	amrex::Vector<amrex::Real> packedspheres_x;
	amrex::Vector<amrex::Real> packedspheres_y;
#if BL_SPACEDIM==3
	amrex::Vector<amrex::Real> packedspheres_z;
#endif
};
}
#endif