#ifndef IC_PACKEDSPHERES_H_
#define IC_PACKEDSPHERES_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>


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
	enum Type {Partition, Values};

	PackedSpheres (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
	PackedSpheres (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains,Set::Scalar a_alpha = 1.0,amrex::Real _volume_fraction =0.05, amrex::Real _R_min=0.1,amrex::Real _R_max=0.2) : IC(_geom) 
	{
    	Define(_number_of_grains,a_alpha,_volume_fraction,_R_min,_R_max);
  	}

  	void Define (int a_number_of_grains, Set::Scalar a_alpha = 1.0,amrex::Real a_volume_fraction=0.05, amrex::Real a_R_min =0.1,amrex::Real a_R_max =0.2)
	{
    	Define(a_number_of_grains, std::vector<Set::Scalar>(a_number_of_grains,a_alpha),a_volume_fraction, a_R_min,a_R_max,Type::Partition);
  	}

	void Define (int a_number_of_grains,
				 std::vector<Set::Scalar> a_alpha,
				 amrex::Real a_volume_fraction,
                 amrex::Real a_R_min,
				 amrex::Real a_R_max,
				 Type a_type = Type::Values)
	{
		number_of_grains = a_number_of_grains;
		alpha = a_alpha;
		vf= a_volume_fraction;
        R_min=a_R_min;
		R_max=a_R_max;
		type = a_type;

		
		
    
		std::vector <amrex::Real> sph_vol;   
		std::vector <amrex::Real> sph_nos; 
        amrex::Real sph_tot_vol=0;
		R.resize(number_of_grains); 
		R_n.resize(number_of_grains); 
        sph_vol.resize(number_of_grains);
		sph_nos.resize(number_of_grains);
		count.resize(number_of_grains);
		R_new.resize(number_of_grains);
		pack_size=0;
		sph_tot_vol=0;
		int inc=0;
	    int jj=0;
		std::vector<double> ins (1,1);	
        packedspheres.resize(number_of_grains);


		int p=0;
	    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	   
		Set::Vector size;
		AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
					 size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
					 size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)	
					 int total=0;
					
        amrex::Real dom_vol=size(0)*size(1)*size(2);
        // converting volume fraction in to no of spheres
        for (int n = 0; n<number_of_grains; n++)
		{
		if (n==0)
		{		
        R[n]=(R_max-R_min)/number_of_grains;
		sph_vol[n]=(4/3)*3.14*(R[n]*R[n]*R[n]);
		sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*vf);
        sph_tot_vol=sph_tot_vol+sph_vol[n]*sph_nos[n];
		pack_size=pack_size+sph_nos[n];
		}

		else if (n>0)
		{	
        R[n]=R[n-1]+((R_max-R_min)/number_of_grains);
		sph_vol[n]=(4/3)*3.14*(R[n]*R[n]*R[n]);
		sph_nos[n]=ceil(ceil(dom_vol/number_of_grains/(sph_vol[n]*1))*vf);
        sph_tot_vol=sph_tot_vol+sph_vol[n]*sph_nos[n];
		pack_size=pack_size+sph_nos[n];
		
		}
		}

//        std::cout<<"sph_tot_vol "<<sph_tot_vol<<endl;
//		std::cout<<"vol frac "<<sph_tot_vol/(size(0)*size(1)*size(2))<<endl;
		
      
	  	
		
  
   
    
		
    packedspheres.resize(pack_size);
	R_new.resize(pack_size);
 
     if (rank==0)			 
		{  
	R_new.assign(R.begin(),R.end());
	
	
// matching the total no of elements with the total no spheres
for (int i=0; i<number_of_grains;i++)
{
	
	if (sph_nos[i]>1)
	{
	   
	   inc=sph_nos[i]-1;	
	 
		
			  if (i==0)
			 {
               it = R_new.begin();	
			   ins.clear();
		       std::vector<double> ins (inc,R[i]);
		       R_new.insert(it+1,ins.begin(),ins.end());
			    jj=jj+inc+1;	
			   
			 }

			else
			{

				it = R_new.begin();
				ins.clear();
		        std::vector<double> ins (inc,R[i]);
		        R_new.insert(it+jj,ins.begin(),ins.end());
				jj=jj+inc+1;	
				
			}
			
		 
	   
	   
	}
	
	else
	{
		
          
	}

	       
}

        std::reverse(R_new.begin(),R_new.end());
        dd.resize(pack_size-1);
		R_n.resize(pack_size-1);
		count.resize(pack_size);
		std::fill(count.begin(), count.end(), 0);



      

	do
		{			 

		for (int n = 0; n<pack_size; n++)
		{
			std::vector<bool> df(pack_size-1);
            std::fill(df.begin(), df.end(), true);

	if(n==0)
	{
            AMREX_D_TERM(packedspheres[0](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random();,
						 packedspheres[0](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random();,
						 packedspheres[0](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(););

	}

	 	
	else if (n>0)								
    {
		    
	int i=0;
while (i<1000)
{


		AMREX_D_TERM(packedspheres[n](0) = geom[0].ProbLo(0) + (geom[0].ProbHi(0)-geom[0].ProbLo(0))*Util::Random();,
						 packedspheres[n](1) = geom[0].ProbLo(1) + (geom[0].ProbHi(1)-geom[0].ProbLo(1))*Util::Random();,
						 packedspheres[n](2) = geom[0].ProbLo(2) + (geom[0].ProbHi(2)-geom[0].ProbLo(2))*Util::Random(););
						 

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
           if (dd[p]>((R_n[p])+0.001))
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
				i=1000;
						
			}
			else if (std::any_of(df.begin(),df.end(),[](bool i) {return i==false;}))
			{
				count[n]=0;
				i++;
			    
			}
       }
	}


	
		total = accumulate(count.begin(), count.end(), 0, plus<int>());
		total=total+1;
		std::cout<<"count  "<<total<<endl;
	    
    }

	}
	while (total<pack_size);
    }
	R_new.resize(pack_size);
	packedspheres.resize(pack_size);
	MPI_Bcast(R_new.data(), R_new.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(packedspheres.data(), packedspheres.size(), MPI_C_LONG_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	
	};
	
	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
	{

 
        Set::Vector size;  
        AMREX_D_TERM(size(0) = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];,
					 size(1) = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];,
					 size(2) = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2];)
		
		
			
		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
                 Set::Vector x;
				 int min_grain_id = -1;
				AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) *geom[lev].CellSize()[0];,
							x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) *geom[lev].CellSize()[1];,
							x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
			    
                packedspheres.resize(pack_size);
				R_new.resize(pack_size);
				
				for (int n = 0; n<pack_size; n++)
				{ 
			    Set::Scalar d ;
				Set::Scalar d1; 
							
				
              // dividing the domain in to 8 x direction and 4 in y direction and none in z direction  
				
				for(int v=-2;v<2;v++)	
				{
				 for (int u=-4;u<4;u++) 
				 {

				 
					if ((packedspheres[n](0)> (size(0)*(u*0.125))) && (packedspheres[n](0) < (size(0)*((u+1)*0.125))) && (packedspheres[n](1)> (size(1)*(v*0.25))) && (packedspheres[n](1)<  (size(1)*((v+1)*0.25))))
					{

					 d = (x - packedspheres[n]).lpNorm<2>();
					 d1 = 10;                
					if (geom[0].isPeriodic(0))
						{
						    d1 = std::min( (x-packedspheres[n] + size(0)*0.25*Set::Vector::Unit(0)).lpNorm<2>(),
										  (x-packedspheres[n] - size(0)*0.25*Set::Vector::Unit(0)).lpNorm<2>());
						 
						}
#if AMREX_SPACEDIM>1
					if (geom[0].isPeriodic(1))
						{
							d1 = std::min( (x-packedspheres[n] + size(1)*0.5*Set::Vector::Unit(1)).lpNorm<2>(),
										  (x-packedspheres[n] - size(1)*0.5*Set::Vector::Unit(1)).lpNorm<2>());
						}
#endif
#if AMREX_SPACEDIM>2
					if (geom[0].isPeriodic(2))
						{
							d1 = std::min( (x-packedspheres[n] + size(2)*Set::Vector::Unit(2)).lpNorm<2>(),
										  (x-packedspheres[n] - size(2)*Set::Vector::Unit(2)).lpNorm<2>());
						}
#endif
                 
						}        
				     
					 

				      d= std::min(d,d1); 

                     

					if (d<R_new[n])
						{
							min_grain_id = 0;
							
						}
					
							
				 }
				}
				}

				if (type == Type::Values) field(i,j,k) = alpha[min_grain_id];
				else if (type == Type::Partition) field(i,j,k,min_grain_id) = alpha[min_grain_id];
			
			});
		}
	}
	
private:
	int number_of_grains;
	int pack_size;
	std::vector<Set::Scalar> alpha;
	std::vector<Set::Vector> packedspheres;
	std::vector <amrex::Real> R_new;
	std::vector <amrex::Real> R;
	std::vector <amrex::Real> dd;
	std::vector <amrex::Real> R_n;
	std::vector <amrex::Real> count;
	std::vector<double>::iterator it;
	//std::vector<double> ins;
    amrex::Real vf;
    amrex::Real R_min;
	amrex::Real R_max;
	int rank;
	Type type;
	amrex::Vector<amrex::Real> packedspheres_x;
	amrex::Vector<amrex::Real> packedspheres_y;
#if BL_SPACEDIM==3
	amrex::Vector<amrex::Real> packedspheres_z;
#endif
};
}
#endif
