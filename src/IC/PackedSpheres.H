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
/// \brief Create a series of patches using the packedspheres Tesselation
///
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
	PackedSpheres (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains,Set::Scalar a_alpha = 1.0,amrex::Real _volume_fraction =0.05, amrex::Real _std_deviation=0) : IC(_geom) 
	{
    	Define(_number_of_grains,a_alpha,_volume_fraction,_std_deviation);
  	}

  	void Define (int a_number_of_grains, Set::Scalar a_alpha = 1.0,amrex::Real a_volume_fraction=0.05, amrex::Real a_std_deviation =0)
	{
    	Define(a_number_of_grains, std::vector<Set::Scalar>(a_number_of_grains,a_alpha),a_volume_fraction, a_std_deviation,Type::Partition);
  	}

	void Define (int a_number_of_grains,
				 std::vector<Set::Scalar> a_alpha,
				 amrex::Real a_volume_fraction,
                 amrex::Real a_std_deviation,
				 Type a_type = Type::Values)
	{
		number_of_grains = a_number_of_grains;
		alpha = a_alpha;
		vf= a_volume_fraction;
        std_deviation=a_std_deviation;
		type = a_type;

		
		
		//double vf=0.2;
		double R_min=0.1;
		double R_max=0.2;
       

	  
		
	     
		std::vector <amrex::Real> sph_vol(1);   
		std::vector <amrex::Real> sph_nos(1); 
		
		int p=0;
	    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	
		Set::Vector size;
		AMREX_D_TERM(size(0) = geom[0].ProbHi()[0] - geom[0].ProbLo()[0];,
					 size(1) = geom[0].ProbHi()[1] - geom[0].ProbLo()[1];,
					 size(2) = geom[0].ProbHi()[2] - geom[0].ProbLo()[2];)	
					 int total=0;
					 R.resize(number_of_grains); 

        for (int n = 0; n<number_of_grains; n++)
		{
		if (n==0)
		{	
		R.resize(n+1);	
        R[n]=(R_max-R_min)/number_of_grains;
		sph_vol.resize(number_of_grains);
		sph_nos.resize(number_of_grains);
		sph_vol[n]=(4/3)*3.14*(R[n]*R[n]*R[n]);
		sph_nos[n]=ceil(ceil((0.6)/(number_of_grains*sph_vol[n]*100))*vf);
        std::cout<<std::fixed<<std::setprecision(8)<<"R  "<<R[n]<<endl;
		std::cout<<std::fixed<<std::setprecision(8)<<"sphere nos  "<<sph_nos[n]<<endl;
		std::cout<<std::fixed<<std::setprecision(8)<<"sphere vol  "<<sph_vol[n]<<endl;
		}

		else if (n>0)
		{	
		R.resize(n+1);	
        R[n]=R[n-1]+((R_max-R_min)/number_of_grains);
		sph_vol.resize(number_of_grains);
		sph_nos.resize(number_of_grains);
		sph_vol[n]=(4/3)*3.14*(R[n]*R[n]*R[n]);
		sph_nos[n]=ceil(ceil((0.6)/(number_of_grains*sph_vol[n]*100))*vf);
        std::cout<<std::fixed<<std::setprecision(8)<<"R  "<<R[n]<<endl;
		std::cout<<std::fixed<<std::setprecision(8)<<"sphere nos  "<<sph_nos[n]<<endl;
		std::cout<<std::fixed<<std::setprecision(8)<<"sphere vol  "<<sph_vol[n]<<endl;
		}
		}
        
		
		
        pack_size=std::accumulate(sph_nos.begin(),sph_nos.end(),0);
		std::cout<<"pack_size   "<<pack_size<<endl;
		
        packedspheres.resize(pack_size);
       std::vector<double>::iterator it;

	std::vector<double> ins (1,1);
	R_new.assign(R.begin(),R.end());
	
	int inc=0;
	int jj=0;

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

for (int i=0; i<pack_size;i++)
{
	std::cout<<std::fixed<<std::setprecision(8)<<"R_new    "<<R_new[i]<<endl;
}

  	if (rank==0)
		{
        std::vector <amrex::Real> dd(pack_size-1);
		std::vector <amrex::Real> R_n(pack_size-1);
		std::vector <amrex::Real> count(pack_size);
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
	for (int n=0;n<pack_size;n++)
	{
	std::cout<<"packedspheres  "<<packedspheres[n]<<endl;
	}

	for (int n=0;n<pack_size;n++)
	{
	std::cout<<"R_new "<<R_new[n]<<endl;
	}
	R_new.resize(pack_size);
	packedspheres.resize(pack_size);
	MPI_Bcast(R_new.data(), R_new.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                 Set::Vector x;
				 int min_grain_id = -1;
				AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) *geom[lev].CellSize()[0];,
							x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) *geom[lev].CellSize()[1];,
							x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
			    

				for (int n = 0; n<pack_size; n++)
				{ 
			    Set::Scalar d ;
				Set::Scalar d1; 
							
				
                 
				
				for(int v=-2;v<2;v++)	
				{
				 for (int u=-4;u<4;u++) 
				 {

				 
					if ((x(0)> (size(0)*(u*0.125))) && (x(0) < (size(0)*((u+1)*0.125))) && (x(1)> (size(1)*(v*0.25))) && (x(1)<  (size(1)*((v+1)*0.25))))
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
    amrex::Real vf;
    amrex::Real std_deviation;
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
