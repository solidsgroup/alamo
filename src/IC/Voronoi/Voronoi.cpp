#include <AMReX_ParallelDescriptor.H>
#include "Voronoi.H"

namespace IC
{
Voronoi::Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains)
  : IC(_geom), number_of_grains(_number_of_grains)
{
  voronoi_x.resize(number_of_grains);
  voronoi_y.resize(number_of_grains);
#if AMREX_SPACEDIM > 2
  voronoi_z.resize(number_of_grains);
#endif

  if(amrex::ParallelDescriptor::IOProcessor())
    {
      for (int n = 0; n<number_of_grains; n++)
  	{
  	  voronoi_x[n] = geom[0].ProbLo(0) +
	    (geom[0].ProbHi(0)-geom[0].ProbLo(0))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);
  	  voronoi_y[n] = geom[0].ProbLo(1) +
	    (geom[0].ProbHi(1)-geom[0].ProbLo(1))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);
#if AMREX_SPACEDIM > 2
  	  voronoi_z[n] = geom[0].ProbLo(2) +
	    (geom[0].ProbHi(2)-geom[0].ProbLo(2))*((amrex::Real)rand())/((amrex::Real)RAND_MAX);
#endif
  	}
    }  

  for (int n = 0; n<number_of_grains; n++)
    {
      amrex::ParallelDescriptor::Bcast<amrex::Real>(&voronoi_x[n],sizeof(amrex::Real));
      amrex::ParallelDescriptor::Bcast<amrex::Real>(&voronoi_y[n],sizeof(amrex::Real));
#if AMREX_SPACEDIM > 2
      amrex::ParallelDescriptor::Bcast<amrex::Real>(&voronoi_z[n],sizeof(amrex::Real));
#endif
    }

}
  
void Voronoi::Initialize(const int lev,
			 amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
{
  for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box& box = mfi.tilebox();

      amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

      for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++) 
	for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++)
#if BL_SPACEDIM==3
	  for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++)
#endif
	    {
	      amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
	      amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
#if BL_SPACEDIM==3
	      amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
#endif
	      amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
	      int min_grain_id = -1;
	      for (int n = 0; n<number_of_grains; n++)
		{
		  field_box(amrex::IntVect(i,j),n) = 0.; // initialize
		  amrex::Real d = sqrt((x-voronoi_x[n])*(x-voronoi_x[n])
				       + (y-voronoi_y[n])*(y-voronoi_y[n])
#if AMREX_SPACEDIM==3
				       + (z-voronoi_z[n])*(z-voronoi_z[n])
#endif
				       );
		  if (d<min_distance)
		    {
		      min_distance = d;
		      min_grain_id = n;
		    }
		}
	      field_box(amrex::IntVect(i,j),min_grain_id) = 1.;
	    }
    }

}
}
