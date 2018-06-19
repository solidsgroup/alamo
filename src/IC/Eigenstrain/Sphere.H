#ifndef IC_EIGENSTRAIN_SPHERE_H_
#define IC_EIGENSTRAIN_SPHERE_H_

#include "IC/IC.H"

namespace IC
{
namespace Eigenstrain
{
class Sphere : public IC::IC
{
public:
	Sphere (amrex::Vector<amrex::Geometry> &_geom) :
		IC::IC(_geom) {}
  
	void Initialize(const int lev,
						 amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eps0)
	{
		Initialize(lev,*eps0[lev].get());
	}

	void Initialize(const int lev,
						 amrex::Vector<amrex::MultiFab> &eps0)
	{
		Initialize(lev,eps0[lev]);
	}
  
	void Initialize(const int lev,
						 amrex::MultiFab &eps0)
	{
		eps0.setVal(0.0);
		const amrex::Real L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
		const amrex::Real L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];
#if AMREX_SPACEDIM==3
		const amrex::Real L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2];
#endif
		for (amrex::MFIter mfi(eps0,true); mfi.isValid(); ++mfi)
			{
				const amrex::Box& box = mfi.tilebox();

				amrex::FArrayBox &eps0fab = (eps0)[mfi];

				for (int i = box.loVect()[0]-eps0.nGrow(); i<=box.hiVect()[0]+eps0.nGrow(); i++) 
					for (int j = box.loVect()[1]-eps0.nGrow(); j<=box.hiVect()[1]+eps0.nGrow(); j++)
#if AMREX_SPACEDIM==3
						for (int k = box.loVect()[2]-eps0.nGrow(); k<=box.hiVect()[2]+eps0.nGrow(); k++)
#endif
							{
								amrex::IntVect m(AMREX_D_DECL(i,j,k));

								amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
								amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
#if AMREX_SPACEDIM==3
								amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
#endif

								amrex::Real eX = (x - (geom[lev].ProbLo()[0] + 0.5*L1))/(0.25*L1);
								amrex::Real eY = (y - (geom[lev].ProbLo()[1] + 0.5*L2))/(0.25*L2);
#if AMREX_SPACEDIM==3
								amrex::Real eZ = (z - (geom[lev].ProbLo()[2] + 0.5*L3))/(0.25*L3);
#endif

								if (AMREX_D_TERM(eX*eX, + eY*eY, + eZ*eZ) < 1.0)
									eps0fab(m,0) = 1.0;     
								else
								 	eps0fab(amrex::IntVect(AMREX_D_DECL(i,j,k))) = 0.0;     
							}
			}

	}

	


private:
	amrex::Real Tin, Tout;
};
}
}
#endif
