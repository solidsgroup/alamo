#ifndef PerturbedInterface_H_
#define PerturbedInterface_H_

#include "GeneralAMRIntegrator/GeneralAMRIntegratorIC.H"

class PerturbedInterface : public GeneralAMRIntegratorIC
{
public:
  PerturbedInterface (amrex::Array<amrex::Geometry> &_geom) :
    GeneralAMRIntegratorIC(_geom)
  {
      amrex::ParmParse pp("ic");
      pp.queryarr("wave_numbers",wave_numbers);
      pp.queryarr("wave_amplitudes",wave_amplitudes);

      if (wave_numbers.size() != wave_amplitudes.size())
	amrex::Abort("Number of wave numbers and amplitudes must match");

  }
  
void Initialize(const int lev,
			  amrex::Array<std::unique_ptr<amrex::MultiFab> > &field)
  {
    const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];

    for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<Real> &field_box = (*field[lev])[mfi];

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
		amrex::Real bdry = 0.;
		for (int n = 0; n < wave_numbers.size(); n++)
		  bdry += wave_amplitudes[n] * sin((amrex::Real)(wave_numbers[n])*x*pi/width);

		// amrex::Real bdry  =
		//   (+ 0.1 * sin(x*(2*pi)/width)
		//    + 0.1 * cos(2.0*x*(2*pi)/width)
		//    + 0.1 * sin(3.0*x*(2*pi)/width)
		//    + 0.1 * cos(4.0*x*(2*pi)/width)
		//    + 0.1 * sin(5.0*x*(2*pi)/width)
		//    + 0.15 * cos(6.0*x*(2*pi)/width)
		//    + 0.15 * sin(7.0*x*(2*pi)/width)
		//    + 0.15 * cos(8.0*x*(2*pi)/width)
		//    + 0.15 * sin(9.0*x*(2*pi)/width)
		//    + 0.15 * cos(10.0*x*(2*pi)/width))*0.1
		//   ;

		if (y < bdry)
		  {
		    field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 1.;     
		    field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.;     
		  }
		else
		  {
		    field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.;     
		    field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.;     
		  }
	      }
      }

  };
  
private:
  static constexpr amrex::Real pi = 3.14159265359;
  amrex::Vector<int> wave_numbers;
  amrex::Vector<amrex::Real> wave_amplitudes;
};

#endif
