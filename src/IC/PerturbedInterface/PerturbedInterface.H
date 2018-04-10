#ifndef IC_PERTURBEDINTERFACE_H_
#define IC_PERTURBEDINTERFACE_H_

#include "IC/IC.H"

/// \class PerturbedInterface
/// \brief Initialize a perturbed interface using Fourier Modes
///
/// \todo Extend this class to allow for 3D perturbations, currently only 2D are allowed
/// \todo Allow for cosine (or complex exponential) expansions rather than just sin.
/// \note This is a **two grain only** initial condition.
/// \note This replaces the depricated "perturbed_bar" initial condition from previous versions
///
/// The interface is defined as the \f$x=0\f$ plane (2D), or the \f$x=0,z=0\f$ plane (3D).
/// The equation for the interface is given by
/// \f[y(x,z) = \sum_{n\in \{n_1,\ldots,n_N\}} A_n \sin(n\pi x/L_x) \f]
/// where \f$A_n\f$ are the amplitudes (stored in #wave_amplitudes),
/// \f$n_1,\ldots,n_N\subset\mathbb{Z}_+\f$ are wave numbers (stored in #wave_numbers),
/// and \f$L_x\f$ is the length in the x direction (obtained using the #geom object).
///
/// Grain 1 is defined as being above \f$y(x,z)\f$, Grain 2 is defined as being below.
namespace IC
{
class PerturbedInterface : public IC
{
public:
  PerturbedInterface (amrex::Array<amrex::Geometry> &_geom) :
    IC(_geom)
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
    const amrex::Real width = geom[lev].ProbLo()[0] - geom[lev].ProbHi()[0];
    //const amrex::Real width = geom[0].ProbHi()[0] - geom[0].ProbHi()[1];

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
		// amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
#endif
		amrex::Real bdry = 0.;
		for (int n = 0; n < wave_numbers.size(); n++)
		  bdry += wave_amplitudes[n] * sin((amrex::Real)(wave_numbers[n])*x*pi/width);

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
  amrex::Vector<int> wave_numbers; ///< Store mode amplitudes \f$A_n\f$
  amrex::Vector<amrex::Real> wave_amplitudes;
};
}
#endif
