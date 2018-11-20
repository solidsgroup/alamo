#ifndef IC_CIRCLE_H_
#define IC_CIRCLE_H_

#include "IC/IC.H"
#include "Util/Util.H"

/// \class Circle
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
class Circle : public IC
{
public:
	Circle (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		amrex::ParmParse pp("ic");
		pp.queryarr("center",center);
		pp.query("radius",radius);
	}
  
	void Initialize(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
	{
		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++),
				     for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++),
				     for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++))
					{
						amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
						amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];

						if ( ((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1])) < radius*radius)
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
	amrex::Vector<amrex::Real> center;
	amrex::Real radius;
};
}
#endif
