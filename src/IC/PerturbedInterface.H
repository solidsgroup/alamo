#ifndef IC_PERTURBEDINTERFACE_H_
#define IC_PERTURBEDINTERFACE_H_

#include "IC/IC.H"
#include "Util/Util.H"

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
	PerturbedInterface (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		amrex::ParmParse pp("ic");

		std::vector<std::string> wave_numbers_str;
		pp.queryarr("wave_numbers",wave_numbers_str);
		for (unsigned int i = 0; i<wave_numbers_str.size(); ++i)
		{
			wave_numbers.push_back(Util::String::Parse<std::complex<int> >(wave_numbers_str[i]));
			phis.push_back(std::atan2(wave_numbers[i].imag(),wave_numbers[i].real()));
		}
		pp.queryarr("wave_amplitudes",wave_amplitudes);

		if (wave_numbers.size() != wave_amplitudes.size())
			Util::Abort(INFO, "Number of wave numbers and amplitudes must match");
	}
  
	void Add(const int lev,
		 amrex::Vector<amrex::MultiFab *> &field)
	{
		bool cellcentered = (field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		Set::Scalar L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];

		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
				     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
				     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Real x1, x2;

				if (cellcentered)
				{
					x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
					x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
				}
				else
				{
					x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
					x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
				}

				amrex::Real bdry = 0.;
				for (int n = 0; n < wave_numbers.size(); n++)
					bdry += wave_amplitudes[n]
						* (fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*x1 / L1) +
						   fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*x1 / L1));

				if (x2 < bdry)
				{
					field_box(m,0) = 1.;     
					field_box(m,1) = 0.;     
				}
				else
				{
					field_box(m,0) = 0.;     
					field_box(m,1) = 1.;     
				}
			}
		}

	};
  
private:
	amrex::Vector<std::complex<int> > wave_numbers; ///< Store mode amplitudes \f$A_n\f$
	amrex::Vector<Set::Scalar> wave_amplitudes;
	std::vector<Set::Scalar> phis;
};
}
#endif
