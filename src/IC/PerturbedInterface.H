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

		std::string normal_str;
		pp.query("normal",normal_str);
		if (normal_str == "x") normal = Direction::X;
		if (normal_str == "y") normal = Direction::Y;
		if (normal_str == "z") normal = Direction::Z;
		std::string offset_str;
		pp.query("offset",offset);
		
		if (wave_numbers.size() != wave_amplitudes.size())
			Util::Abort(INFO, "Number of wave numbers and amplitudes must match");
	}
  
	void Add(const int lev,
		 amrex::Vector<amrex::MultiFab *> &field)
	{
		bool cellcentered = (field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		Set::Scalar AMREX_D_DECL(L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0],
								 L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1],
								 L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2]);

		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
						 for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
						 for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Real x1=NAN, x2=NAN, x3=NAN;

				if (cellcentered)
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
								 x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
								 x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
				}
				else
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
								 x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
								 x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
				}

				amrex::Real bdry = 0.;
#if AMREX_SPACEDIM == 2				
				Set::Scalar s1=NAN, l1=NAN;
				if (normal==Direction::X) {s1 = x2; l1 = L2;}
				if (normal==Direction::Y) {s1 = x1; l1 = L1;}
				for (int n = 0; n < wave_numbers.size(); n++)
					bdry += wave_amplitudes[n]
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s1 / l1) +
							fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s1 / l1));
#elif AMREX_SPACEDIM == 3				
				Set::Scalar s1=NAN, s2=NAN, l1=NAN, l2=NAN;
				if (normal==Direction::X) {s1 = x2; s2 = x3; l1 = L2; l2 = L3;}
				if (normal==Direction::Y) {s1 = x3; s2 = x1; l1 = L3; l2 = L1;}
				if (normal==Direction::Z) {s1 = x1; s2 = x1; l1 = L1; l2 = L2;}
				for (int n = 0; n < wave_numbers.size(); n++)
					bdry += wave_amplitudes[n]
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s1 / l1) +
							fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s1 / l1))
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s2 / l2) +
					 	    fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s2 / l2))
					;
#endif
				
				if ((normal == Direction::X && x1 < bdry + offset)||
					(normal == Direction::Y && x2 < bdry + offset)||
					(normal == Direction::Z && x3 < bdry + offset))
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
	enum Direction {X,Y,Z};
	Direction normal = Direction::Y;
	Set::Scalar offset = 0.0;
	amrex::Vector<std::complex<int> > wave_numbers; ///< Store mode amplitudes \f$A_n\f$
	amrex::Vector<Set::Scalar> wave_amplitudes;
	std::vector<Set::Scalar> phis;
};
}
#endif
