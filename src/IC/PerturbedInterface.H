#ifndef IC_PERTURBEDINTERFACE_H_
#define IC_PERTURBEDINTERFACE_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

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
	enum Mollifier {Dirac, Gaussian};
	PerturbedInterface (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
	}
  
	void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
	{
		Set::Scalar AMREX_D_DECL(L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0],
								 L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1],
								 L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2]);

		for (amrex::MFIter mfi(*a_field[lev],true); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::IndexType type = a_field[lev]->ixType();

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				
				Set::Vector x;
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
                }
                else if (type == amrex::IndexType::TheCellType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                }

				Set::Scalar bdry = 0.;
#if AMREX_SPACEDIM == 2				
				Set::Scalar s1=NAN, l1=NAN;
				if (normal==Direction::X) {s1 = x(1); l1 = L2;}
				if (normal==Direction::Y) {s1 = x(0); l1 = L1;}
				for (int n = 0; n < wave_numbers.size(); n++)
					bdry += wave_amplitudes[n]
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s1 / l1) +
							fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s1 / l1));
#elif AMREX_SPACEDIM == 3				
				Set::Scalar s1=NAN, s2=NAN, l1=NAN, l2=NAN;
				if (normal==Direction::X) {s1 = x(1); s2 = x(2); l1 = L2; l2 = L3;}
				if (normal==Direction::Y) {s1 = x(2); s2 = x(0); l1 = L3; l2 = L1;}
				if (normal==Direction::Z) {s1 = x(0); s2 = x(0); l1 = L1; l2 = L2;}
				for (int n = 0; n < wave_numbers.size(); n++)
					bdry += wave_amplitudes[n]
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s1 / l1) +
							fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s1 / l1))
						* ( fabs(std::cos(phis[n]))*std::cos(wave_numbers[n].real()*Set::Constant::Pi*s2 / l2) +
					 	    fabs(std::sin(phis[n]))*std::sin(wave_numbers[n].imag()*Set::Constant::Pi*s2 / l2))
					;
#endif
				if (mol == Mollifier::Dirac)
				{
					// Util::Message(INFO);
					if ((normal == Direction::X && x(0) < bdry + offset)||
						(normal == Direction::Y && x(1) < bdry + offset)||
						(normal == Direction::Z && x(2) < bdry + offset))
					{
						field(i,j,k,0) = 1.;     
						field(i,j,k,1) = 0.;     
					}
					else
					{
						field(i,j,k,0) = 0.;     
						field(i,j,k,1) = 1.;     
					}
				}
				else
				{
					Set::Scalar t = 0.0;
					if (normal == Direction::X) t = x(0) - bdry - offset;
					else if (normal == Direction::Y) t = x(1) - bdry - offset;
					else if (normal == Direction::Z) t = x(2) - bdry - offset;

					Set::Scalar value = 0.5 + 0.5*std::erf(t/eps);
					// Util::Message(INFO, "value = ", value);
					field(i,j,k,0) = value;
					field(i,j,k,1) = 1. - value;

					if (field(i,j,k,0) < 0.0) field(i,j,k,0) = 0.0;
					if (field(i,j,k,0) > 1.0) field(i,j,k,0) = 1.0;
					if (field(i,j,k,1) < 0.0) field(i,j,k,1) = 0.0;
					if (field(i,j,k,1) > 1.0) field(i,j,k,1) = 1.0;
				}
			});
		}

	};
  
private:
	enum Direction {X,Y,Z};
	Direction normal = Direction::Y;
	Set::Scalar offset = 0.0;
	amrex::Vector<std::complex<int> > wave_numbers; ///< Store mode amplitudes \f$A_n\f$
	amrex::Vector<Set::Scalar> wave_amplitudes;
	std::vector<Set::Scalar> phis;
	Mollifier mol = Mollifier::Gaussian;
	Set::Scalar eps;

public:
	static void Parse(PerturbedInterface & value, IO::ParmParse & pp)
    {
		std::vector<std::string> wave_numbers_str;
		pp.queryarr("wave_numbers",wave_numbers_str);
		for (unsigned int i = 0; i<wave_numbers_str.size(); ++i)
		{
			value.wave_numbers.push_back(Util::String::Parse<std::complex<int> >(wave_numbers_str[i]));
			value.phis.push_back(std::atan2(value.wave_numbers[i].imag(),value.wave_numbers[i].real()));
		}
		pp.queryarr("wave_amplitudes",value.wave_amplitudes);

		std::string normal_str;
		pp.query("normal",normal_str);
		if (normal_str == "x") value.normal = Direction::X;
		if (normal_str == "y") value.normal = Direction::Y;
		if (normal_str == "z") value.normal = Direction::Z;
		std::string offset_str;
		pp.query("offset",value.offset);
		
		if (value.wave_numbers.size() != value.wave_amplitudes.size())
			Util::Abort(INFO, "Number of wave numbers and amplitudes must match");

		std::string mollifier;
		pp.query("mollifier", mollifier);

		if(mollifier == "dirac")
			value.mol = Mollifier::Dirac;
		else
		{
			value.mol = Mollifier::Gaussian;
			pp.query("eps",value.eps);
			Util::Message(INFO, "eps = ", value.eps);
		}
    }

};
}
#endif
