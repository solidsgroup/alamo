#ifndef IC_SPHERE_H_
#define IC_SPHERE_H_

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"

/// \class Sphere
/// \brief Initialize a spherical inclusion
namespace IC
{
class Sphere : public IC
{
public:
/// Sphere is also used for a cylindrical case, but switching from the default type to the desired orientation.
/// For instance, XY refers to a cylinder with the circular faces in the XY plane.
	enum Type {XYZ, XY, YZ, XZ};

	Sphere (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
	/// Constructor defining radius, center, dimension/orientation, and field values within and outside of the area.
	Sphere (amrex::Vector<amrex::Geometry> &_geom, Set::Scalar _radius, amrex::Vector<Set::Scalar> _center, Type _type = Type::XYZ ,Set::Scalar _alpha_in = 1, Set::Scalar _alpha_out = 0) : IC(_geom)
	{
		Define(_radius, _center, _type,_alpha_in,_alpha_out);
	}

	void Define(Set::Scalar _radius,
				amrex::Vector<amrex::Real> _center,
				Type _type,
				Set::Scalar _alpha_in,
				Set::Scalar _alpha_out)
	{
		radius = _radius;
		center = _center;
		type = _type;
		alpha_in = _alpha_in;
		alpha_out = _alpha_out;
	}
  
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab * > &a_field)
	{
		bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());

			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

					Set::Scalar AMREX_D_DECL(x,y,z);
					if (cellcentered)
					{
						AMREX_D_TERM(x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
									 y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
									 z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
					}
					else
					{
						AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
									 y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
									 z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];);
					}

					Set::Scalar rsq = NAN;

					if (type == Type::XYZ)
					{
						// 3D Sphere
						rsq =
						AMREX_D_TERM(  (x-center[0])*(x-center[0]) , 
									 + (y-center[1])*(y-center[1]) ,
									 + (z-center[2])*(z-center[2]));
					}
					else if (type == Type::XY)
					{
						// Cylinder along Z axis
						rsq =
						AMREX_D_TERM(  (x-center[0])*(x-center[0]) , 
									 + (y-center[1])*(y-center[1]) ,
									 );
					}
					else if (type == Type::YZ)
					{
						// Cylinder along X axis
						rsq = 
						AMREX_D_TERM( , 
									 + (y-center[1])*(y-center[1]) ,
									 + (z-center[2])*(z-center[2]));
					}
					else if (type == Type::XZ)
					{
						// Cylinder along Y axis
						rsq =
						AMREX_D_TERM(  (x-center[0])*(x-center[0]) , 
									 ,
									 + (z-center[2])*(z-center[2]));
					}

					if ( rsq < radius*radius)
						{
							field(i,j,k,0) = alpha_in;     
							field(i,j,k,1) = 0.;     
						}
						else
						{
							field(i,j,k,0) = alpha_out;     
							field(i,j,k,1) = 1;     
						}
				});
		}

	};
  
private:
	int dim = AMREX_SPACEDIM;
	amrex::Vector<Set::Scalar> center;
	Set::Scalar radius,alpha_in,alpha_out;
	Type type;
};
}
#endif
