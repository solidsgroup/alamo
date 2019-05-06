#ifndef BC_OPERATOR_ELASTIC_H
#define BC_OPERATOR_ELASTIC_H

#include "Operator/Elastic.H"

namespace BC
{
namespace Operator
{
template<class T>
class Elastic
{
public:
	static const int min = 1, max = 2;

 	friend     ::Operator::Elastic<T>;
	friend Test::Operator::Elastic;

public:
	enum BCType {Displacement, Traction, Periodic, Neumann}; 

 	Elastic()
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					bndry[i][j][k] = BCType::Traction;
	};
 	~Elastic() {};
	
	void
	SetBC(int x, int y, int z,
	      BCType type, Set::Vector value,
	      amrex::Vector<amrex::MultiFab> &a_rhs,
	      const amrex::Vector<Geometry> &a_geom)
	{
		
		for (int ilev = 0; ilev < a_rhs.size(); ilev++)
		{
			amrex::Box domain(a_geom[ilev].Domain());
			domain.convert(amrex::IntVect::TheNodeVector());

			for (MFIter mfi(a_rhs[ilev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
			{
				Box bx = mfi.tilebox();
				//bx.grow(1);        // Expand to cover first layer of ghost nodes
				//Util::Message(INFO,bx);
				

				amrex::Array4<amrex::Real> const& rhs       = a_rhs[ilev].array(mfi);

				const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			
				amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					});

			}					

		}
	}
	

 	AMREX_FORCE_INLINE
 	Set::Vector operator () (const int &amrlev, const int &mglev,
 				 const Set::Vector &u,
 				 const Set::Vector &gradu,
 				 const T & C,
 				 const int &i, const int &j, const int &k,
				 const Geometry &geom)
 	{
 		return u;
 	}


private:
	BCType bndry[3][3][3];

	// BODY
	// bndry[0][0][0]

	// FACE
	// bndry[xmin][0][0] bndry[0][ymin][0] bndry[0][0][zmin]
	// bndry[xmax][0][0] bndry[0][ymax][0] bndry[0][0][zmax]

	// EDGE
	// bndry[0][ymin][zmin] bndry[0][ymin][zmax] bndry[0][ymax][zmin] bndry[0][ymax][zmax]
	// bndry[ymin][0][zmin] bndry[ymin][0][zmax] bndry[ymax][0][zmin] bndry[ymax][0][zmax]
	// bndry[ymin][zmin][0] bndry[ymin][zmax][0] bndry[ymax][zmin][0] bndry[ymax][zmax][0]

	// CORNER
	// bndry[xmin][ymin][zmin] bndry[xmin][ymin][zmax] bndry[xmin][ymax][zmin] bndry[xmin][ymax][zax]
	// bndry[xmax][ymin][zmin] bndry[xmax][ymin][zmax] bndry[xmax][ymax][zmin] bndry[xmax][ymax][zax]

	// faces
	// edges
};

}
}
#endif
