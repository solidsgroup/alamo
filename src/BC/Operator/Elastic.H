#ifndef BC_OPERATOR_ELASTIC_H
#define BC_OPERATOR_ELASTIC_H

// #include "Operator/Elastic.H"
#include "IO/ParmParse.H"

namespace BC
{
namespace Operator
{
template<class T>
class Elastic
{
public:
	static const int min = 1, max = 2;

 	// friend     ::Operator::Elastic<T>;
	// friend Test::Operator::Elastic;

public:
	/// The four types of boundary conditions are listed and evauluated here
	///
	enum Type {Displacement, Traction, Periodic, Neumann}; 
	enum Face {AMREX_D_DECL(XLO=0,YLO=1,ZLO=2),
		   AMREX_D_DECL(XHI=AMREX_SPACEDIM+0,YHI=AMREX_SPACEDIM+1,ZHI=AMREX_SPACEDIM+2)}; 
	enum Direction {AMREX_D_DECL(X=0,Y=1,Z=2)}; 

 	Elastic()
	{
		// By default, all boundary conditions are displacement
		for (int face = 0; face < AMREX_SPACEDIM; face++)
			for (int direction = 0; direction < AMREX_SPACEDIM; direction++)
			{
				m_bc_lo[face][direction] = Type::Displacement;
				m_bc_lo_val[face][direction] = 0.0;
				m_bc_hi[face][direction] = Type::Displacement;
				m_bc_hi_val[face][direction] = 0.0;
			}
	};
 	~Elastic() {};
	
	void
	Set(const Face face,
	      const Direction direction,
	      const Type type,
	      const Set::Scalar value)
	{
		if (AMREX_D_TERM(face == Face::XLO, || face == Face::YLO, ||  face == Face::ZLO))
		{
			m_bc_lo[face%AMREX_SPACEDIM][direction] = type;
			m_bc_lo_val[face%AMREX_SPACEDIM][direction] = value;
		}
		if (AMREX_D_TERM(face == Face::XHI, || face == Face::YHI, ||  face == Face::ZHI))
		{
			m_bc_hi[face%AMREX_SPACEDIM][direction] = type;
			m_bc_hi_val[face%AMREX_SPACEDIM][direction] = value;
		}
	}

	void
	Set(const Face face,
	      const Direction direction,
	      const Type type,
	      const Set::Scalar value,
	      amrex::Vector<amrex::MultiFab *> &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		if (AMREX_D_TERM(face == Face::XLO, || face == Face::YLO, ||  face == Face::ZLO))
		{
			m_bc_lo[face%AMREX_SPACEDIM][direction] = type;
			m_bc_lo_val[face%AMREX_SPACEDIM][direction] = value;
		}
		if (AMREX_D_TERM(face == Face::XHI, || face == Face::YHI, ||  face == Face::ZHI))
		{
			m_bc_hi[face%AMREX_SPACEDIM][direction] = type;
			m_bc_hi_val[face%AMREX_SPACEDIM][direction] = value;
		}

		for (int ilev = 0; ilev < a_rhs.size(); ilev++)
		{
			amrex::Box domain(a_geom[ilev].Domain());
			domain.convert(amrex::IntVect::TheNodeVector());
			const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			for (amrex::MFIter mfi(*a_rhs[ilev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
			{
				amrex::Box bx = mfi.tilebox();
				bx.grow(2);
				bx = bx & domain;
				amrex::Array4<amrex::Real> const& rhs       = a_rhs[ilev]->array(mfi);
				amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
						bool cont = true;
						if (type != Type::Displacement)
						{
							AMREX_D_PICK(
								if(i == lo.x || i == hi.x) //end points of a line
									cont = false;
								,
								if((i == lo.x || i == hi.x) && (j == lo.y || j == hi.y)) // corners of the square
									cont = false;
								,
								if ((i == lo.x || i == hi.x) && (j == lo.y || j == hi.y) && (k == lo.z || k == hi.z))
									cont = false;	// corners of the cube
								else if ((j == lo.y || j == hi.y) && (k == lo.z || k == hi.z))
						 			cont = false; // edge of the cube
						 		else if ((i == lo.x || i == hi.x) && (k == lo.z || k == hi.z))
						 			cont = false; // edge of the cube
						 		else if ((i == lo.x || i == hi.x) && (j == lo.y || j == hi.y))
						 			cont = false; // edge of the cube
							);
						 	//if ((i == lo.x || i == hi.x) && (j == lo.y || j == hi.y) && (k == lo.z || k == hi.z))
						 	//	cont = false; // corner
						 	//else if ((j == lo.y || j == hi.y) && (k == lo.z || k == hi.z))
						 	//	cont = false; // edge
						 	//else if ((i == lo.x || i == hi.x) && (k == lo.z || k == hi.z))
						 	//	cont = false; // edge
						 	//else if ((i == lo.x || i == hi.x) && (j == lo.y || j == hi.y))
						 	//	cont = false; // edge
						}
						if (cont)
						{
							AMREX_D_TERM(if (i == lo.x && face == Face::XLO) rhs(i,j,k,direction) = value;,
								     if (j == lo.y && face == Face::YLO) rhs(i,j,k,direction) = value;,
								     if (k == lo.z && face == Face::ZLO) rhs(i,j,k,direction) = value;);
							AMREX_D_TERM(if (i == hi.x && face == Face::XHI) rhs(i,j,k,direction) = value;,
								     if (j == hi.y && face == Face::YHI) rhs(i,j,k,direction) = value;,
								     if (k == hi.z && face == Face::ZHI) rhs(i,j,k,direction) = value;);
						}
					});
			}					
		}
	}

	void
	Set(const Face face,
	    const Direction direction,
	    const Type type,
	    const Set::Scalar value,
	    amrex::Vector<amrex::MultiFab> &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
		Set(face,direction,type,value,pa_rhs,a_geom);
	}

	void
	Set(const Face face,
	    const Direction direction,
	    const Type type,
	    const Set::Scalar value,
	    amrex::Vector<std::unique_ptr<amrex::MultiFab> > &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
		Set(face,direction,type,value,pa_rhs,a_geom);
	}

	void
	Init(amrex::Vector<amrex::MultiFab *> &a_rhs,
	     const amrex::Vector<amrex::Geometry> &a_geom)
	{
		#if AMREX_SPACEDIM > 0
		Set(Face::XLO,Direction::X,m_bc_lo[Face::XLO%AMREX_SPACEDIM][Direction::X],m_bc_lo_val[Face::XLO%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		Set(Face::XHI,Direction::X,m_bc_hi[Face::XHI%AMREX_SPACEDIM][Direction::X],m_bc_hi_val[Face::XHI%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		#endif
		#if AMREX_SPACEDIM > 1
		Set(Face::XLO,Direction::Y,m_bc_lo[Face::XLO%AMREX_SPACEDIM][Direction::Y],m_bc_lo_val[Face::XLO%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		Set(Face::XHI,Direction::Y,m_bc_hi[Face::XHI%AMREX_SPACEDIM][Direction::Y],m_bc_hi_val[Face::XHI%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		Set(Face::YLO,Direction::X,m_bc_lo[Face::YLO%AMREX_SPACEDIM][Direction::X],m_bc_lo_val[Face::YLO%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		Set(Face::YLO,Direction::Y,m_bc_lo[Face::YLO%AMREX_SPACEDIM][Direction::Y],m_bc_lo_val[Face::YLO%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		Set(Face::YHI,Direction::X,m_bc_hi[Face::YHI%AMREX_SPACEDIM][Direction::X],m_bc_hi_val[Face::YHI%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		Set(Face::YHI,Direction::Y,m_bc_hi[Face::YHI%AMREX_SPACEDIM][Direction::Y],m_bc_hi_val[Face::YHI%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		#endif
		#if AMREX_SPACEDIM > 2
		Set(Face::XLO,Direction::Z,m_bc_lo[Face::XLO%AMREX_SPACEDIM][Direction::Z],m_bc_lo_val[Face::XLO%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		Set(Face::XHI,Direction::Z,m_bc_hi[Face::XHI%AMREX_SPACEDIM][Direction::Z],m_bc_hi_val[Face::XHI%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		Set(Face::YLO,Direction::Z,m_bc_lo[Face::YLO%AMREX_SPACEDIM][Direction::Z],m_bc_lo_val[Face::YLO%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		Set(Face::YHI,Direction::Z,m_bc_hi[Face::YHI%AMREX_SPACEDIM][Direction::Z],m_bc_hi_val[Face::YHI%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		Set(Face::ZLO,Direction::X,m_bc_lo[Face::ZLO%AMREX_SPACEDIM][Direction::X],m_bc_lo_val[Face::ZLO%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		Set(Face::ZLO,Direction::Y,m_bc_lo[Face::ZLO%AMREX_SPACEDIM][Direction::Y],m_bc_lo_val[Face::ZLO%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		Set(Face::ZLO,Direction::Z,m_bc_lo[Face::ZLO%AMREX_SPACEDIM][Direction::Z],m_bc_lo_val[Face::ZLO%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		Set(Face::ZHI,Direction::X,m_bc_hi[Face::ZHI%AMREX_SPACEDIM][Direction::X],m_bc_hi_val[Face::ZHI%AMREX_SPACEDIM][Direction::X],a_rhs,a_geom);
		Set(Face::ZHI,Direction::Y,m_bc_hi[Face::ZHI%AMREX_SPACEDIM][Direction::Y],m_bc_hi_val[Face::ZHI%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);
		Set(Face::ZHI,Direction::Z,m_bc_hi[Face::ZHI%AMREX_SPACEDIM][Direction::Z],m_bc_hi_val[Face::ZHI%AMREX_SPACEDIM][Direction::Z],a_rhs,a_geom);
		#endif
	}

	void
	Init(amrex::Vector<amrex::MultiFab> &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
		Init(pa_rhs,a_geom);
	}

	void
	Init(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
		Init(pa_rhs,a_geom);
	}


 	AMREX_FORCE_INLINE
 	Set::Vector operator () (const Set::Vector &u,
 				 const Set::Matrix &gradu,
 				 const Set::Matrix &sigma,
 				 const int &i, const int &j, const int &k,
				 const amrex::Box &domain)
 	{
		(void)i; (void)j; (void)k; // Suppress "unused variable" warnings
		Set::Vector f;

		const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
		



		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		for (int p = 0; p < AMREX_SPACEDIM; p++) // iterate over DIMENSIONS
		{
			// corners
			AMREX_D_TERM(
			if (AMREX_D_TERM( i==lo.x, && j==lo.y, && k==lo.z) &&
			    (AMREX_D_TERM(m_bc_lo[0][p] == Type::Displacement, || m_bc_lo[1][p] == Type::Displacement, || m_bc_lo[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (AMREX_D_TERM(i==hi.x,&& j==lo.y,&& k==lo.z) &&
			    (AMREX_D_TERM(m_bc_hi[0][p] == Type::Displacement, || m_bc_lo[1][p] == Type::Displacement, || m_bc_lo[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			},
			if (AMREX_D_TERM(i==lo.x, && j==hi.y, && k==lo.z) &&
			    (AMREX_D_TERM(m_bc_lo[0][p] == Type::Displacement, || m_bc_hi[1][p] == Type::Displacement, || m_bc_lo[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (AMREX_D_TERM(i==hi.x,&& j==hi.y, && k==lo.z) &&
			    (AMREX_D_TERM(m_bc_hi[0][p] == Type::Displacement, || m_bc_hi[1][p] == Type::Displacement, || m_bc_lo[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			},
			if (AMREX_D_TERM(i==lo.x, && j==lo.y, && k==hi.z) &&
			    (AMREX_D_TERM(m_bc_lo[0][p] == Type::Displacement, || m_bc_lo[1][p] == Type::Displacement, || m_bc_hi[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (AMREX_D_TERM(i==lo.x, && j==hi.y, && k==hi.z) &&
			    (AMREX_D_TERM(m_bc_lo[0][p] == Type::Displacement, || m_bc_hi[1][p] == Type::Displacement, || m_bc_hi[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (AMREX_D_TERM(i==hi.x, && j==lo.y, && k==hi.z) &&
			    (AMREX_D_TERM(m_bc_hi[0][p] == Type::Displacement, || m_bc_lo[1][p] == Type::Displacement, || m_bc_hi[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (AMREX_D_TERM(i==hi.x, && j==hi.y, && k==hi.z) &&
			    (AMREX_D_TERM(m_bc_hi[0][p] == Type::Displacement, || m_bc_hi[1][p] == Type::Displacement, || m_bc_hi[2][p] == Type::Displacement)))
			{
			 	f(p) = u(p);
			 	continue;
			});

			// Displacement edges
			AMREX_D_TERM(
			,
			if (i==lo.x                                && j==lo.y && 
			    (m_bc_lo[0][p] == Type::Displacement || m_bc_lo[1][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==lo.x                                && j==hi.y && 
			    (m_bc_lo[0][p] == Type::Displacement || m_bc_hi[1][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==hi.x                                && j==lo.y && 
			    (m_bc_hi[0][p] == Type::Displacement || m_bc_lo[1][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==hi.x                                && j==hi.y && 
			    (m_bc_hi[0][p] == Type::Displacement || m_bc_hi[1][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			,
			if (i==lo.x                                && k==hi.z && 
			    (m_bc_lo[0][p] == Type::Displacement || m_bc_hi[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==lo.x                                && k==lo.z && 
			    (m_bc_lo[0][p] == Type::Displacement || m_bc_lo[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==hi.x                                && k==hi.z && 
			    (m_bc_hi[0][p] == Type::Displacement || m_bc_hi[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (i==hi.x                                && k==lo.z && 
			    (m_bc_hi[0][p] == Type::Displacement || m_bc_lo[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (j==lo.y                                && k==lo.z && 
			    (m_bc_lo[1][p] == Type::Displacement || m_bc_lo[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (j==lo.y                                && k==hi.z && 
			    (m_bc_lo[1][p] == Type::Displacement || m_bc_hi[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (j==hi.y                                && k==lo.z && 
			    (m_bc_hi[1][p] == Type::Displacement || m_bc_lo[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			if (j==hi.y                                && k==hi.z && 
			    (m_bc_hi[1][p] == Type::Displacement || m_bc_hi[2][p] == Type::Displacement))
			{
			 	f(p) = u(p);
			 	continue;
			}
			);

			AMREX_D_TERM(
			,
			// X-Y traction edges...
			if (i == hi.x                       && j == hi.y &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_hi[1][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(1.0,1.0,0.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == lo.x                       && j == hi.y &&
			    m_bc_lo[0][p] == Type::Traction && m_bc_hi[1][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(-1.0,1.0,0.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == hi.x                       && j == lo.y &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_lo[1][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(1.0,-1.0,0.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == lo.x                       && j == lo.y &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_lo[1][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(-1.0,-1.0,0.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			,
			// Y-Z traction edges...
			if (j == hi.y                       && k == hi.z &&
			    m_bc_hi[1][p] == Type::Traction && m_bc_hi[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(0,1.0,1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (j == lo.y                       && k == hi.z &&
			    m_bc_lo[1][p] == Type::Traction && m_bc_hi[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(0,-1.0,1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (j == hi.y                       && k == lo.z &&
			    m_bc_hi[1][p] == Type::Traction && m_bc_lo[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(0,1.0,-1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (j == lo.y                       && k == lo.z &&
			    m_bc_hi[1][p] == Type::Traction && m_bc_lo[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(0,-1.0,-1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}

			// X-Z traction edges...
			if (i == hi.x                       && k == hi.z &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_hi[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(1.0,0.0,1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == lo.x                       && k == hi.z &&
			    m_bc_lo[0][p] == Type::Traction && m_bc_hi[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(-1.0,0.0,1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == hi.x                       && k == lo.z &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_lo[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(1.0,0.0,-1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			if (i == lo.x                       && k == lo.z &&
			    m_bc_hi[0][p] == Type::Traction && m_bc_lo[2][p] == Type::Traction)
			{
				Set::Vector n(AMREX_D_DECL(-1.0,0.0,-1.0));
				Set::Vector t = sigma*n;
				f(p) = t(p);
				continue;
			}
			);

			for (int q = 0; q < AMREX_SPACEDIM; q++) // iterate over FACES
			{
				if (m[q] == domain.loVect()[q])
				{
					if      (m_bc_lo[q][p] == Type::Displacement)
						f(p) =   u(p);
					else if (m_bc_lo[q][p] == Type::Traction)
						f(p) = -sigma(p,q);
					else if (m_bc_lo[q][p] == Type::Neumann)
						f(p) = -gradu(p,q);
					else if (m_bc_lo[q][p] == Type::Periodic)
						continue;
					else
						Util::Abort(INFO, "Invalid BC");
				}
				if (m[q] == domain.hiVect()[q])
				{
					if      (m_bc_hi[q][p] == Type::Displacement)
						f(p) = u(p);
					else if (m_bc_hi[q][p] == Type::Traction)
						f(p) = +sigma(p,q);
					else if (m_bc_hi[q][p] == Type::Neumann)
						f(p) = +gradu(p,q);
					else if (m_bc_hi[q][p] == Type::Periodic)
						continue;
					else Util::Abort(INFO, "Invalid BC");
				}
			}
		}
		// bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
		// 	AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));
		
 		return f;
 	}


private:
	// m_bc_lo[face][dimension]
	std::array<std::array<Type,AMREX_SPACEDIM>, AMREX_SPACEDIM>        m_bc_lo; 
	std::array<std::array<Set::Scalar,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_lo_val; 
	// m_bc_hi[face][dimension]
	std::array<std::array<Type,AMREX_SPACEDIM>, AMREX_SPACEDIM>        m_bc_hi; 
	std::array<std::array<Set::Scalar,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_hi_val;

public:
    static void Parse(Elastic & value, IO::ParmParse & pp)
    {
		std::map<std::string, Type> bcmap;
		bcmap["displacement"] = Type::Displacement;
		bcmap["disp"]         = Type::Displacement;
		bcmap["neumann"]      = Type::Neumann;
		bcmap["traction"]     = Type::Traction;
		bcmap["trac"]         = Type::Traction;
		bcmap["periodic"]     = Type::Periodic;
		
		std::vector<std::string> str;
		if (pp.contains("type.xlo")) { pp.queryarr("type.xlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo[Face::XLO%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhi")) { pp.queryarr("type.xhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi[Face::XHI%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		#if AMREX_SPACEDIM > 1
		if (pp.contains("type.ylo")) { pp.queryarr("type.ylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo[Face::YLO%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		if (pp.contains("type.yhi")) { pp.queryarr("type.yhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi[Face::YHI%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		#endif
		#if AMREX_SPACEDIM > 2
		if (pp.contains("type.zlo")) { pp.queryarr("type.zlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo[Face::ZLO%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		if (pp.contains("type.zhi")) { pp.queryarr("type.zhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi[Face::ZHI%AMREX_SPACEDIM][i] = bcmap[str[i]]; }
		#endif

		std::vector<Set::Scalar> val;
		if (pp.contains("val.xlo"))  { pp.queryarr("val.xlo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo_val[Face::XLO%AMREX_SPACEDIM][i] = val[i];}
		if (pp.contains("val.xhi"))  { pp.queryarr("val.xhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi_val[Face::XHI%AMREX_SPACEDIM][i] = val[i];}
		#if AMREX_SPACEDIM > 1
		if (pp.contains("val.ylo"))  { pp.queryarr("val.ylo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo_val[Face::YLO%AMREX_SPACEDIM][i] = val[i];}
		if (pp.contains("val.yhi"))  { pp.queryarr("val.yhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi_val[Face::YHI%AMREX_SPACEDIM][i] = val[i];}
		#endif
		#if AMREX_SPACEDIM > 2
		if (pp.contains("val.zlo"))  { pp.queryarr("val.zlo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_lo_val[Face::ZLO%AMREX_SPACEDIM][i] = val[i];}
		if (pp.contains("val.zhi"))  { pp.queryarr("val.zhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_hi_val[Face::ZHI%AMREX_SPACEDIM][i] = val[i];}
		#endif
		//if (pp.contains("val.xlo")) for (int i = 0; i < AMREX_SPACEDIM; i++) 

    }
};

}
}
#endif

