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
	//enum Face {AMREX_D_DECL(XLO=0,YLO=1,ZLO=2),
	//	   AMREX_D_DECL(XHI=AMREX_SPACEDIM+0,YHI=AMREX_SPACEDIM+1,ZHI=AMREX_SPACEDIM+2)}; 

	enum Face{
		XLO, YLO, ZLO, XHI, YHI, ZHI, // 6

		YLO_ZLO, YLO_ZHI, YHI_ZLO, YHI_ZHI, // 12
		ZLO_XLO, ZLO_XHI, ZHI_XLO, ZHI_XHI,
		XLO_YLO, XLO_YHI, XHI_YLO, XHI_YHI,

		XLO_YLO_ZLO, XLO_YLO_ZHI, XLO_YHI_ZLO, XLO_YHI_ZHI, // 8
		XHI_YLO_ZLO, XHI_YLO_ZHI, XHI_YHI_ZLO, XHI_YHI_ZHI
	};

	enum Direction {AMREX_D_DECL(X=0,Y=1,Z=2)}; 

 	Elastic()
	{
		// By default, all boundary conditions are displacement
		for (int face = 0; face < 26; face++)
			for (int direction = 0; direction < 3; direction++)
			{
				m_bc_type[face][direction] = Type::Displacement;
				m_bc_val [face][direction] = 0.0;
			}
	};
 	~Elastic() {};
	
	void
	Set(const Face face,
	      const Direction direction,
	      const Type type,
	      const Set::Scalar value)
	{
		m_bc_type[face][direction] = type;
		m_bc_val [face][direction] = value;
	}

	void
	Set(const Face face,
	      const Direction direction,
	      const Type type,
	      const Set::Scalar value,
	      amrex::Vector<amrex::MultiFab *> &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		Util::Abort(INFO,"This has been depricated");
		m_bc_type[face][direction] = type;
		m_bc_val [face][direction] = value;

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
						AMREX_D_TERM(if (i == lo.x && face == Face::XLO)
								     rhs(i,j,k,direction) = value;,
							     if (j == lo.y && face == Face::YLO)
								     rhs(i,j,k,direction) = value;,
							     if (k == lo.z && face == Face::ZLO)
								     rhs(i,j,k,direction) = value;)
						AMREX_D_TERM(if (i == hi.x && face == Face::XHI)
								     rhs(i,j,k,direction) = value;,
							     if (j == hi.y && face == Face::YHI)
								     rhs(i,j,k,direction) = value;,
							     if (k == hi.z && face == Face::ZHI)
								     rhs(i,j,k,direction) = value;)
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
					
					for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
					{
						if      (i==lo.x && j==lo.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YLO_ZLO][dir];
						else if (i==lo.x && j==lo.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YLO_ZHI][dir];
						else if (i==lo.x && j==hi.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YHI_ZLO][dir];
						else if (i==lo.x && j==hi.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YHI_ZHI][dir];
						else if (i==hi.x && j==lo.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YLO_ZLO][dir];
						else if (i==hi.x && j==lo.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YLO_ZHI][dir];
						else if (i==hi.x && j==hi.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YHI_ZLO][dir];
						else if (i==hi.x && j==hi.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YHI_ZHI][dir];

						else if (j==lo.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::YLO_ZLO][dir];
						else if (j==lo.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::YLO_ZHI][dir];
						else if (j==hi.y && k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::YHI_ZLO][dir];
						else if (j==hi.y && k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::YHI_ZHI][dir];
						else if (k==lo.z && i==lo.x) rhs(i,j,k,dir) = m_bc_val[Face::ZLO_XLO][dir];
						else if (k==lo.z && i==hi.x) rhs(i,j,k,dir) = m_bc_val[Face::ZLO_XHI][dir];
						else if (k==hi.z && i==lo.x) rhs(i,j,k,dir) = m_bc_val[Face::ZHI_XLO][dir];
						else if (k==hi.z && i==hi.x) rhs(i,j,k,dir) = m_bc_val[Face::ZHI_XHI][dir];
						else if (i==lo.x && j==lo.y) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YLO][dir];
						else if (i==lo.x && j==hi.y) rhs(i,j,k,dir) = m_bc_val[Face::XLO_YHI][dir];
						else if (i==hi.x && j==lo.y) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YLO][dir];
						else if (i==hi.x && j==hi.y) rhs(i,j,k,dir) = m_bc_val[Face::XHI_YHI][dir];

						else if (i==lo.x) rhs(i,j,k,dir) = m_bc_val[Face::XLO][dir];
						else if (i==hi.x) rhs(i,j,k,dir) = m_bc_val[Face::XHI][dir];

						else if (j==lo.y) rhs(i,j,k,dir) = m_bc_val[Face::YLO][dir];
						else if (j==hi.y) rhs(i,j,k,dir) = m_bc_val[Face::YHI][dir];

						else if (k==lo.z) rhs(i,j,k,dir) = m_bc_val[Face::ZLO][dir];
						else if (k==hi.z) rhs(i,j,k,dir) = m_bc_val[Face::ZHI][dir];

//						AMREX_D_TERM(if (i == lo.x && face == Face::XLO)
//								     rhs(i,j,k,direction) = value;,
//							     if (j == lo.y && face == Face::YLO)
//								     rhs(i,j,k,direction) = value;,
//							     if (k == lo.z && face == Face::ZLO)
//								     rhs(i,j,k,direction) = value;)
//						AMREX_D_TERM(if (i == hi.x && face == Face::XHI)
//								     rhs(i,j,k,direction) = value;,
//							     if (j == hi.y && face == Face::YHI)
//								     rhs(i,j,k,direction) = value;,
//							     if (k == hi.z && face == Face::ZHI)
//								     rhs(i,j,k,direction) = value;)
					}
				});
			}					
		}


		//#if AMREX_SPACEDIM > 0
		//Set(Face::XLO,Direction::X,m_bc_type[Face::XLO][Direction::X],m_bc_val[Face::XLO][Direction::X],a_rhs,a_geom);
		//Set(Face::XHI,Direction::X,m_bc_type[Face::XHI][Direction::X],m_bc_val[Face::XHI][Direction::X],a_rhs,a_geom);
		//#endif
		//#if AMREX_SPACEDIM > 1
		//Set(Face::XLO,Direction::Y,m_bc_type[Face::XLO][Direction::Y],m_bc_val[Face::XLO][Direction::Y],a_rhs,a_geom);
		//Set(Face::XHI,Direction::Y,m_bc_type[Face::XHI][Direction::Y],m_bc_val[Face::XHI][Direction::Y],a_rhs,a_geom);
		//Set(Face::YLO,Direction::X,m_bc_type[Face::YLO][Direction::X],m_bc_val[Face::YLO][Direction::X],a_rhs,a_geom);
		//Set(Face::YLO,Direction::Y,m_bc_type[Face::YLO][Direction::Y],m_bc_val[Face::YLO][Direction::Y],a_rhs,a_geom);
		//Set(Face::YHI,Direction::X,m_bc_type[Face::YHI][Direction::X],m_bc_val[Face::YHI][Direction::X],a_rhs,a_geom);
		//Set(Face::YHI,Direction::Y,m_bc_type[Face::YHI][Direction::Y],m_bc_val[Face::YHI][Direction::Y],a_rhs,a_geom);
		//#endif
		//#if AMREX_SPACEDIM > 2
		//Set(Face::XLO,Direction::Z,m_bc_type[Face::XLO][Direction::Z],m_bc_val[Face::XLO][Direction::Z],a_rhs,a_geom);
		//Set(Face::XHI,Direction::Z,m_bc_type[Face::XHI][Direction::Z],m_bc_val[Face::XHI][Direction::Z],a_rhs,a_geom);
		//Set(Face::YLO,Direction::Z,m_bc_type[Face::YLO][Direction::Z],m_bc_val[Face::YLO][Direction::Z],a_rhs,a_geom);
		//Set(Face::YHI,Direction::Z,m_bc_type[Face::YHI][Direction::Z],m_bc_val[Face::YHI][Direction::Z],a_rhs,a_geom);
		//Set(Face::ZLO,Direction::X,m_bc_type[Face::ZLO][Direction::X],m_bc_val[Face::ZLO][Direction::X],a_rhs,a_geom);
		//Set(Face::ZLO,Direction::Y,m_bc_type[Face::ZLO][Direction::Y],m_bc_val[Face::ZLO][Direction::Y],a_rhs,a_geom);
		//Set(Face::ZLO,Direction::Z,m_bc_type[Face::ZLO][Direction::Z],m_bc_val[Face::ZLO][Direction::Z],a_rhs,a_geom);
		//Set(Face::ZHI,Direction::X,m_bc_type[Face::ZHI][Direction::X],m_bc_val[Face::ZHI][Direction::X],a_rhs,a_geom);
		//Set(Face::ZHI,Direction::Y,m_bc_type[Face::ZHI][Direction::Y],m_bc_val[Face::ZHI][Direction::Y],a_rhs,a_geom);
		//Set(Face::ZHI,Direction::Z,m_bc_type[Face::ZHI][Direction::Z],m_bc_val[Face::ZHI][Direction::Z],a_rhs,a_geom);
		//#endif

		//TMP
		//Set(Face::XHI,Direction::Y,m_bc_hi[Face::XHI%AMREX_SPACEDIM][Direction::Y],m_bc_hi_val[Face::XHI%AMREX_SPACEDIM][Direction::Y],a_rhs,a_geom);

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

#define SQRT3INV 0.57735026919
#define SQRT2INV 0.70710678118


 	AMREX_FORCE_INLINE
 	Set::Vector operator () (const Set::Vector &u,
 				 const Set::Matrix &gradu,
 				 const Set::Matrix &sigma,
 				 const int &i, const int &j, const int &k,
				 const amrex::Box &domain)
 	{
		(void)i; (void)j; (void)k; // Suppress "unused variable" warnings
		//Set::Vector f;

		const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
		
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		 /*

		if (m[0] > domain.loVect()[0] && m[0] < domain.hiVect()[0])
		{
			if (m[1] == domain.hiVect()[1] && m[2] == domain.hiVect()[2])
			{
				Set::Vector n(0,0.70710678118,0.70710678118);
				f = sigma*n;
				return f;
			}

			if (m[1] == domain.loVect()[1] && m[2] == domain.loVect()[2])
			{
				Set::Vector n(0,-0.70710678118,-0.70710678118);
				f = sigma*n;
				return f;
			}

			if (m[1] == domain.loVect()[1] && m[2] == domain.hiVect()[2])
			{
				Set::Vector n(0,0.70710678118,-0.70710678118);
				f = sigma*n;
				return f;
			}

			if (m[1] == domain.hiVect()[1] && m[2] == domain.loVect()[2])
			{
				Set::Vector n(0,-0.70710678118,0.70710678118);
				f = sigma*n;
				return f;
			}
		}
		else if (m[0] == domain.loVect()[0] && m[1] == domain.loVect()[0])
		{
			return u;
		} 
		else if (m[0] == domain.hiVect()[0])
		{
			return u;
		}
		//if (m[0] == domain.loVect()[0] && m[1] == domain.loVect()[1])
		//	{
 		//		f(0) = u(0);
		//		f(1) = u(1);
		//		f(2) = u(2);
		//		return f;
		//	}
		*/

		// Corners
		if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YLO_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,-SQRT3INV));
		if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YLO_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,+SQRT3INV));
		if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YHI_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,-SQRT3INV));
		if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YHI_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,+SQRT3INV));
		if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YLO_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,-SQRT3INV));
		if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YLO_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,+SQRT3INV));
		if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YHI_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,-SQRT3INV));
		if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YHI_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,+SQRT3INV));
		
		if (m[1] == lo.y && m[2] == lo.z)                  return set(m_bc_type[Face::YLO_ZLO],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,-SQRT2INV));
		if (m[1] == lo.y && m[2] == hi.z)                  return set(m_bc_type[Face::YLO_ZHI],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,+SQRT2INV));
		if (m[1] == hi.y && m[2] == lo.z)                  return set(m_bc_type[Face::YHI_ZLO],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,-SQRT2INV));
		if (m[1] == hi.y && m[2] == hi.z)                  return set(m_bc_type[Face::YHI_ZHI],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,+SQRT2INV));
		if (m[2] == lo.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZLO_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       -SQRT2INV));
		if (m[2] == lo.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZLO_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       -SQRT2INV));
		if (m[2] == hi.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZHI_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       +SQRT2INV));
		if (m[2] == hi.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZHI_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       +SQRT2INV));
		if (m[0] == lo.x && m[1] == lo.y)                  return set(m_bc_type[Face::XLO_YLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV,0       ));
		if (m[0] == lo.x && m[1] == hi.y)                  return set(m_bc_type[Face::XLO_YHI],     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV,0       ));
		if (m[0] == hi.x && m[1] == lo.y)                  return set(m_bc_type[Face::XHI_YLO],     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV,0       ));
		if (m[0] == hi.x && m[1] == hi.y)                  return set(m_bc_type[Face::XHI_YHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV,0       ));
		
		if (m[0] == lo.x)                                  return set(m_bc_type[Face::XLO],         u, gradu, sigma, Set::Vector(-1, 0, 0));
		if (m[0] == hi.x)                                  return set(m_bc_type[Face::XHI],         u, gradu, sigma, Set::Vector(+1, 0, 0));
		if (m[1] == lo.y)                                  return set(m_bc_type[Face::YLO],         u, gradu, sigma, Set::Vector( 0,-1, 0));
		if (m[1] == hi.y)                                  return set(m_bc_type[Face::YHI],         u, gradu, sigma, Set::Vector( 0,+1, 0));
		if (m[2] == lo.z)                                  return set(m_bc_type[Face::ZLO],         u, gradu, sigma, Set::Vector( 0, 0,-1));
		if (m[2] == hi.z)                                  return set(m_bc_type[Face::ZHI],         u, gradu, sigma, Set::Vector( 0, 0,+1));
		
		Util::Abort(INFO,"Boundary condition error");
		return Set::Vector::Zero();
		/*
		for (int p = 0; p < AMREX_SPACEDIM; p++) // iterate over DIMENSIONS
		{
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
		}*/
		// bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
		// 	AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));
		
 		//return f;
 	}

	AMREX_FORCE_INLINE
	Set::Vector set(std::array<Type,3> &bc_type, 
					const Set::Vector &u, const Set::Matrix &gradu, const Set::Matrix &sigma, Set::Vector n) const
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
		{
			if      (bc_type[i] == Type::Displacement) 
				f(i) = u(i);
			else if (bc_type[i] == Type::Traction)
				f(i) = (sigma*n)(i);
			else if (bc_type[i] == Type::Neumann)
				f(i) = (gradu*n)(i);
			else if (bc_type[i] == Periodic)
				continue;
		}
		return f;
	}

private:
	// m_bc_lo[face][dimension]
	std::array<std::array<Type,       3>, 26> m_bc_type; 
	std::array<std::array<Set::Scalar,3>, 26> m_bc_val; 

	//Type xloylo = Type::Displacement;
	//Set::Scalar xloylo_val = 0.0;

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
		
		// TYPES
		
		if (pp.contains("type.xloylozlo")) { pp.queryarr("type.xloylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xloylozhi")) { pp.queryarr("type.xloylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.xloyhizlo")) { pp.queryarr("type.xloyhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xloyhizhi")) { pp.queryarr("type.xloyhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiylozlo")) { pp.queryarr("type.xhiylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiylozhi")) { pp.queryarr("type.xhiylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiyhizlo")) { pp.queryarr("type.xhiyhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiyhizhi")) { pp.queryarr("type.xhiyhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZHI][i] = bcmap[str[i]]; }

		if (pp.contains("type.ylozlo")) { pp.queryarr("type.ylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.ylozhi")) { pp.queryarr("type.ylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.yhizlo")) { pp.queryarr("type.yhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.yhizhi")) { pp.queryarr("type.yhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.zloxlo")) { pp.queryarr("type.zloxlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.zloxhi")) { pp.queryarr("type.zloxhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.zhixlo")) { pp.queryarr("type.zhixlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.zhixhi")) { pp.queryarr("type.zhixhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.xloylo")) { pp.queryarr("type.xloylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xloyhi")) { pp.queryarr("type.xloyhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiylo")) { pp.queryarr("type.xhiylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhiyhi")) { pp.queryarr("type.xhiyhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI][i] = bcmap[str[i]]; }

		if (pp.contains("type.xlo")) { pp.queryarr("type.xlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.xhi")) { pp.queryarr("type.xhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.ylo")) { pp.queryarr("type.ylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.yhi")) { pp.queryarr("type.yhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI][i] = bcmap[str[i]]; }
		if (pp.contains("type.zlo")) { pp.queryarr("type.zlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO][i] = bcmap[str[i]]; }
		if (pp.contains("type.zhi")) { pp.queryarr("type.zhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI][i] = bcmap[str[i]]; }


		// VALS
		std::vector<Set::Scalar> val;
		
		if (pp.contains("val.xloylozlo")) { pp.queryarr("val.xloylozlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO_ZLO][i] = val[i]; }
		if (pp.contains("val.xloylozhi")) { pp.queryarr("val.xloylozhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO_ZHI][i] = val[i]; }
		if (pp.contains("val.xloyhizlo")) { pp.queryarr("val.xloyhizlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI_ZLO][i] = val[i]; }
		if (pp.contains("val.xloyhizhi")) { pp.queryarr("val.xloyhizhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI_ZHI][i] = val[i]; }
		if (pp.contains("val.xhiylozlo")) { pp.queryarr("val.xhiylozlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO_ZLO][i] = val[i]; }
		if (pp.contains("val.xhiylozhi")) { pp.queryarr("val.xhiylozhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO_ZHI][i] = val[i]; }
		if (pp.contains("val.xhiyhizlo")) { pp.queryarr("val.xhiyhizlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI_ZLO][i] = val[i]; }
		if (pp.contains("val.xhiyhizhi")) { pp.queryarr("val.xhiyhizhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI_ZHI][i] = val[i]; }

		if (pp.contains("val.ylozlo")) { pp.queryarr("val.ylozlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO_ZLO][i] = val[i]; }
		if (pp.contains("val.ylozhi")) { pp.queryarr("val.ylozhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO_ZHI][i] = val[i]; }
		if (pp.contains("val.yhizlo")) { pp.queryarr("val.yhizlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI_ZLO][i] = val[i]; }
		if (pp.contains("val.yhizhi")) { pp.queryarr("val.yhizhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI_ZHI][i] = val[i]; }
		if (pp.contains("val.zloxlo")) { pp.queryarr("val.zloxlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO_XLO][i] = val[i]; }
		if (pp.contains("val.zloxhi")) { pp.queryarr("val.zloxhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO_XHI][i] = val[i]; }
		if (pp.contains("val.zhixlo")) { pp.queryarr("val.zhixlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI_XLO][i] = val[i]; }
		if (pp.contains("val.zhixhi")) { pp.queryarr("val.zhixhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI_XHI][i] = val[i]; }
		if (pp.contains("val.xloylo")) { pp.queryarr("val.xloylo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO][i] = val[i]; }
		if (pp.contains("val.xloyhi")) { pp.queryarr("val.xloyhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI][i] = val[i]; }
		if (pp.contains("val.xhiylo")) { pp.queryarr("val.xhiylo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO][i] = val[i]; }
		if (pp.contains("val.xhiyhi")) { pp.queryarr("val.xhiyhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI][i] = val[i]; }

		if (pp.contains("val.xlo")) { pp.queryarr("val.xlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO][i] = val[i]; }
		if (pp.contains("val.xhi")) { pp.queryarr("val.xhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI][i] = val[i]; }
		if (pp.contains("val.ylo")) { pp.queryarr("val.ylo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO][i] = val[i]; }
		if (pp.contains("val.yhi")) { pp.queryarr("val.yhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI][i] = val[i]; }
		if (pp.contains("val.zlo")) { pp.queryarr("val.zlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO][i] = val[i]; }
		if (pp.contains("val.zhi")) { pp.queryarr("val.zhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI][i] = val[i]; }

//		
//		if (pp.contains("val.xlo"))  { pp.queryarr("val.xlo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO][i] = val[i];}
//		if (pp.contains("val.xhi"))  { pp.queryarr("val.xhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI][i] = val[i];}
//		#if AMREX_SPACEDIM > 1
//		if (pp.contains("val.ylo"))  { pp.queryarr("val.ylo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO][i] = val[i];}
//		if (pp.contains("val.yhi"))  { pp.queryarr("val.yhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI][i] = val[i];}
//		#endif
//		#if AMREX_SPACEDIM > 2
//		if (pp.contains("val.zlo"))  { pp.queryarr("val.zlo", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO][i] = val[i];}
//		if (pp.contains("val.zhi"))  { pp.queryarr("val.zhi", val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI][i] = val[i];}
//		#endif
		//if (pp.contains("val.xlo")) for (int i = 0; i < AMREX_SPACEDIM; i++) 

    }
};

}
}
#endif
