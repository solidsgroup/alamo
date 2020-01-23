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
		for (int face = 0; face < AMREX_SPACEDIM; face++)
			for (int direction = 0; direction < AMREX_SPACEDIM; direction++)
				Util::Message(INFO,"face = ", face, " direction = ", direction, " val_lo = ", m_bc_lo_val[face][direction], " val_hi = ", m_bc_hi_val[face][direction]);

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
							
						for (int face = 0; face < 2*AMREX_SPACEDIM; face++)
						{
							for (int direction = 0; direction < AMREX_SPACEDIM; direction++)
							{
							
							
								AMREX_D_TERM(if (i == lo.x && face == Face::XLO)
										     rhs(i,j,k,direction) = m_bc_lo_val[face%AMREX_SPACEDIM][direction];,
									     if (j == lo.y && face == Face::YLO)
										     rhs(i,j,k,direction) = m_bc_lo_val[face%AMREX_SPACEDIM][direction];,
									     if (k == lo.z && face == Face::ZLO)
										     rhs(i,j,k,direction) = m_bc_lo_val[face%AMREX_SPACEDIM][direction];)
								AMREX_D_TERM(if (i == hi.x && face == Face::XHI)
										     rhs(i,j,k,direction) = m_bc_hi_val[face%AMREX_SPACEDIM][direction];,
									     if (j == hi.y && face == Face::YHI)
										     rhs(i,j,k,direction) = m_bc_hi_val[face%AMREX_SPACEDIM][direction];,
									     if (k == hi.z && face == Face::ZHI)
										     rhs(i,j,k,direction) = m_bc_hi_val[face%AMREX_SPACEDIM][direction];)

							}
						}
											 
											 
					});
			}					
		}
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

		// const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
		
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
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

		
//		std::array<Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xlo, bctype_ylo, bctype_zlo);
//		std::array<Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xhi, bctype_yhi, bctype_zhi);
//		
//		Set::Vector AMREX_D_DECL(bc_xlo,bc_ylo,bc_zlo);
//		Set::Vector AMREX_D_DECL(bc_xhi,bc_yhi,bc_zhi);
//
//		// Read in boundary types as strings, then convert to Operator::Elastic::BC types and store for future use.
//		std::vector<std::string> AMREX_D_DECL(bctype_xlo_str, bctype_ylo_str, bctype_zlo_str);
//		std::vector<std::string> AMREX_D_DECL(bctype_xhi_str, bctype_yhi_str, bctype_zhi_str);
//		AMREX_D_TERM(pp.queryarr("type_xlo", bctype_xlo_str);, pp.queryarr("type_ylo", bctype_ylo_str);, pp.queryarr("type_zlo", bctype_zlo_str););
//		AMREX_D_TERM(pp.queryarr("type_xhi", bctype_xhi_str);, pp.queryarr("type_yhi", bctype_yhi_str);, pp.queryarr("type_zhi", bctype_zhi_str););
//		if (AMREX_D_TERM(bctype_xlo_str.size() < AMREX_SPACEDIM, || bctype_ylo_str.size() < AMREX_SPACEDIM, || bctype_zlo_str.size() < AMREX_SPACEDIM) ||
//			AMREX_D_TERM(bctype_xhi_str.size() < AMREX_SPACEDIM, || bctype_yhi_str.size() < AMREX_SPACEDIM, || bctype_zhi_str.size() < AMREX_SPACEDIM))
////			Util::Abort(INFO, "incorrect number of terms specified in bctype");
//
//		AMREX_D_TERM(
//			bctype_xlo = {AMREX_D_DECL(bc[bctype_xlo_str[0]], bc[bctype_xlo_str[1]], bc[bctype_xlo_str[2]])};,
//			bctype_ylo = {AMREX_D_DECL(bc[bctype_ylo_str[0]], bc[bctype_ylo_str[1]], bc[bctype_ylo_str[2]])};,
//			bctype_zlo = {AMREX_D_DECL(bc[bctype_zlo_str[0]], bc[bctype_zlo_str[1]], bc[bctype_zlo_str[2]])};);
//		AMREX_D_TERM(
//			bctype_xhi = {AMREX_D_DECL(bc[bctype_xhi_str[0]], bc[bctype_xhi_str[1]], bc[bctype_xhi_str[2]])};,
//			bctype_yhi = {AMREX_D_DECL(bc[bctype_yhi_str[0]], bc[bctype_yhi_str[1]], bc[bctype_yhi_str[2]])};,
//			bctype_zhi = {AMREX_D_DECL(bc[bctype_zhi_str[0]], bc[bctype_zhi_str[1]], bc[bctype_zhi_str[2]])};);
//		AMREX_D_TERM(pp.queryarr("xlo", bc_xlo);, pp.queryarr("ylo", bc_ylo);, pp.queryarr("zlo", bc_zlo););
//		AMREX_D_TERM(pp.queryarr("xhi", bc_xhi);, pp.queryarr("yhi", bc_yhi);, pp.queryarr("zhi", bc_zhi););
		

		
		
//		value.Set(Face::XLO, Direction::X, bctype_xlo[Direction::X], bc_xlo[Direction::X]);
//		value.Set(Face::XHI, Direction::X, bctype_xhi[Direction::X], bc_xhi[Direction::X]);
//		#if AMREX_SPACEDIM > 1
//		value.Set(Face::XLO, Direction::Y, bctype_xlo[Direction::Y], bc_xlo[Direction::Y]);
//		value.Set(Face::XHI, Direction::Y, bctype_xhi[Direction::Y], bc_xhi[Direction::Y]);
//		value.Set(Face::YLO, Direction::X, bctype_ylo[Direction::X], bc_ylo[Direction::X]);
//		value.Set(Face::YLO, Direction::Y, bctype_ylo[Direction::Y], bc_ylo[Direction::Y]);
//		value.Set(Face::YHI, Direction::X, bctype_yhi[Direction::X], bc_yhi[Direction::X]);
//		value.Set(Face::YHI, Direction::Y, bctype_yhi[Direction::Y], bc_yhi[Direction::Y]);
//		#endif
//		#if AMREX_SPACEDIM > 2
//		value.Set(Face::XLO, Direction::Z, bctype_xlo[Direction::Z], bc_xlo[Direction::Z]);
//		value.Set(Face::XHI, Direction::Z, bctype_xhi[Direction::Z], bc_xhi[Direction::Z]);
//		value.Set(Face::YLO, Direction::Z, bctype_ylo[Direction::Z], bc_ylo[Direction::Z]);
//		value.Set(Face::YHI, Direction::Z, bctype_yhi[Direction::Z], bc_yhi[Direction::Z]);
//		value.Set(Face::ZLO, Direction::X, bctype_zlo[Direction::X], bc_zlo[Direction::X]);
//		value.Set(Face::ZLO, Direction::Y, bctype_zlo[Direction::Y], bc_zlo[Direction::Y]);
//		value.Set(Face::ZLO, Direction::Z, bctype_zlo[Direction::Z], bc_zlo[Direction::Z]);
//		value.Set(Face::ZHI, Direction::X, bctype_zhi[Direction::X], bc_zhi[Direction::X]);
//		value.Set(Face::ZHI, Direction::Y, bctype_zhi[Direction::Y], bc_zhi[Direction::Y]);
//		value.Set(Face::ZHI, Direction::Z, bctype_zhi[Direction::Z], bc_zhi[Direction::Z]);
//		#endif

    }
};

}
}
#endif
