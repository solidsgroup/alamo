#ifndef BC_OPERATOR_ELASTIC_H
#define BC_OPERATOR_ELASTIC_H

// #include "Operator/Elastic.H"

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
				m_bc_hi[face][direction] = Type::Displacement;
			}
	};
 	~Elastic() {};
	
	void
	Set(const Face face,
	      const Direction direction,
	      const Type type,
	      const Set::Scalar value,
	      amrex::Vector<amrex::MultiFab *> &a_rhs,
	    const amrex::Vector<amrex::Geometry> &a_geom)
	{
		if (AMREX_D_TERM(face == Face::XLO, || face == Face::YLO, ||  face == Face::ZLO))
			m_bc_lo[face%AMREX_SPACEDIM][direction] = type;
		if (AMREX_D_TERM(face == Face::XHI, || face == Face::YHI, ||  face == Face::ZHI))
			m_bc_hi[face%AMREX_SPACEDIM][direction] = type;

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
	std::array<std::array<Type,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_lo; 
	// m_bc_hi[face][dimension]
	std::array<std::array<Type,AMREX_SPACEDIM>, AMREX_SPACEDIM> m_bc_hi; 
};

}
}
#endif
