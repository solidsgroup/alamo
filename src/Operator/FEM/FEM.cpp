#include <AMReX_MultiFabUtil.H>
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

#include "FEM.H"

namespace Operator
{
namespace FEM
{

FEM::FEM (Model::Solid::Solid &_model)
  : model(_model)
{
}

FEM::~FEM ()
{}

template<> 
void FEM::Energy<1> (int /*amrlev*/,		///<[in] AMR Level
		     int /*mglev*/,		///<[in]
		     MultiFab& /*dw*/,	///<[out] The force vector
		     const MultiFab& /*u*/	///<[in] The displacements vector
		     ) const
{
  BL_PROFILE("Operator::FEM::FEM::apply()");

  // const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  // Element<Q4> element(DX[0], DX[1]);

  // // DPhi[shape function #][quadrature point #][dimension #]
  // std::array<std::array<std::array<amrex::Real,2>,4>,4> DPhi;// = element.DPhis();
  // std::array<amrex::Real,4> W = element.Ws();


  // static amrex::IntVect dx(1,0), dy(0,1);
  // for (MFIter mfi(dw, true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();
  //     const FArrayBox &ufab  = u[mfi];
  //     FArrayBox       &dwfab  = dw[mfi];

  //     dwfab.setVal(0.0);

  //     for (int my = bx.loVect()[1]-1; my<=bx.hiVect()[1]; my++)
  // 	for (int mx = bx.loVect()[0]-1; mx<=bx.hiVect()[0]; mx++)
  // 	  {
  // 	    amrex::IntVect O(mx,my);
  // 	    for (int _m=0; _m<4; _m++)
  // 	      {
  // 		amrex::IntVect m = O + dx*(_m%2) + dy*((_m/2)%2); // dz*((_n/4)%2)

  // 		if (m[0] < bx.loVect()[0] || m[0] > bx.hiVect()[0] ||
  // 		    m[1] < bx.loVect()[1] || m[1] > bx.hiVect()[1]) continue;

  // 		for (int i=0; i < 2; i++)
  // 		  {
  // 		    amrex::Real K_minj = 0.0;
  // 		    //K[_m][i][_n][j] = 0.0;

  // 		    for (int p=0; p < 2; p++)
  // 		      for (int Q=0; Q < 4; Q++)
  // 			{
  // 			  K_minj +=
  // 			    W[Q] *
  // 			    //C(i,p,j,q,m,amrlev,mglev,mfi) * 
  // 			    DPhi[_m][Q][p];// *
  // 			    //DPhi[_n][Q][q];
  // 			}
  // 		    //dwfab(m,i) -= K_minj * ufab(n,j);

  // 		  }
		
  // 	      }
  // 	  }
  //   }
}



template<>
void FEM::Energy<2> (int amrlev,		///<[in] AMR Level
		     int mglev,		///<[in]
		     MultiFab& f,	///<[out] The force vector
		     const MultiFab& u	///<[in] The displacements vector
		     ) const
{
  BL_PROFILE("Operator::FEM::FEM::Energy<2>()");

  const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  Element<Q4> element(DX[0], DX[1]);

  // DPhi[shape function #][quadrature point #][dimension #]
  std::array<std::array<Set::Vector,4>,4> DPhi = element.DPhis();
  std::array<amrex::Real,4> W = element.Ws();

  // Grid dimensions
  static amrex::IntVect dx(1,0), dy(0,1);
  // Iterating over patches in a single AMR level
  //   *** mfi = "multi fab iterator"
  for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
      // Variable for each patch
      const Box& bx = mfi.tilebox();
      const FArrayBox &ufab  = u[mfi];
      FArrayBox       &ffab  = f[mfi];

      // Initialize the result to zero
      ffab.setVal(0.0);

      // Iterate over the grid
      // - INCLUDE x=0,y=0,(z=0) ghost cells
      // - DON'T INCLUDE x=xmax, y=ymax, (z=zmax) ghost cells
      // actually iterating over CELLS
      for (int my = bx.loVect()[1]-1; my<=bx.hiVect()[1]; my++)
	for (int mx = bx.loVect()[0]-1; mx<=bx.hiVect()[0]; mx++)
	  {
	    // O is the origin for the cell that we're on
	    amrex::IntVect O(mx,my);

	    // Iterate over every quadrature point
	    for (int Q=0; Q < 4; Q++)
	      {
		// Compute gradient of u at this quadrature point
		// using shape functions to interpolate
		Set::Matrix gradu = Set::Matrix::Zero();
		for (int i=0; i<2; i++)
		  for (int j=0; i<2; i++)
		    gradu(i,j) +=
		      ufab(O,i) * DPhi[0][Q][j] +
		      ufab(O+dx,i) * DPhi[1][Q][j] + 
		      ufab(O+dy,i) * DPhi[2][Q][j] + 
		      ufab(O+dx+dy,i) * DPhi[3][Q][j];

		// Iterate over every node in the cell
		for (int _m=0; _m<4; _m++)
		  {
		    // Computes the IntVect for the node we're on
		    amrex::IntVect m = O + dx*(_m%2) + dy*((_m/2)%2); 

		    // Safety: abort if we're outside the domain
		    // because m indexes the node whose value we're setting
		    // and the output FAB does not have ghost cells
		    if (m[0] < bx.loVect()[0] || m[0] > bx.hiVect()[0] ||
			m[1] < bx.loVect()[1] || m[1] > bx.hiVect()[1]) continue;

		    // Iterate over every node in the cell
		    for (int _n=0; _n<4; _n++)
		      {
			// Compute the IntVect for node n
			amrex::IntVect n = O + dx*(_n%2) + dy*((_n/2)%2);

			// Compute the quantity
			// ddw(i,j) = [quadrature weight]  *  DDW(i,p,j,q)*DPhi(m,p)*DPhi(n,q)
			Set::Matrix ddw = W[Q]*model.DDW(gradu,DPhi[_m][Q],DPhi[_n][Q]);

			// Matrix multiplcation
			// f(m,i) = ddw(i,j) * ufab(n,j)
			// NOTE: here we are actually APPLYING the stiffness matrix to
			// the input, ufab.
			// (as opposed to some FEM implementation where you calculate
			// the stiffness matrix only)
			for (int i=0; i < 2; i++)
			  for (int j=0; j < 2; j++)
			    ffab(m,i) += ddw(i,j) * ufab(n,j);
		      }
		  }
	      }
	  }
	      
    }
}
void FEM::smooth (int amrlev,          ///<[in] AMR level
		  int mglev,           ///<[in]
		  MultiFab& u,         ///<[inout] Solution (displacement field)
		  const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
		  int redblack         ///<[in] Smooth even vs. odd modes
		  ) const
{
  BL_PROFILE("Operator::FEM::FEM::smooth()");
  return;

  amrex::Abort("Does not work if using smooth");

  const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  Element<Q4> element(DX[0], DX[1]);

  // // DPhi[shape function #][quadrature point #][dimension #]
  std::array<std::array<Set::Vector,4>,4> DPhi = element.DPhis();
  std::array<amrex::Real,4> W = element.Ws();

  MultiFab rho, aa;
  rho.define(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], getNComp(),0);
  aa.define(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], getNComp(),0);

  static amrex::IntVect dx(1,0), dy(0,1);
  for (MFIter mfi(u, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const FArrayBox &rhsfab  = rhs[mfi];
      FArrayBox       &ufab  = u[mfi];
      FArrayBox &rhofab = rho[mfi];
      FArrayBox &aafab = aa[mfi];

      aafab.setVal(0.0);
      rhofab.setVal(0.0);

      for (int mx = bx.loVect()[0]-1; mx<=bx.hiVect()[0]; mx++)
	for (int my = bx.loVect()[1]-1; my<=bx.hiVect()[1]; my++)
	  {
	    amrex::IntVect O(mx,my);

	    for (int Q=0; Q < 4; Q++)
	      {
		Set::Matrix gradu = Set::Matrix::Zero();
		for (int i=0; i<2; i++)
		  for (int j=0; i<2; i++)
		    gradu(i,j) +=
		      ufab(O,i) * DPhi[0][Q][j] +
		      ufab(O+dx,i) * DPhi[1][Q][j] + 
		      ufab(O+dy,i) * DPhi[2][Q][j] + 
		      ufab(O+dx+dy,i) * DPhi[3][Q][j];

		for (int _m=0; _m<4; _m++)
		  {
		    amrex::IntVect m = O + dx*((_m/1)%2) + dy*((_m/2)%2); // dz*((_m/4)%2)
		    if (m[0] < bx.loVect()[0] || m[0] > bx.hiVect()[0] ||
			m[1] < bx.loVect()[1] || m[1] > bx.hiVect()[1]) continue;

		    if ( (m[0] + m[1])%2 == redblack) continue; // Red-Black bit


		    for (int _n=0; _n<4; _n++)
		      {
			amrex::IntVect n = O + dx*(_n%2) + dy*((_n/2)%2); // dz*((_n/4)%2)

			Set::Matrix ddw = W[Q]*model.DDW(gradu,DPhi[_m][Q],DPhi[_n][Q]);

			for (int i=0; i < 2; i++)
			  for (int j=0; j < 2; j++)
			    if (m==n && i==j)
			      aafab(m,i) = ddw(i,j);
			    else
			      rhofab(m,i) += ddw(i,j) * ufab(n,j);
		      }
		  }
	      }
	  }
      for (int mx = bx.loVect()[0]; mx<=bx.hiVect()[0]; mx++)
	for (int my = bx.loVect()[1]; my<=bx.hiVect()[1]; my++)
   	  {
   	    amrex::IntVect m(mx,my);
	    if ( (m[0] + m[1])%2 == redblack) continue;
	    for (int i = 0; i<AMREX_SPACEDIM; i++)
	      ufab(m,i) = (rhsfab(m,i) - rhofab(m,i))/aafab(m,i);
	  }
    }
}


void
FEM::Fapply (int amrlev, ///<[in] AMR Level
	     int mglev,  ///<[in]
	     MultiFab& f,///<[out] The force vector
	     const MultiFab& u ///<[in] The displacements vector
	     ) const
{
  BL_PROFILE("Operator::FEM::FEM::Fapply()");
  Energy<2>(amrlev,mglev,f,u);
}


/// \fn Operator::FEM::Fsmooth
///
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in Operator::FEM::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
///
/// \todo Extend to 3D
///
void
FEM::Fsmooth (int amrlev,          ///<[in] AMR level
			     int mglev,           ///<[in]
			     MultiFab& u,         ///<[inout] Solution (displacement field)
			     const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			     int redblack         ///<[in] Smooth even vs. odd modes
			     ) const
{
  BL_PROFILE("Operator::FEM::FEM::Fsmooth()");
  //std::cout << "IN FSMOOTH"<< std::endl;
  smooth(amrlev,mglev,u,rhs,redblack);
}

/// \fn Operator::FEM::FFlux
///
/// The fluxes are simply set to zero and returned.
///
/// \todo Extend to 3D
///
void
FEM::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
			   const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
			   const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
  // THIS ONLY HAPPENS WHEN MULTIPLE AMR LEVELS ARE USED...

  amrex::BaseFab<amrex::Real> &fxfab = *sigmafab[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigmafab[1];
  fxfab.setVal(0.0);
  fyfab.setVal(0.0);
}


void
FEM::Stress (FArrayBox& sigmafab,
			    const FArrayBox& ufab,
			    int amrlev, const MFIter& mfi) const
{
  /// \todo add assert for sigmafab.ncomp=3 (for SPACEDIM=2) and ncomp=6 (for SPACEDIM=3)

  const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  amrex::IntVect dx(1,0);
  amrex::IntVect dy(0,1);

  const Box& bx = mfi.tilebox();
  for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
      {
	amrex::IntVect m(m1,m2);

	amrex::Real du1_dx1 = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
	amrex::Real du1_dx2 = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
	amrex::Real du2_dx1 = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
	amrex::Real du2_dx2 = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);

	for (int i=0; i<AMREX_SPACEDIM; i++)
	  for (int j=0; j<AMREX_SPACEDIM; j++)
	    {
	      int voigt;
	      if (i==0 && j==0) voigt = 0;
	      if (i==1 && j==1) voigt = 1;
	      if (i==0 && j==1) voigt = 2;
	      if (i==1 && j==0) continue;

	      sigmafab(m,voigt) =
		C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 +
		C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 +
		C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 +
		C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;
	    }
      }

}

template<>
void FEM::Energy<0> (int /*amrlev*/,
		     int /*mglev*/,
		     MultiFab & /*w*/,
		     const MultiFab &/*u*/) const
// FArrayBox& energyfab,
// 			    const FArrayBox& ufab,
// 			    int amrlev, const MFIter& mfi) const
{
  // const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  // amrex::IntVect dx(1,0);
  // amrex::IntVect dy(0,1);

  // const Box& bx = mfi.tilebox();
  // for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
  //   for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
  //     {
  // 	amrex::IntVect m(m1,m2);
  // 	amrex::Real gradu[2][2] = {{(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]),
  // 				    (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1])},
  // 				   {(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]),
  // 				    (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1])}};

  // 	energyfab(m) = 0.0;

  // 	for (int i=0; i<AMREX_SPACEDIM; i++)
  // 	  for (int j=0; j<AMREX_SPACEDIM; j++)
  // 	    for (int k=0; k<AMREX_SPACEDIM; k++)
  // 	      for (int l=0; l<AMREX_SPACEDIM; l++)
  // 		energyfab(m) += gradu[i][j] * C(i,j,k,l,m,amrlev,0,mfi) * gradu[k][l];
  // 	energyfab(m) *= 0.5;
  //     }

}

amrex::Real
FEM::C(const int i, const int j, const int k, const int l,
		      const amrex::IntVect /*loc*/, const int /*amrlev*/,const int /*mglev*/, const MFIter & /*mfi*/) const
{
  amrex::Real mu = 1.0, lambda=1.0;
  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  return ret;
}
}
}
