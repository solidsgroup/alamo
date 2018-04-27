#include <AMReX_MultiFabUtil.H>
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

#include "FEM.H"

/// \fn Operator::FEM::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
Operator::FEM::FEM::FEM (const Vector<Geometry>& a_geom,
			 const Vector<BoxArray>& a_grids,
			 const Vector<DistributionMapping>& a_dmap,
			 const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);

  // rho.resize(m_num_amr_levels);
  // aa.resize(m_num_amr_levels);
  // for (int amrlev=0; amrlev < m_num_amr_levels; amrlev++)
  //   {
  //     rho[amrlev].resize(m_num_mg_levels[amrlev]);
  //     aa[amrlev].resize(m_num_mg_levels[amrlev]);
  //     for (int mglev=0; mglev < m_num_mg_levels[amrlev]; mglev++)
  // 	{
  // 	  rho[amrlev][mglev].define(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], getNComp(),0);
  // 	  aa[amrlev][mglev].define(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], getNComp(),0);
  // 	}
  //   }
}

Operator::FEM::FEM::~FEM ()
{}



void
Operator::FEM::FEM::apply (int amrlev, ///<[in] AMR Level
			    int mglev,  ///<[in]
			    MultiFab& f,///<[out] The force vector
			    const MultiFab& u ///<[in] The displacements vector
			    ) const
{
  BL_PROFILE("Operator::FEM::FEM::apply()");

  const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  Element<Q4> element(DX[0], DX[1]);

  // DPhi[shape function #][quadrature point #][dimension #]
  std::array<std::array<std::array<amrex::Real,2>,4>,4> DPhi = 
    {{{element.DPhi<1>(element.QPoint<1>()),
       element.DPhi<1>(element.QPoint<2>()),
       element.DPhi<1>(element.QPoint<3>()),
       element.DPhi<1>(element.QPoint<4>())},
      {element.DPhi<2>(element.QPoint<1>()),
       element.DPhi<2>(element.QPoint<2>()),
       element.DPhi<2>(element.QPoint<3>()),
       element.DPhi<2>(element.QPoint<4>())},
      {element.DPhi<3>(element.QPoint<1>()),
       element.DPhi<3>(element.QPoint<2>()),
       element.DPhi<3>(element.QPoint<3>()),
       element.DPhi<3>(element.QPoint<4>())},
      {element.DPhi<4>(element.QPoint<1>()),
       element.DPhi<4>(element.QPoint<2>()),
       element.DPhi<4>(element.QPoint<3>()),
       element.DPhi<4>(element.QPoint<4>())}}};
  std::array<amrex::Real,4> W =
    {element.QWeight<1>(), 
     element.QWeight<2>(), 
     element.QWeight<3>(), 
     element.QWeight<4>()};


  static amrex::IntVect dx(1,0), dy(0,1);
  for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const FArrayBox &ufab  = u[mfi];
      FArrayBox       &ffab  = f[mfi];

      ffab.setVal(0.0);

      for (int mx = bx.loVect()[0]-1; mx<=bx.hiVect()[0]; mx++)
	for (int my = bx.loVect()[1]-1; my<=bx.hiVect()[1]; my++)
	  {
	    amrex::IntVect O(mx,my);

	    for (int _m=0; _m<4; _m++)
	      {
		amrex::IntVect m = O + dx*(_m%2) + dy*((_m/2)%2); // dz*((_n/4)%2)

		if (m[0] < bx.loVect()[0] || m[0] > bx.hiVect()[0] ||
		    m[1] < bx.loVect()[1] || m[1] > bx.hiVect()[1]) continue;

		for (int _n=0; _n<4; _n++)
		  {
		    amrex::IntVect n = O + dx*(_n%2) + dy*((_n/2)%2); // dz*((_n/4)%2)

		    //if (fabs(ufab(n,0)) > 0) std::cout << ufab(n,0) << std::endl;

		    for (int i=0; i < 2; i++)
		      for (int j=0; j < 2; j++)
			for (int p=0; p < 2; p++)
			  for (int q=0; q < 2; q++)
			    for (int Q=0; Q<4; Q++)
			      {			    
				ffab(m,i) -=
				  W[Q] *
				  C(i,p,j,q,m,amrlev,mglev,mfi) * //TODO need averaged C
				  DPhi[_m][Q][p] *
				  DPhi[_n][Q][q] *
				  ufab(n,j);
			      }
		    // if (fabs(ffab(m,0)) > 0 || fabs(ffab(m,1)) > 0)
		    //   std::cout << "Fapply, ffab:" << ffab(m,0)<< "," << ffab(m,1) << std::endl;
		  }
	      }
	  }
    }
}

void
Operator::FEM::FEM::smooth (int amrlev,          ///<[in] AMR level
			     int mglev,           ///<[in]
			     MultiFab& u,         ///<[inout] Solution (displacement field)
			     const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			     int redblack         ///<[in] Smooth even vs. odd modes
			     ) const
{
  BL_PROFILE("Operator::FEM::FEM::smooth()");

  //if (!redblack) return;

  const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  Element<Q4> element(DX[0], DX[1]);

  // DPhi[shape function #][quadrature point #][dimension #]
  std::array<std::array<std::array<amrex::Real,2>,4>,4> DPhi = 
    {{{element.DPhi<1>(element.QPoint<1>()),
       element.DPhi<1>(element.QPoint<2>()),
       element.DPhi<1>(element.QPoint<3>()),
       element.DPhi<1>(element.QPoint<4>())},
      {element.DPhi<2>(element.QPoint<1>()),
       element.DPhi<2>(element.QPoint<2>()),
       element.DPhi<2>(element.QPoint<3>()),
       element.DPhi<2>(element.QPoint<4>())},
      {element.DPhi<3>(element.QPoint<1>()),
       element.DPhi<3>(element.QPoint<2>()),
       element.DPhi<3>(element.QPoint<3>()),
       element.DPhi<3>(element.QPoint<4>())},
      {element.DPhi<4>(element.QPoint<1>()),
       element.DPhi<4>(element.QPoint<2>()),
       element.DPhi<4>(element.QPoint<3>()),
       element.DPhi<4>(element.QPoint<4>())}}};
  std::array<amrex::Real,4> W =
    {element.QWeight<1>(), 
     element.QWeight<2>(), 
     element.QWeight<3>(), 
     element.QWeight<4>()};


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

	    for (int _m=0; _m<4; _m++)
	      {
		amrex::IntVect m = O + dx*(_m%2) + dy*((_m/2)%2); // dz*((_n/4)%2)
		if (m[0] < bx.loVect()[0] || m[0] > bx.hiVect()[0] ||
		    m[1] < bx.loVect()[1] || m[1] > bx.hiVect()[1]) continue;


		if ( (m[0] + m[1])%2 == redblack) continue; // Red-Black bit

		for (int _n=0; _n<4; _n++)
		  {
		    amrex::IntVect n = O + dx*(_n%2) + dy*((_n/2)%2); // dz*((_n/4)%2)

		    // if (n[0] < bx.loVect()[0] || n[0] > bx.hiVect()[0] ||
		    // 	n[1] < bx.loVect()[1] || n[1] > bx.hiVect()[1]) continue;

		    for (int i=0; i < 2; i++)
		      for (int j=0; j < 2; j++)
			{
			  for (int p=0; p < 2; p++)
			    for (int q=0; q < 2; q++)
			      for (int Q=0; Q < 4; Q++)
				{			    
				  if (m==n && i==j)
				    aafab(m,i) +=
				      -W[Q] *
				      C(i,p,j,q,m,amrlev,mglev,mfi) * //TODO need averaged C  
				      DPhi[_m][Q][p] *
				      DPhi[_n][Q][q];
				  else 
				    rhofab(m,i) +=
				      -W[Q] *
				      C(i,p,j,q,m,amrlev,mglev,mfi) * //TODO need averaged C  
				      DPhi[_m][Q][p] *
				      DPhi[_n][Q][q] *
				      ufab(n,j);
				}
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
  // MultiFab::Copy(u, rhs, 0, 0, u.nComp(), 0);
  // u.minus(rho,0,u.nComp(),0);
  // u.divide(aa,0,u.nComp(),0);

}


void
Operator::FEM::FEM::Fapply (int amrlev, ///<[in] AMR Level
			    int mglev,  ///<[in]
			    MultiFab& f,///<[out] The force vector
			    const MultiFab& u ///<[in] The displacements vector
			    ) const
{
  BL_PROFILE("Operator::FEM::FEM::Fapply()");
  apply(amrlev,mglev,f,u);
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
Operator::FEM::FEM::Fsmooth (int amrlev,          ///<[in] AMR level
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
Operator::FEM::FEM::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
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
Operator::FEM::FEM::Stress (FArrayBox& sigmafab,
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

void
Operator::FEM::FEM::Energy (FArrayBox& energyfab,
			    const FArrayBox& ufab,
			    int amrlev, const MFIter& mfi) const
{
  const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  amrex::IntVect dx(1,0);
  amrex::IntVect dy(0,1);

  const Box& bx = mfi.tilebox();
  for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
      {
	amrex::IntVect m(m1,m2);
	amrex::Real gradu[2][2] = {{(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]),
				    (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1])},
				   {(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]),
				    (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1])}};

	energyfab(m) = 0.0;

	for (int i=0; i<AMREX_SPACEDIM; i++)
	  for (int j=0; j<AMREX_SPACEDIM; j++)
	    for (int k=0; k<AMREX_SPACEDIM; k++)
	      for (int l=0; l<AMREX_SPACEDIM; l++)
		energyfab(m) += gradu[i][j] * C(i,j,k,l,m,amrlev,0,mfi) * gradu[k][l];
	energyfab(m) *= 0.5;
      }

}

amrex::Real
Operator::FEM::FEM::C(const int i, const int j, const int k, const int l,
				const amrex::IntVect /*loc*/, const int /*amrlev*/,const int /*mglev*/, const MFIter & /*mfi*/) const
{
  amrex::Real mu = 1.0, lambda=1.0;
  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  return ret;
}
