#include "Elastic.H"

namespace BC
{

Elastic::Elastic(amrex::Vector<std::string> bc_hi_str
		    ,amrex::Vector<std::string> bc_lo_str
		    ,amrex::Vector<amrex::Real> _bc_lo_1
		    ,amrex::Vector<amrex::Real> _bc_hi_1
		    ,amrex::Vector<amrex::Real> _bc_lo_2
		    ,amrex::Vector<amrex::Real> _bc_hi_2
#if AMREX_SPACEDIM > 2
		    ,amrex::Vector<amrex::Real> _bc_lo_3
		    ,amrex::Vector<amrex::Real> _bc_hi_3
#endif
	  )
	: 
	bc_lo_1(_bc_lo_1), bc_hi_1(_bc_hi_1)
	,bc_lo_2(_bc_lo_2), bc_hi_2(_bc_hi_2)
#if AMREX_SPACEDIM > 2
	,bc_lo_3(_bc_lo_3), bc_hi_3(_bc_hi_3)
#endif
{
	for (int i=0;i<BL_SPACEDIM;i++)
	{
		bc_hi[i] = BCUtil::ReadString(bc_hi_str[i]);
		bc_lo[i] = BCUtil::ReadString(bc_lo_str[i]);

		if (BCUtil::IsPeriodic(bc_lo[i]) != BCUtil::IsPeriodic(bc_hi[i]))
			Util::Abort("Invalid BCs cannot be periodic on one side and not the other");
	}
}

void
Elastic::FillBoundary (amrex::MultiFab& mf, amrex::Real time)
{
	/* 
		We want to fill all dirichlet BC first and then Neumann.
		This is to ensure that when the stencil function is called, 
		cells that can be filled, have been filled.
		It is further assumed that ncomp = 3.
	*/
	
	amrex::Box domain(geom[lev].Domain());

	mf.FillBoundary(geom[lev].periodicity());

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
	{
		const amrex::Box& box = mfi.tilebox();

		amrex::BaseFab<amrex::Real> &mf_box = mf[mfi];

		if (BCUtil::IsDirichlet(bc_lo[0]))
		{
			int i = box.loVect()[0] - 1;
			AMREX_D_TERM(	,
		     			for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++),
					for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_lo_1[n];
			}
		}
		if (BCUtil::IsDirichlet(bc_hi[0]))
		{
			int i = box.hiVect()[0] + 1;
			AMREX_D_TERM(	,
		     			for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++),
					for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_hi_1[n];
			}
		}
		if (BCUtil::IsDirichlet(bc_lo[1]))
		{
			int j = box.loVect()[1] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++),
		     			,
					for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_lo_2[n];
			}
		}
		if (BCUtil::IsDirichlet(bc_hi[1]))
		{
			int j = box.hiVect()[1] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++),
		     			,
					for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_hi_2[n];
			}
		}
#if AMREX_SPACEDIM > 2
		if (BCUtil::IsDirichlet(bc_lo[2]))
		{
			int k = box.loVect()[2] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++),
		     			for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_lo_3[n];
			}
		}
		if (BCUtil::IsDirichlet(bc_hi[2]))
		{
			int k = box.hiVect()[2] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++),
		     			for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = bc_hi_3[n];
			}
		}
#endif
		/* Now that the dirichlet is taken care of, neumann bc should be implemented.
		Think of a 2D case - with a rectangular grid. Let's call the top-left 'real'
		cell as A. It will be surrounded by following ghost cells.
			1. Left ghost cell (say B)
			2. Top ghost cell (say C)
			3. Diagonal top left cell (say D)
			4. Diagonal bottom left (say E)
			5. Diagonal top right cell (say F)
		The next piece of code classifies these cells as following:
			1. End ghost cells: B and C
			2. Non-end ghost cells: E and F
			3. Corner ghost cell: D
		The next piece of code does the following:
			1. Fill non-end ghost cells next to all the faces
			2. Fill end, but non-corner ghost cell
			3. Fill corner ghost cell by averaging.
		*/
		if(BCUtil::IsNeumann(bc_lo[0]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int i = box.loVect()[0] - 1;
			AMREX_D_TERM(	,
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					for (int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				// the first entry of stencil is the center - so we have to shift right in this case
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+2*dx,0),mf_box(m+2*dx,1),mf_box(m+2*dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dx,0),mf_box(m-dy+dx,1),mf_box(m-dy+dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dx,0),mf_box(m+dy+dx,1),mf_box(m+dy+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dx,0),mf_box(m-dz+dx,1),mf_box(m-dz+dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dx,0),mf_box(m+dz+dx,1),mf_box(m+dz+dx,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));

				amrex::Vector<int> points;
				points.push_back(1);
				StencilFill(stencil, traction, points, m+dx, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[1](n);
			}
#if AMREX_SPACEDIM>1
			if(BCUtil::IsNeumann(bc_lo[1]))
			{	// This is the case when bottom boundary is Neumann. So stencil has two 
				// points missing.
				int j = box.loVect()[1];
				AMREX_D_TERM(	,
						,
						for(int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					amrex::Vector<Set::Vector> stencil;
					stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
					AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dx,0),mf_box(m+dx+dx,1),mf_box(m+dx+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dx,0),mf_box(m-dy+dx,1),mf_box(m-dy+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dx,0),mf_box(m+dy+dx,1),mf_box(m+dy+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dx,0),mf_box(m-dz+dx,1),mf_box(m-dz+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dx,0),mf_box(m+dz+dx,1),mf_box(m+dz+dx,2))));
						);
					amrex::Vector<Set::Vector> traction;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

					amrex::Vector<int> points;
					points.push_back(1);
					points.push_back(3);
					StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

					for (int n = 0; n<AMREX_SPACEDIM; n++)
					{
						mf_box(m,n) = stencil[1](n);
						mf_box(m+dx-dy,n) = stencil[3](n);
					}
				}
			}
			else
			{	// This is the case when bottom boundary is not Neumann. So stencil only has 
				// one point missing.
				int j = box.loVect()[1];
				AMREX_D_TERM(	,
						,
						for(int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					amrex::Vector<Set::Vector> stencil;
					stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
					AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dx,0),mf_box(m+dx+dx,1),mf_box(m+dx+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dx,0),mf_box(m-dy+dx,1),mf_box(m-dy+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dx,0),mf_box(m+dy+dx,1),mf_box(m+dy+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dx,0),mf_box(m-dz+dx,1),mf_box(m-dz+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dx,0),mf_box(m+dz+dx,1),mf_box(m+dz+dx,2))));
						);
					amrex::Vector<Set::Vector> traction;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
					
					amrex::Vector<int> points;
					points.push_back(1);
					
					StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

					for (int n = 0; n<AMREX_SPACEDIM; n++)
						mf_box(m,n) = stencil[1](n);
				}
			}
			if(BCUtil::IsNeumann(bc_hi[1]))
			{	// This is the case when top boundary is Neumann. So stencil has two 
				// points missing.
				int j = box.hiVect()[1];
				AMREX_D_TERM(	,
						,
						for(int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					amrex::Vector<Set::Vector> stencil;
					stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
					AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dx,0),mf_box(m+dx+dx,1),mf_box(m+dx+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dx,0),mf_box(m-dy+dx,1),mf_box(m-dy+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dx,0),mf_box(m+dy+dx,1),mf_box(m+dy+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dx,0),mf_box(m-dz+dx,1),mf_box(m-dz+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dx,0),mf_box(m+dz+dx,1),mf_box(m+dz+dx,2))));
						);
					amrex::Vector<Set::Vector> traction;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

					amrex::Vector<int> points;
					points.push_back(1);
					points.push_back(4);
					StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

					for (int n = 0; n<AMREX_SPACEDIM; n++)
					{
						mf_box(m,n) = stencil[1](n);
						mf_box(m+dx+dy,n) = stencil[4](n);
					}
				}
			}
			else
			{	// This is the case when bottom boundary is not Neumann. So stencil only has 
				// one point missing.
				int j = box.hiVect()[1];
				AMREX_D_TERM(	,
						,
						for(int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					amrex::Vector<Set::Vector> stencil;
					stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
					AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dx,0),mf_box(m+dx+dx,1),mf_box(m+dx+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dx,0),mf_box(m-dy+dx,1),mf_box(m-dy+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dx,0),mf_box(m+dy+dx,1),mf_box(m+dy+dx,2))));
							,
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dx,0),mf_box(m-dz+dx,1),mf_box(m-dz+dx,2))));
							stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dx,0),mf_box(m+dz+dx,1),mf_box(m+dz+dx,2))));
						);
					amrex::Vector<Set::Vector> traction;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
					
					amrex::Vector<int> points;
					points.push_back(1);
					
					StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

					for (int n = 0; n<AMREX_SPACEDIM; n++)
						mf_box(m,n) = stencil[1](n);
				}
			}
#endif
		}
		if(BCUtil::IsNeumann(bc_hi[0]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int i = box.hiVect()[0] + 1;
			AMREX_D_TERM(	,
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					for (int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));

				amrex::Vector<int> points;
				points.push_back(2);
				StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[2](n);
			}
		}
		if(BCUtil::IsNeumann(bc_lo[1]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int j = box.loVect()[1] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					,
					for (int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

				amrex::Vector<int> points;
				points.push_back(3);
				StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[3](n);
			}
		}
		if(BCUtil::IsNeumann(bc_hi[1]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int j = box.hiVect()[1] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					,
					for (int k = box.loVect()[2] + 1; k <= box.hiVect()[2] - 1; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

				amrex::Vector<int> points;
				points.push_back(4);
				StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[4](n);
			}
		}
#if AMREX_SPACEDIM > 2
		if(BCUtil::IsNeumann(bc_lo[2]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int k = box.loVect()[2] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

				amrex::Vector<int> points;
				points.push_back(5);
				StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[5](n);
			}
		}
		if(BCUtil::IsNeumann(bc_hi[2]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int k = box.hiVect()[2] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m,0),mf_box(m,1),mf_box(m,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

				amrex::Vector<int> points;
				points.push_back(6);
				StencilFill(stencil, traction, points, m, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[6](n);
			}
		}
#endif
		/* Non end ghost cells have been taken care of. Now we solve for end ghost cells */
	}

}

amrex::BCRec
Elastic::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
Elastic::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(bc_lo[0]),BCUtil::IsPeriodic(bc_lo[1]),BCUtil::IsPeriodic(bc_lo[2]))};
}

}
