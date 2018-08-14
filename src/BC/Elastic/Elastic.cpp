#include "Elastic.H"

namespace BC
{

Elastic::Elastic(amrex::Vector<std::string> bc_hi_str,
		    amrex::Vector<std::string> bc_lo_str,
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
				 amrex::Vector<amrex::Real> _bc_lo_2,
				 amrex::Vector<amrex::Real> _bc_lo_3),
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
				 amrex::Vector<amrex::Real> _bc_hi_2,
				 amrex::Vector<amrex::Real> _bc_hi_3))
	: 
	AMREX_D_DECL(bc_lo_1(_bc_lo_1),bc_lo_2(_bc_lo_2),bc_lo_3(_bc_lo_3)),
	AMREX_D_DECL(bc_hi_1(_bc_hi_1),bc_hi_2(_bc_hi_2),bc_hi_3(_bc_hi_3))
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
	
	amrex::Box domain(m_geom.Domain());

	mf.FillBoundary(m_geom.periodicity());

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
	{
		const amrex::Box& box = mfi.tilebox();

		amrex::BaseFab<amrex::Real> &mf_box = mf[mfi];

		/* Dirichlet boundaries are first */
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
#if AMREX_SPACEDIM > 1
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

		/* Step 1: Fill the left face ghost cells */
		if (BCUtil::IsNeumann(bc_lo[0]))
		{
			int i = box.loVect()[0] - 1;

			AMREX_D_TERM(	,
					for (int j = box.loVect()[1]; j <= box.hiVect()[1]; j++),
					for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));

				amrex::Vector<Set::Vector> stencil;
				// the first entry of stencil is the center - so we have to shift right in this case
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx,0),mf_box(m+dx,1),mf_box(m+dx,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx+dx,0),mf_box(m-dx+dx,1),mf_box(m-dx+dx,2))));
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
#if AMREX_SPACEDIM > 1
				if (j == box.loVect()[1] && BCUtil::IsNeumann(bc_lo[1]))
				{
					/* Solving for end points */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
					points.push_back(3);
#if AMREX_SPACEDIM > 2
					/* Solving for triple end point */ 
					if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
						points.push_back(5);
					}
					else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
						points.push_back(6);
					}
#endif
				}
				else if (j == box.hiVect()[1] && BCUtil::IsNeumann(bc_hi[1]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
					points.push_back(4);
#if AMREX_SPACEDIM > 2
					if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
						points.push_back(5);
					}
					else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
						points.push_back(6);
					}
#endif
				}
#if AMREX_SPACEDIM > 2
				else if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
				{
					/* Solving for end point. No need to consider triple end point.
					   They have already been solved. */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
#endif
				StencilFill(stencil, traction, points, m+dx, m_amrlev, m_mglev, mfi);

				for (int n = 0; n < AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[1](n);
#if AMREX_SPACEDIM > 1
				for (int p = 1; p < points.size(); p++)
				{
					switch (points[p])
					{
						case 3: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dy+dx,n) = stencil[3](n);
							break;
						case 4: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dy+dx,n) = stencil[4](n);
							break;
#if AMREX_SPACEDIM > 2
						case 5: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dz+dx,n) = stencil[5](n);
							break;
						case 6: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dz+dx,n) = stencil[6](n);
							break;
#endif
					}
				}
#endif
			}
		}

		
		/* Step 2: Fill right face of ghost cells */
		if(BCUtil::IsNeumann(bc_hi[0]))
		{
			int i = box.hiVect()[0] + 1;
			AMREX_D_TERM(	,
					for (int j = box.loVect()[1]; j <= box.hiVect()[1]; j++),
					for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx,0),mf_box(m-dx,1),mf_box(m-dx,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx-dx,0),mf_box(m-dx-dx,1),mf_box(m-dx-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx-dx,0),mf_box(m+dx-dx,1),mf_box(m+dx-dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy-dx,0),mf_box(m-dy-dx,1),mf_box(m-dy-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy-dx,0),mf_box(m+dy-dx,1),mf_box(m+dy-dx,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz-dx,0),mf_box(m-dz-dx,1),mf_box(m-dz-dx,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz-dx,0),mf_box(m+dz-dx,1),mf_box(m+dz-dx,2))));
					);

				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));

				amrex::Vector<int> points;
				points.push_back(2);

#if AMREX_SPACEDIM > 1
				if (j == box.loVect()[1] && BCUtil::IsNeumann(bc_lo[1]))
				{
					/* Solving for end points */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
					points.push_back(3);
#if AMREX_SPACEDIM > 2
					/* Solving for triple end point */ 
					if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
						points.push_back(5);
					}
					else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
						points.push_back(6);
					}
#endif
				}
				else if (j == box.hiVect()[1] && BCUtil::IsNeumann(bc_hi[1]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
					points.push_back(4);
#if AMREX_SPACEDIM > 2
					if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
						points.push_back(5);
					}
					else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
					{
						traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
						points.push_back(6);
					}
#endif
				}
#if AMREX_SPACEDIM > 2
				else if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
				{
					/* Solving for end point. No need to consider triple end point.
					   They have already been solved. */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
#endif

				StencilFill(stencil, traction, points, m-dx, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[2](n);
#if AMREX_SPACEDIM > 1
				for (int p = 1; p < points.size(); p++)
				{
					switch (points[p])
					{
						case 3: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dy-dx,n) = stencil[3](n);
							break;
						case 4: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dy-dx,n) = stencil[4](n);
							break;
#if AMREX_SPACEDIM > 2
						case 5: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dz-dx,n) = stencil[5](n);
							break;
						case 6: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dz-dx,n) = stencil[6](n);
							break;
#endif
					}
				}
#endif
			}
		}
#if AMREX_SPACEDIM > 1
		/* Step 3: Fill bottom face of ghost cells */
		if(BCUtil::IsNeumann(bc_lo[1]))
		{
			int j = box.loVect()[1] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0]+1; i <= box.hiVect()[0]-1; i++),
					,
					for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy,0),mf_box(m+dy,1),mf_box(m+dy,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx+dy,0),mf_box(m-dx+dy,1),mf_box(m-dx+dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dy,0),mf_box(m+dx+dy,1),mf_box(m+dx+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dy,0),mf_box(m-dy+dy,1),mf_box(m-dy+dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dy,0),mf_box(m+dy+dy,1),mf_box(m+dy+dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dy,0),mf_box(m-dz+dy,1),mf_box(m-dz+dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dy,0),mf_box(m+dz+dy,1),mf_box(m+dz+dy,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));

				amrex::Vector<int> points;
				points.push_back(3);
#if AMREX_SPACEDIM > 2
				if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
				{
					/* Solving for end points */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
				StencilFill(stencil, traction, points, m+dy, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[3](n);

				for (int p = 1; p < points.size(); p++)
				{
					switch (points[p])
					{
#if AMREX_SPACEDIM > 2
						case 5: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dz+dy,n) = stencil[5](n);
							break;
						case 6: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dz+dy,n) = stencil[6](n);
							break;
#endif
					}
				}
			}
		}

		/* Step 4: Fill top face of ghost cells */
		if(BCUtil::IsNeumann(bc_hi[1]))
		{
			// non-end ghost cells. End ghost cells will be treated separately.
			int j = box.hiVect()[1] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0]+1; i <= box.hiVect()[0]-1; i++),
					,
					for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy,0),mf_box(m-dy,1),mf_box(m-dy,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx-dy,0),mf_box(m-dx-dy,1),mf_box(m-dx-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx-dy,0),mf_box(m+dx-dy,1),mf_box(m+dx-dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy-dy,0),mf_box(m-dy-dy,1),mf_box(m-dy-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy-dy,0),mf_box(m+dy-dy,1),mf_box(m+dy-dy,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz-dy,0),mf_box(m-dz-dy,1),mf_box(m-dz-dy,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz-dy,0),mf_box(m+dz-dy,1),mf_box(m+dz-dy,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));

				amrex::Vector<int> points;
				points.push_back(4);
#if AMREX_SPACEDIM > 2
				if (k == box.loVect()[2] && BCUtil::IsNeumann(bc_lo[2]))
				{
					/* Solving for end points */
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == box.hiVect()[2] && BCUtil::IsNeumann(bc_hi[2]))
				{
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
				StencilFill(stencil, traction, points, m-dy, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[4](n);

				for (int p = 1; p < points.size(); p++)
				{
					switch (points[p])
					{
#if AMREX_SPACEDIM > 2
						case 5: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m-dz-dy,n) = stencil[5](n);
							break;
						case 6: for (int n = 0; n < AMREX_SPACEDIM; n++)
								mf_box(m+dz-dy,n) = stencil[6](n);
							break;
#endif
					}
				}
			}
		}

#if AMREX_SPACEDIM > 2
		/* Step 5: Fill back face of ghost cells */
		if(BCUtil::IsNeumann(bc_lo[2]))
		{
			// End and triple end cells have already been filled - so no need to iterate over those
			int k = box.loVect()[2] - 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz,0),mf_box(m+dz,1),mf_box(m+dz,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx+dz,0),mf_box(m-dx+dz,1),mf_box(m-dx+dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx+dz,0),mf_box(m+dx+dz,1),mf_box(m+dx+dz,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy+dz,0),mf_box(m-dy+dz,1),mf_box(m-dy+dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy+dz,0),mf_box(m+dy+dz,1),mf_box(m+dy+dz,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz+dz,0),mf_box(m-dz+dz,1),mf_box(m-dz+dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz+dz,0),mf_box(m+dz+dz,1),mf_box(m+dz+dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));

				amrex::Vector<int> points;
				points.push_back(5);
				StencilFill(stencil, traction, points, m+dz, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[5](n);
			}
		}

		/* Step 6: Fill front face of ghost cells */
		if(BCUtil::IsNeumann(bc_hi[2]))
		{
			int k = box.hiVect()[2] + 1;
			AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
					for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				amrex::Vector<Set::Vector> stencil;
				stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz,0),mf_box(m-dz,1),mf_box(m-dz,2))));
				AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dx-dz,0),mf_box(m-dx-dz,1),mf_box(m-dx-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dx-dz,0),mf_box(m+dx-dz,1),mf_box(m+dx-dz,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dy-dz,0),mf_box(m-dy-dz,1),mf_box(m-dy-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dy-dz,0),mf_box(m+dy-dz,1),mf_box(m+dy-dz,2))));
						,
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m-dz-dz,0),mf_box(m-dz-dz,1),mf_box(m-dz-dz,2))));
						stencil.push_back(Set::Vector(AMREX_D_DECL(mf_box(m+dz-dz,0),mf_box(m+dz-dz,1),mf_box(m+dz-dz,2))));
					);
				amrex::Vector<Set::Vector> traction;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));

				amrex::Vector<int> points;
				points.push_back(6);
				StencilFill(stencil, traction, points, m-dz, m_amrlev, m_mglev, mfi);

				for (int n = 0; n<AMREX_SPACEDIM; n++)
					mf_box(m,n) = stencil[6](n);
			}
		}
#endif
#endif
	}
}

amrex::BCRec
Elastic::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
Elastic::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(bc_lo[0]),BCUtil::IsPeriodic(bc_lo[1]),BCUtil::IsPeriodic(bc_lo[2]))};
}


// Stencil Fill routine - takes in a stencil, a list of unknown points and fills the unknown values
// in the stencil.
#define C(i,j,k,l,m,amrlev,mglev,mfi) Operator::Elastic::Elastic::C(i,j,k,l,m,amrlev,mglev,mfi)
void 
Elastic::StencilFill(	amrex::Vector<Set::Vector> &stencil,
			const amrex::Vector<Set::Vector> &traction,
			const amrex::Vector<int> &points,
			const amrex::IntVect &m,
			const int amrlev,
			const int mglev,
			const amrex::MFIter &mfi)
{

/*
Description of arguments

stencil:	list of 7 points, each with three components
traction:	list of tractions, each with three components
points:		list of points that are empty (can't be more than 3)
m:		Position of the middle point
amrlev:		AMR level
mglev:		MG level
mfi:		MFI box

point nomenclature:
	0 = mid
	1 = left
	2 = right
	3 = bottom
	4 = top
	5 = back
	6 = front

Restriction: If there are more than one unknown points, 
		the list must be in asceding order. 
*/

	if(points.size() > 3 || points.size() < 1)
		Util::Abort("Number of unknown points can not be greater than 3");

	if(points.size() != traction.size())
		Util::Abort("Mismatch between number of unknown points and tractions");

	const Real* DX = m_geom[amrlev][mglev].CellSize();
  
#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	if(points.size() == 1)		//Need three equations - non corner cases
	{
		int mul = 0;
		if(points[0] == 1 || points[0] == 2)		// left or right faces
		{
			mul = points[0] == 1 ? -1 : 1;
			//gradu_k,2
			Set::Vector gradu_2;
			AMREX_D_TERM(	gradu_2(0) = (stencil[4](0) - stencil[3](0))/(2.0*DX[1]);
					, // 2D
					gradu_2(1) = (stencil[4](1) - stencil[3](1))/(2.0*DX[1]);
					, // 3D
					gradu_2(2) = (stencil[4](2) - stencil[3](2))/(2.0*DX[1]);
					//gradu_k,3
					Set::Vector gradu_3;
					gradu_3(0) = (stencil[6](0) - stencil[5](0))/(2.0*DX[2]);
					gradu_3(1) = (stencil[6](1) - stencil[5](1))/(2.0*DX[2]);
					gradu_3(2) = (stencil[6](2) - stencil[5](2))/(2.0*DX[2]););

			Set::Matrix left;
			Set::Vector right; 
			Set::Vector sol;

			AMREX_D_TERM( 	left(0,0) = C(0,0,0,0,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = C(0,0,1,0,m,amrlev,mglev,mfi);
					left(1,0) = C(1,0,0,0,m,amrlev,mglev,mfi);
					left(1,1) = C(1,0,1,0,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = C(0,0,2,0,m,amrlev,mglev,mfi);
					left(1,2) = C(1,0,2,0,m,amrlev,mglev,mfi);
					left(2,0) = C(2,0,0,0,m,amrlev,mglev,mfi);
					left(2,1) = C(2,0,1,0,m,amrlev,mglev,mfi);
					left(2,2) = C(2,0,2,0,m,amrlev,mglev,mfi););
			AMREX_D_TERM(	right(0) = AMREX_D_DECL(mul*traction[0](0)
								, 
								- C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								);
					, // 2D
					right(1) = AMREX_D_DECL(mul*traction[0](1)
								, 
								- C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								);
					, // 3D
					right(2) = AMREX_D_DECL(mul*traction[0](2)
								, 
								- C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
					););
			sol = left.ldlt().solve(right); // we can change this solver as needed
					
			AMREX_D_TERM(	stencil[points[0]](0) = stencil[0](0) + mul*DX[0]*sol(0);,
					stencil[points[0]](1) = stencil[0](1) + mul*DX[0]*sol(1);,
					stencil[points[0]](2) = stencil[0](2) + mul*DX[0]*sol(2););
		}

		else if(points[0] == 3 || points[0] == 4)	// bottom or top faces
		{
			mul = points[0] == 3 ? -1 : 1;
			//gradu_k,1
			Set::Vector gradu_1;
			AMREX_D_TERM(	gradu_1(0) = (stencil[2](0) - stencil[1](0))/(2.0*DX[0]);
					, // 2D
					gradu_1(1) = (stencil[2](1) - stencil[1](1))/(2.0*DX[0]);
					, // 3D
					gradu_1(2) = (stencil[2](2) - stencil[1](2))/(2.0*DX[0]);
					//gradu_k,3
					Set::Vector gradu_3;
					gradu_3(0) = (stencil[6](0) - stencil[5](0))/(2.0*DX[2]);
					gradu_3(1) = (stencil[6](1) - stencil[5](1))/(2.0*DX[2]);
					gradu_3(2) = (stencil[6](2) - stencil[5](2))/(2.0*DX[2]););

			Set::Matrix left;
			Set::Vector right; 
			Set::Vector sol;
			AMREX_D_TERM( 	left(0,0) = C(0,1,0,1,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = C(0,1,1,1,m,amrlev,mglev,mfi);
					left(1,0) = C(1,1,0,1,m,amrlev,mglev,mfi);
					left(1,1) = C(1,1,1,1,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = C(0,1,2,1,m,amrlev,mglev,mfi);
					left(1,2) = C(1,1,2,1,m,amrlev,mglev,mfi);
					left(2,0) = C(2,1,0,1,m,amrlev,mglev,mfi);
					left(2,1) = C(2,1,1,1,m,amrlev,mglev,mfi);
					left(2,2) = C(2,1,2,1,m,amrlev,mglev,mfi););
			AMREX_D_TERM(	right(0) = AMREX_D_DECL(mul*traction[0](0)
								, 
								- C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 2D
					right(1) = AMREX_D_DECL(mul*traction[0](1)
								, 
								- C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 3D
					right(2) = AMREX_D_DECL(mul*traction[0](2)
								, 
								- C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
					););
			sol = left.ldlt().solve(right); // we can change this solver as needed
			
			AMREX_D_TERM(	stencil[points[0]](0) = stencil[0](0) + mul*DX[1]*sol(0);,
					stencil[points[0]](1) = stencil[0](1) + mul*DX[1]*sol(1);,
					stencil[points[0]](2) = stencil[0](2) + mul*DX[1]*sol(2););
		}

		else if(points[0] == 5 || points[0] == 6)	// back and front faces
		{
			mul = points[0] == 5 ? -1 : 1;
			//gradu_k,1
			Set::Vector gradu_1;
			AMREX_D_TERM(	gradu_1(0) = (stencil[2](0) - stencil[1](0))/(2.0*DX[0]);
					, // 2D
					gradu_1(1) = (stencil[2](1) - stencil[1](1))/(2.0*DX[0]);
					, // 3D
					gradu_1(2) = (stencil[2](2) - stencil[1](2))/(2.0*DX[0]);
					//gradu_k,2
					Set::Vector gradu_2;
					gradu_2(0) = (stencil[4](0) - stencil[3](0))/(2.0*DX[1]);
					gradu_2(1) = (stencil[4](1) - stencil[3](1))/(2.0*DX[1]);
					gradu_2(2) = (stencil[4](2) - stencil[3](2))/(2.0*DX[1]););

			Set::Matrix left;
			Set::Vector right; 
			Set::Vector sol;
			AMREX_D_TERM( 	left(0,0) = C(0,2,0,2,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = C(0,2,1,2,m,amrlev,mglev,mfi);
					left(1,0) = C(1,2,0,2,m,amrlev,mglev,mfi);
					left(1,1) = C(1,2,1,2,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = C(0,2,2,2,m,amrlev,mglev,mfi);
					left(1,2) = C(1,2,2,2,m,amrlev,mglev,mfi);
					left(2,0) = C(2,2,0,2,m,amrlev,mglev,mfi);
					left(2,1) = C(2,2,1,2,m,amrlev,mglev,mfi);
							left(2,2) = C(2,2,2,2,m,amrlev,mglev,mfi););
					AMREX_D_TERM(	right(0) = AMREX_D_DECL(mul*traction[0](0)
										, 
										- C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
										- C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
										,
										- C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
										- C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
										- C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
										- C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
										);
							, // 2D
							right(1) = AMREX_D_DECL(mul*traction[0](1)
										, 
										- C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
										- C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
										,
										- C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
										- C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
										- C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
										- C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
										);
							, // 3D
							right(2) = AMREX_D_DECL(mul*traction[0](2)
										, 
										- C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
										- C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
										,
										- C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
										- C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
										- C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
										- C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
										););
					sol = left.ldlt().solve(right); // we can change this solver as needed
					
					AMREX_D_TERM(	stencil[points[0]](0) = stencil[0](0) + mul*DX[2]*sol(0);,
							stencil[points[0]](1) = stencil[0](1) + mul*DX[2]*sol(1);,
							stencil[points[0]](2) = stencil[0](2) + mul*DX[2]*sol(2););
		}

		else
			Util::Abort("Incorrect values of points");
	}

	else if (points.size() == 2)	//Need six equations - corner case
	{
		int mul1 = 0, mul2 = 0;
		if ((points[0] == 1 || points[0] == 2) && (points[1] == 3 || points[1] == 4))
		{
			mul1 = points[0] == 1 ? -1 : 1;
			mul2 = points[1] == 3 ? -1 : 1;
#if AMREX_SPACEDIM > 2
			//gradu_k,3
			Set::Vector gradu_3;
			gradu_3(0) = (stencil[6](0) - stencil[5](0))/(2.0*DX[2]);
			gradu_3(1) = (stencil[6](1) - stencil[5](1))/(2.0*DX[2]);
			gradu_3(2) = (stencil[6](2) - stencil[5](2))/(2.0*DX[2]);
#endif
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,2*AMREX_SPACEDIM> left;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> right;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> sol; 

			AMREX_D_TERM(	left(0,0) = C(0,0,0,0,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = C(0,0,1,0,m,amrlev,mglev,mfi);
					left(0,2) = C(0,0,0,1,m,amrlev,mglev,mfi);
					left(0,3) = C(0,0,1,1,m,amrlev,mglev,mfi);
					left(1,0) = C(1,0,0,0,m,amrlev,mglev,mfi);
					left(1,1) = C(1,0,1,0,m,amrlev,mglev,mfi);
					left(1,2) = C(1,0,0,1,m,amrlev,mglev,mfi);
					left(1,3) = C(1,0,1,1,m,amrlev,mglev,mfi);
					left(2,0) = C(0,1,0,0,m,amrlev,mglev,mfi);
					left(2,1) = C(0,1,1,0,m,amrlev,mglev,mfi);
					left(2,2) = C(0,1,0,1,m,amrlev,mglev,mfi);
					left(2,3) = C(0,1,1,1,m,amrlev,mglev,mfi);
					left(3,0) = C(1,1,0,0,m,amrlev,mglev,mfi);
					left(3,1) = C(1,1,1,0,m,amrlev,mglev,mfi);
					left(3,2) = C(1,1,0,1,m,amrlev,mglev,mfi);
					left(3,3) = C(1,1,1,1,m,amrlev,mglev,mfi);
					, // 3D
					left(0,4) = C(0,0,2,0,m,amrlev,mglev,mfi);
					left(0,5) = C(0,0,2,1,m,amrlev,mglev,mfi);
					left(1,4) = C(1,0,2,0,m,amrlev,mglev,mfi);
					left(1,5) = C(1,0,2,1,m,amrlev,mglev,mfi);
					left(2,4) = C(0,1,2,0,m,amrlev,mglev,mfi);
					left(2,5) = C(0,1,2,1,m,amrlev,mglev,mfi);
					left(3,4) = C(1,1,2,0,m,amrlev,mglev,mfi);
					left(3,5) = C(1,1,2,1,m,amrlev,mglev,mfi);
					left(4,0) = C(2,0,0,0,m,amrlev,mglev,mfi);
					left(4,1) = C(2,0,1,0,m,amrlev,mglev,mfi);
					left(4,2) = C(2,0,0,1,m,amrlev,mglev,mfi);
					left(4,3) = C(2,0,1,1,m,amrlev,mglev,mfi);
					left(4,4) = C(2,0,2,0,m,amrlev,mglev,mfi);
					left(4,5) = C(2,0,2,1,m,amrlev,mglev,mfi);
					left(5,0) = C(2,1,0,0,m,amrlev,mglev,mfi);
					left(5,1) = C(2,1,1,0,m,amrlev,mglev,mfi);
					left(5,2) = C(2,1,0,1,m,amrlev,mglev,mfi);
					left(5,3) = C(2,1,1,1,m,amrlev,mglev,mfi);
					left(5,4) = C(2,1,2,0,m,amrlev,mglev,mfi);
					left(5,5) = C(2,1,2,1,m,amrlev,mglev,mfi););

			AMREX_D_TERM(	right(0) = AMREX_D_DECL(	mul1*traction[0](0)
									, // 2D
									+ 0.0
									, // 3D
									- C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								);
					, // 2D
					right(1) = AMREX_D_DECL(	mul1*traction[0](1)
									, // 2D
									+ 0.0
									, // 3D
									- C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								);
					right(2) = AMREX_D_DECL(	mul2*traction[1](0)
									, // 2D
									+ 0.0
									, // 3D
									- C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								);
					right(3) = AMREX_D_DECL(	mul2*traction[1](1)
									, // 2D
									+ 0.0
									, // 3D
									- C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								);
					, // 3D
					right(4) = AMREX_D_DECL(	mul1*traction[0](2)
									, // 2D
									+ 0.0
									, // 3D
									- C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								);
					right(5) = AMREX_D_DECL(	mul2*traction[1](2)
									, // 2D
									+ 0.0
									, // 3D
									- C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								););

			sol = left.ldlt().solve(right); // we can change this solver as needed
			AMREX_D_TERM(	stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
					,
					stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
					stencil[points[1]](0) = stencil[0](0) + mul2*DX[1]*sol(2);
					stencil[points[1]](1) = stencil[0](1) + mul2*DX[1]*sol(3);
					,
					stencil[points[0]](2) = stencil[0](2) + mul1*DX[2]*sol(4);
					stencil[points[1]](2) = stencil[0](2) + mul2*DX[2]*sol(5););
		} 

#if AMREX_SPACEDIM > 2
		else if ((points[0] == 1 || points[0] == 2) && (points[1] == 5 || points[1] == 6)) 
		{
			mul1 = points[0] == 1 ? -1 : 1;
			mul2 = points[1] == 5 ? -1 : 1;
			//gradu_k,2
			Set::Vector gradu_2;
			gradu_2(0) = (stencil[4](0) - stencil[3](0))/(2.0*DX[1]);
			gradu_2(1) = (stencil[4](1) - stencil[3](1))/(2.0*DX[1]);
			gradu_2(2) = (stencil[4](2) - stencil[3](2))/(2.0*DX[1]);

			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,2*AMREX_SPACEDIM> left;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> right;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> sol; 

			left(0,0) = C(0,0,0,0,m,amrlev,mglev,mfi);
			left(0,1) = C(0,0,1,0,m,amrlev,mglev,mfi);
			left(0,2) = C(0,0,2,0,m,amrlev,mglev,mfi);
			left(0,3) = C(0,0,0,2,m,amrlev,mglev,mfi);
			left(0,4) = C(0,0,1,2,m,amrlev,mglev,mfi);
			left(0,5) = C(0,0,2,2,m,amrlev,mglev,mfi);
			left(1,0) = C(1,0,0,0,m,amrlev,mglev,mfi);
			left(1,1) = C(1,0,1,0,m,amrlev,mglev,mfi);
			left(1,2) = C(1,0,2,0,m,amrlev,mglev,mfi);
			left(1,3) = C(1,0,0,2,m,amrlev,mglev,mfi);
			left(1,4) = C(1,0,1,2,m,amrlev,mglev,mfi);
			left(1,5) = C(1,0,2,2,m,amrlev,mglev,mfi);
			left(2,0) = C(2,0,0,0,m,amrlev,mglev,mfi);
			left(2,1) = C(2,0,1,0,m,amrlev,mglev,mfi);
			left(2,2) = C(2,0,2,0,m,amrlev,mglev,mfi);
			left(2,3) = C(2,0,0,2,m,amrlev,mglev,mfi);
			left(2,4) = C(2,0,1,2,m,amrlev,mglev,mfi);
			left(2,5) = C(2,0,2,2,m,amrlev,mglev,mfi);
			left(3,0) = C(0,2,0,0,m,amrlev,mglev,mfi);
			left(3,1) = C(0,2,1,0,m,amrlev,mglev,mfi);
			left(3,2) = C(0,2,2,0,m,amrlev,mglev,mfi);
			left(3,3) = C(0,2,0,2,m,amrlev,mglev,mfi);
			left(3,4) = C(0,2,1,2,m,amrlev,mglev,mfi);
			left(3,5) = C(0,2,2,2,m,amrlev,mglev,mfi);
			left(4,0) = C(1,2,0,0,m,amrlev,mglev,mfi);
			left(4,1) = C(1,2,1,0,m,amrlev,mglev,mfi);
			left(4,2) = C(1,2,2,0,m,amrlev,mglev,mfi);
			left(4,3) = C(1,2,0,2,m,amrlev,mglev,mfi);
			left(4,4) = C(1,2,1,2,m,amrlev,mglev,mfi);
			left(4,5) = C(1,2,2,2,m,amrlev,mglev,mfi);
			left(5,0) = C(2,2,0,0,m,amrlev,mglev,mfi);
			left(5,1) = C(2,2,1,0,m,amrlev,mglev,mfi);
			left(5,2) = C(2,2,2,0,m,amrlev,mglev,mfi);
			left(5,3) = C(2,2,0,2,m,amrlev,mglev,mfi);
			left(5,4) = C(2,2,1,2,m,amrlev,mglev,mfi);
			left(5,5) = C(2,2,2,2,m,amrlev,mglev,mfi);

			right(0) = 	mul1*traction[0](0)
					- C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(1) = 	mul1*traction[0](1)
					- C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(2) = 	mul1*traction[0](2)
					- C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(3) = 	mul2*traction[1](0)
					- C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(4) = 	mul2*traction[1](1)
					- C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(5) = 	mul2*traction[1](2)
					- C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
					- C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
					- C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);

			sol = left.ldlt().solve(right); // we can change this solver as needed
			
			stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
			stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
			stencil[points[0]](2) = stencil[0](2) + mul1*DX[0]*sol(2);

			stencil[points[1]](0) = stencil[0](0) + mul2*DX[2]*sol(0);
			stencil[points[1]](1) = stencil[0](1) + mul2*DX[2]*sol(1);
			stencil[points[1]](2) = stencil[0](2) + mul2*DX[2]*sol(2);
		}

		else if ((points[0] == 3 || points[0] == 4) && (points[1] == 5 || points[1] == 6)) 
		{
			mul1 = points[0] == 3 ? -1 : 1;
			mul2 = points[1] == 5 ? -1 : 1;

			//gradu_k,1
			Set::Vector gradu_1;
			gradu_1(0) = (stencil[2](0) - stencil[1](0))/(2.0*DX[1]);
			gradu_1(1) = (stencil[2](1) - stencil[1](1))/(2.0*DX[1]);
			gradu_1(2) = (stencil[2](2) - stencil[1](2))/(2.0*DX[1]);

			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,2*AMREX_SPACEDIM> left;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> right;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> sol; 

			left(0,0) = C(0,1,0,1,m,amrlev,mglev,mfi);
			left(0,1) = C(0,1,1,1,m,amrlev,mglev,mfi);
			left(0,2) = C(0,1,2,1,m,amrlev,mglev,mfi);
			left(0,3) = C(0,1,0,2,m,amrlev,mglev,mfi);
			left(0,4) = C(0,1,1,2,m,amrlev,mglev,mfi);
			left(0,5) = C(0,1,2,2,m,amrlev,mglev,mfi);
			left(1,0) = C(1,1,0,1,m,amrlev,mglev,mfi);
			left(1,1) = C(1,1,1,1,m,amrlev,mglev,mfi);
			left(1,2) = C(1,1,2,1,m,amrlev,mglev,mfi);
			left(1,3) = C(1,1,0,2,m,amrlev,mglev,mfi);
			left(1,4) = C(1,1,1,2,m,amrlev,mglev,mfi);
			left(1,5) = C(1,1,2,2,m,amrlev,mglev,mfi);
			left(2,0) = C(2,1,0,1,m,amrlev,mglev,mfi);
			left(2,1) = C(2,1,1,1,m,amrlev,mglev,mfi);
			left(2,2) = C(2,1,2,1,m,amrlev,mglev,mfi);
			left(2,3) = C(2,1,0,2,m,amrlev,mglev,mfi);
			left(2,4) = C(2,1,1,2,m,amrlev,mglev,mfi);
			left(2,5) = C(2,1,2,2,m,amrlev,mglev,mfi);
			left(3,0) = C(0,2,0,1,m,amrlev,mglev,mfi);
			left(3,1) = C(0,2,1,1,m,amrlev,mglev,mfi);
			left(3,2) = C(0,2,2,1,m,amrlev,mglev,mfi);
			left(3,3) = C(0,2,0,2,m,amrlev,mglev,mfi);
			left(3,4) = C(0,2,1,2,m,amrlev,mglev,mfi);
			left(3,5) = C(0,2,2,2,m,amrlev,mglev,mfi);
			left(4,0) = C(1,2,0,1,m,amrlev,mglev,mfi);
			left(4,1) = C(1,2,1,1,m,amrlev,mglev,mfi);
			left(4,2) = C(1,2,2,1,m,amrlev,mglev,mfi);
			left(4,3) = C(1,2,0,2,m,amrlev,mglev,mfi);
			left(4,4) = C(1,2,1,2,m,amrlev,mglev,mfi);
			left(4,5) = C(1,2,2,2,m,amrlev,mglev,mfi);
			left(5,0) = C(2,2,0,1,m,amrlev,mglev,mfi);
			left(5,1) = C(2,2,1,1,m,amrlev,mglev,mfi);
			left(5,2) = C(2,2,2,1,m,amrlev,mglev,mfi);
			left(5,3) = C(2,2,0,2,m,amrlev,mglev,mfi);
			left(5,4) = C(2,2,1,2,m,amrlev,mglev,mfi);
			left(5,5) = C(2,2,2,2,m,amrlev,mglev,mfi);

			right(0) = 	mul1*traction[0](0)
					- C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(1) = 	mul1*traction[0](1)
					- C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(2) = 	mul1*traction[0](2)
					- C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(3) = 	mul2*traction[1](0)
					- C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(4) = 	mul2*traction[1](1)
					- C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(5) = 	mul2*traction[1](2)
					- C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
					- C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
					- C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);

			sol = left.ldlt().solve(right); // we can change this solver as needed
			
			stencil[points[0]](0) = stencil[0](0) + mul1*DX[1]*sol(0);
			stencil[points[0]](1) = stencil[0](1) + mul1*DX[1]*sol(1);
			stencil[points[0]](2) = stencil[0](2) + mul1*DX[1]*sol(2);

			stencil[points[1]](0) = stencil[0](0) + mul2*DX[2]*sol(0);
			stencil[points[1]](1) = stencil[0](1) + mul2*DX[2]*sol(1);
			stencil[points[1]](2) = stencil[0](2) + mul2*DX[2]*sol(2);
		}
#endif
		else
			Util::Abort("Corner case: Incorrect list/values of points");
	}

#if AMREX_SPACEDIM > 2
	else if (points.size() == 3)	//Need 9 equations - triple point case
	{
		int mul1 = points[0] == 1 ? -1 : 1;
		int mul2 = points[1] == 3 ? -1 : 1;
		int mul3 = points[2] == 5 ? -1 : 1;

		Eigen::Matrix<amrex::Real,3*AMREX_SPACEDIM,3*AMREX_SPACEDIM> left;
		Eigen::Matrix<amrex::Real,3*AMREX_SPACEDIM,1> right;
		Eigen::Matrix<amrex::Real,3*AMREX_SPACEDIM,1> sol; 

		left(0,0) = C(0,0,0,0,m,amrlev,mglev,mfi);
		left(0,1) = C(0,0,1,0,m,amrlev,mglev,mfi);
		left(0,2) = C(0,0,2,0,m,amrlev,mglev,mfi);
		left(0,3) = C(0,0,0,1,m,amrlev,mglev,mfi);
		left(0,4) = C(0,0,1,1,m,amrlev,mglev,mfi);
		left(0,5) = C(0,0,2,1,m,amrlev,mglev,mfi);
		left(0,6) = C(0,0,0,2,m,amrlev,mglev,mfi);
		left(0,7) = C(0,0,2,2,m,amrlev,mglev,mfi);
		left(0,8) = C(0,0,2,2,m,amrlev,mglev,mfi);
		left(1,0) = C(1,0,0,0,m,amrlev,mglev,mfi);
		left(1,1) = C(1,0,1,0,m,amrlev,mglev,mfi);
		left(1,2) = C(1,0,2,0,m,amrlev,mglev,mfi);
		left(1,3) = C(1,0,0,1,m,amrlev,mglev,mfi);
		left(1,4) = C(1,0,1,1,m,amrlev,mglev,mfi);
		left(1,5) = C(1,0,2,1,m,amrlev,mglev,mfi);
		left(1,6) = C(1,0,0,2,m,amrlev,mglev,mfi);
		left(1,7) = C(1,0,2,2,m,amrlev,mglev,mfi);
		left(1,8) = C(1,0,2,2,m,amrlev,mglev,mfi);
		left(2,0) = C(2,0,0,0,m,amrlev,mglev,mfi);
		left(2,1) = C(2,0,1,0,m,amrlev,mglev,mfi);
		left(2,2) = C(2,0,2,0,m,amrlev,mglev,mfi);
		left(2,3) = C(2,0,0,1,m,amrlev,mglev,mfi);
		left(2,4) = C(2,0,1,1,m,amrlev,mglev,mfi);
		left(2,5) = C(2,0,2,1,m,amrlev,mglev,mfi);
		left(2,6) = C(2,0,0,2,m,amrlev,mglev,mfi);
		left(2,7) = C(2,0,1,2,m,amrlev,mglev,mfi);
		left(2,8) = C(2,0,2,2,m,amrlev,mglev,mfi);
		left(3,0) = C(0,1,0,0,m,amrlev,mglev,mfi);
		left(3,1) = C(0,1,1,0,m,amrlev,mglev,mfi);
		left(3,2) = C(0,1,2,0,m,amrlev,mglev,mfi);
		left(3,3) = C(0,1,0,1,m,amrlev,mglev,mfi);
		left(3,4) = C(0,1,1,1,m,amrlev,mglev,mfi);
		left(3,5) = C(0,1,2,1,m,amrlev,mglev,mfi);
		left(3,6) = C(0,1,0,2,m,amrlev,mglev,mfi);
		left(3,7) = C(0,1,1,2,m,amrlev,mglev,mfi);
		left(3,8) = C(0,1,2,2,m,amrlev,mglev,mfi);
		left(4,0) = C(1,1,0,0,m,amrlev,mglev,mfi);
		left(4,1) = C(1,1,1,0,m,amrlev,mglev,mfi);
		left(4,2) = C(1,1,2,0,m,amrlev,mglev,mfi);
		left(4,3) = C(1,1,0,1,m,amrlev,mglev,mfi);
		left(4,4) = C(1,1,1,1,m,amrlev,mglev,mfi);
		left(4,5) = C(1,1,2,1,m,amrlev,mglev,mfi);
		left(4,6) = C(1,1,0,2,m,amrlev,mglev,mfi);
		left(4,7) = C(1,1,1,2,m,amrlev,mglev,mfi);
		left(4,8) = C(1,1,2,2,m,amrlev,mglev,mfi);
		left(5,0) = C(2,1,0,0,m,amrlev,mglev,mfi);
		left(5,1) = C(2,1,1,0,m,amrlev,mglev,mfi);
		left(5,2) = C(2,1,2,0,m,amrlev,mglev,mfi);
		left(5,3) = C(2,1,0,1,m,amrlev,mglev,mfi);
		left(5,4) = C(2,1,1,1,m,amrlev,mglev,mfi);
		left(5,5) = C(2,1,2,1,m,amrlev,mglev,mfi);
		left(5,6) = C(2,1,0,2,m,amrlev,mglev,mfi);
		left(5,7) = C(2,1,1,2,m,amrlev,mglev,mfi);
		left(5,8) = C(2,1,2,2,m,amrlev,mglev,mfi);
		left(6,0) = C(0,2,0,0,m,amrlev,mglev,mfi);
		left(6,1) = C(0,2,1,0,m,amrlev,mglev,mfi);
		left(6,2) = C(0,2,2,0,m,amrlev,mglev,mfi);
		left(6,3) = C(0,2,0,1,m,amrlev,mglev,mfi);
		left(6,4) = C(0,2,1,1,m,amrlev,mglev,mfi);
		left(6,5) = C(0,2,2,1,m,amrlev,mglev,mfi);
		left(6,6) = C(0,2,0,2,m,amrlev,mglev,mfi);
		left(6,7) = C(0,2,1,2,m,amrlev,mglev,mfi);
		left(6,8) = C(0,2,2,2,m,amrlev,mglev,mfi);
		left(7,0) = C(1,2,0,0,m,amrlev,mglev,mfi);
		left(7,1) = C(1,2,1,0,m,amrlev,mglev,mfi);
		left(7,2) = C(1,2,2,0,m,amrlev,mglev,mfi);
		left(7,3) = C(1,2,0,1,m,amrlev,mglev,mfi);
		left(7,4) = C(1,2,1,1,m,amrlev,mglev,mfi);
		left(7,5) = C(1,2,2,1,m,amrlev,mglev,mfi);
		left(7,6) = C(1,2,0,2,m,amrlev,mglev,mfi);
		left(7,7) = C(1,2,1,2,m,amrlev,mglev,mfi);
		left(7,8) = C(1,2,2,2,m,amrlev,mglev,mfi);
		left(7,0) = C(2,2,0,0,m,amrlev,mglev,mfi);
		left(8,1) = C(2,2,1,0,m,amrlev,mglev,mfi);
		left(8,2) = C(2,2,2,0,m,amrlev,mglev,mfi);
		left(8,3) = C(2,2,0,1,m,amrlev,mglev,mfi);
		left(8,4) = C(2,2,1,1,m,amrlev,mglev,mfi);
		left(8,5) = C(2,2,2,1,m,amrlev,mglev,mfi);
		left(8,6) = C(2,2,0,2,m,amrlev,mglev,mfi);
		left(8,7) = C(2,2,1,2,m,amrlev,mglev,mfi);
		left(8,8) = C(2,2,2,2,m,amrlev,mglev,mfi);

		right(0) = mul1*traction[0](0);
		right(1) = mul1*traction[0](1);
		right(2) = mul1*traction[0](2);
		right(3) = mul2*traction[1](0);
		right(4) = mul2*traction[1](1);
		right(5) = mul2*traction[1](2);
		right(6) = mul3*traction[2](0);
		right(7) = mul3*traction[2](1);
		right(8) = mul3*traction[2](2);

		sol = left.ldlt().solve(right); // we can change this solver as needed

		stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
		stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
		stencil[points[0]](2) = stencil[0](2) + mul1*DX[0]*sol(2);

		stencil[points[1]](0) = stencil[0](0) + mul2*DX[1]*sol(0);
		stencil[points[1]](1) = stencil[0](1) + mul2*DX[1]*sol(1);
		stencil[points[1]](2) = stencil[0](2) + mul2*DX[1]*sol(2);

		stencil[points[2]](0) = stencil[0](0) + mul3*DX[2]*sol(0);
		stencil[points[2]](1) = stencil[0](1) + mul3*DX[2]*sol(1);
		stencil[points[2]](2) = stencil[0](2) + mul3*DX[2]*sol(2);
	}
#endif
}
}
