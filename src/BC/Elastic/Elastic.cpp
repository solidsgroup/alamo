#include "Elastic.H"

namespace BC
{

Elastic::Elastic(amrex::Vector<std::string> _bc_hi_str,
		 amrex::Vector<std::string> _bc_lo_str,
		 AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
			      amrex::Vector<amrex::Real> _bc_lo_2,
			      amrex::Vector<amrex::Real> _bc_lo_3),
		 AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
			      amrex::Vector<amrex::Real> _bc_hi_2,
			      amrex::Vector<amrex::Real> _bc_hi_3))
	:
	bc_hi_str(_bc_hi_str),
	bc_lo_str(_bc_lo_str), 
	AMREX_D_DECL(bc_lo_1(_bc_lo_1),bc_lo_2(_bc_lo_2),bc_lo_3(_bc_lo_3)),
	AMREX_D_DECL(bc_hi_1(_bc_hi_1),bc_hi_2(_bc_hi_2),bc_hi_3(_bc_hi_3))
{
	for (int i=0;i<AMREX_SPACEDIM;i++)
	{
		bc_hi[i] = BCUtil::ReadString(bc_hi_str[i]);
		bc_lo[i] = BCUtil::ReadString(bc_lo_str[i]);

		if (BCUtil::IsPeriodic(bc_lo[i]) != BCUtil::IsPeriodic(bc_hi[i]))
			Util::Abort("Invalid BCs cannot be periodic on one side and not the other");
	}
}

void
Elastic::SetElasticOperator(Operator::Elastic::Elastic* a_operator)
{
	m_operator = a_operator;
}

//void
//Elastic::FillBoundary (amrex::MultiFab& mf, amrex::Real time)
void Elastic::FillBoundary (amrex::FArrayBox &mf_box,
			    const amrex::Box &box,
			    int ngrow, int dcomp, int ncomp,
			    amrex::Real time,
			    amrex::MFIter &mfi,
			    Orientation face,
			    const amrex::Mask *mask)
{
	/* 
	   We want to fill all dirichlet BC first and then Neumann.
	   This is to ensure that when the stencil function is called, 
	   cells that can be filled, have been filled.
	   It is further assumed that ncomp = 3.
	*/
	
	amrex::Box domain(m_geom.Domain());

	std::cout << "Box loVect = (" << box.loVect()[0] << "," << box.loVect()[1] << "," << box.loVect()[2] << ")" << std::endl;
	std::cout << "Box hiVect = (" << box.hiVect()[0] << "," << box.hiVect()[1] << "," << box.hiVect()[2] << ")" << std::endl;
	std::cout << "Domain loVect = (" << domain.loVect()[0] << "," << domain.loVect()[1] << "," << domain.loVect()[2] << ")" << std::endl;
	std::cout << "Domain hiVect = (" << domain.hiVect()[0] << "," << domain.hiVect()[1] << "," << domain.hiVect()[2] << ")" << std::endl;
	std::cout << "Amrlev = " << m_amrlev << ". Mglev = " << m_mglev << std::endl;

	//mf.FillBoundary(m_geom.periodicity());

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	// for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
	// {
	// 	const amrex::Box& box = mfi.tilebox();

	//amrex::BaseFab<amrex::Real> &mf_box = mf[mfi];

	/* Dirichlet boundaries are first */
	AMREX_D_TERM(	for(int i = box.loVect()[0] - ngrow; i<=box.hiVect()[0] + ngrow; i++),
					for (int j = box.loVect()[1] - ngrow; j<=box.hiVect()[1] + ngrow; j++),
					for (int k = box.loVect()[2] - ngrow; k<=box.hiVect()[2] + ngrow; k++))
	{
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		if (i == domain.loVect()[0]-1 && (face == Orientation::xlo || face == Orientation::All) && bc_lo_str[0] == "dirichlet") // Left boundary
		{
			//std::cout << "Dirichet boundary in left face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_lo_1[n];
		}
		if (i == domain.hiVect()[0]+1 && (face == Orientation::xhi || face == Orientation::All) && bc_hi_str[0] == "dirichlet") // Right boundary
		{
			//std::cout << "Dirichet boundary in right face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_hi_1[n];
		}
#if AMREX_SPACEDIM>1
		if (j == domain.loVect()[1]-1 && (face == Orientation::ylo || face == Orientation::All) && bc_lo_str[1] == "dirichlet") // Bottom boundary
		{
			//std::cout << "Dirichet boundary in bottom face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_lo_2[n];
		}
		if (j == domain.hiVect()[1]+1 && (face == Orientation::yhi || face == Orientation::All) && bc_hi_str[1] == "dirichlet") // Top boundary
		{
			//std::cout << "Dirichet boundary in top face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_hi_2[n];
		}
#if AMREX_SPACEDIM>2
		if (k == domain.loVect()[2]-1 && (face == Orientation::zlo || face == Orientation::All) && bc_lo_str[2] == "dirichlet") // Back boundary
		{
			//std::cout << "Dirichet boundary in back face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_lo_3[n];
		}
		if (k == domain.hiVect()[2]+1 && (face == Orientation::zhi || face == Orientation::All) && bc_hi_str[2] == "dirichlet") // Front boundary
		{
			//std::cout << "Dirichet boundary in front face" << std::endl;
			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = bc_hi_3[n];
		}
#endif
#endif
	}

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

	/* The following steps are for debugging purposes.
		They can be disabled by setting debug = false*/
	bool debug = true;

	amrex::IntVect point_left(-1,3,3);
	amrex::IntVect point_right(8,3,3);
	amrex::IntVect point_bottom(3,-1,3);
	amrex::IntVect point_top(3,8,3);
	amrex::IntVect point_back(3,3,-1);
	amrex::IntVect point_front(3,3,8);

	/* Step 1: Fill the left face ghost cells */
	
	if (bc_lo_str[0] == "traction" && box.loVect()[0]==domain.loVect()[0])
	{
		int i = box.loVect()[0] - 1;
		//std::cout << "Neumann boundary in left face. i = " << i << std::endl;

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
			if (j == domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				//std::cout << "End point: Neumann boundary in bottom face. i,j,k = "<< i << "," << j << "," << k << std::endl;
				/* Solving for end points */
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
#if AMREX_SPACEDIM > 2
				/* Solving for triple end point */ 
				if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
				{
					//std::cout << "Triple End point: Neumann boundary in back face. i,j,k = "<< i << "," << j << "," << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
				{
					//std::cout << "Triple End point: Neumann boundary in front face. i,j,k = "<< i << "," << j << "," << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
			}
			else if (j == domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				//std::cout << "End point: Neumann boundary in bottom face. i,j,k = "<< i << "," << j << "," << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
#if AMREX_SPACEDIM > 2
				if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
				{
					//std::cout << "Triple End point: Neumann boundary in back face. i,j,k = "<< i << "," << j << "," << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
				{
					//std::cout << "Triple End point: Neumann boundary in front face. i,j,k = "<< i << "," << j << "," << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
			}
#if AMREX_SPACEDIM > 2
			else if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				/* Solving for end point. No need to consider triple end point.
				   They have already been solved. */
				//std::cout << "End point: Neumann boundary in back face. i,j,k = "<< i << "," << j << "," << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				//std::cout << "End point: Neumann boundary in front face. i,j,k = "<< i << "," << j << "," << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
#endif
			//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
			if(debug)
			{
				if(m == point_left)
				{
					std::cout << "Looking at point on the left (" << point_left[0] << "," << point_left[1] << "," << point_left[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m+dx, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				if(m == point_left)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			
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
	if(bc_hi_str[0] == "traction" && box.hiVect()[0]==domain.hiVect()[0])
	{
		int i = box.hiVect()[0] + 1;
		//std::cout << "Neumann BC on the right face. i = " << i << std::endl;
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
			if (j == domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				/* Solving for end points */
				//std::cout << "End point: Neumann BC on the bottom face. j = " << j << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
#if AMREX_SPACEDIM > 2
				/* Solving for triple end point */ 
				if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
				{
					//std::cout << "Triple end point. Neumann BC on the back face. k = " << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
				{
					//std::cout << "Triple end point. Neumann BC on the front face. k = " << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
			}
			else if (j == domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				//std::cout << "End point. Neumann BC on the top face. j = " << j << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
#if AMREX_SPACEDIM > 2
				if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
				{
					//std::cout << "Triple end point. Neumann BC on the back face. k = " << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
					points.push_back(5);
				}
				else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
				{
					//std::cout << "Triple end point. Neumann BC on the front face. k = " << k << std::endl;
					traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
					points.push_back(6);
				}
#endif
			}
#if AMREX_SPACEDIM > 2
			else if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the back face. k = " << k << std::endl;
				/* Solving for end point. No need to consider triple end point.
				   They have already been solved. */
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the front face. k = " << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
#endif
			//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
			if(debug)
			{
				if(m == point_right)
				{
					std::cout << "Looking at point on the right (" << point_right[0] << "," << point_right[1] << "," << point_right[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m-dx, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				if(m == point_right)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}

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
	if(bc_lo_str[1] == "traction" && box.loVect()[1]==domain.loVect()[1])
	{
		int j = box.loVect()[1] - 1;
		//std::cout << "Neumann BC on bottom face. j = " << j << std::endl;
		AMREX_D_TERM(	for (int i = box.loVect()[0]; i <= box.hiVect()[0]; i++),
						,
						for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
		{
			if(i == domain.loVect()[0] && bc_lo_str[0] == "traction") continue;
			if(i == domain.hiVect()[0] && bc_hi_str[0] == "traction") continue;
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
			if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the back face. k = " << k << std::endl;
				/* Solving for end points */
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the front face. k = " << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
			//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
			if(debug)
			{
				if(m == point_bottom)
				{
					std::cout << "Looking at point on the bottom (" << point_bottom[0] << "," << point_bottom[1] << "," << point_bottom[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m+dy, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				//std::cout << "Looking at point on the bottom (" << point_bottom[0] << "," << point_bottom[1] << "," << point_bottom[2] <<")" <<std::endl; 
				if(m == point_bottom)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}

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
	if(bc_hi_str[1] == "traction" && box.hiVect()[1]==domain.hiVect()[1])
	{
		// non-end ghost cells. End ghost cells will be treated separately.
		int j = box.hiVect()[1] + 1;
		//std::cout << "Neumann BC on the top face. j = " << j << std::endl;
		AMREX_D_TERM(	for (int i = box.loVect()[0]; i <= box.hiVect()[0]; i++),
				,
				for (int k = box.loVect()[2]; k <= box.hiVect()[2]; k++))
		{
			if(i == domain.loVect()[0] && bc_lo_str[0] == "traction") continue;
			if(i == domain.hiVect()[0] && bc_hi_str[0] == "traction") continue;
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
			if (k == domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the back face. k = " << k << std::endl;
				/* Solving for end points */
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k == domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				//std::cout << "End point: Neumann BC on the front face. k = " << k << std::endl;
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
			//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
			if(debug)
			{
				if(m == point_top)
				{
					std::cout << "Looking at point on the top (" << point_top[0] << "," << point_top[1] << "," << point_top[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m-dy, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				//std::cout << "Looking at point on the top (" << point_top[0] << "," << point_top[1] << "," << point_top[2] <<")" <<std::endl; 
				if(m == point_top)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}

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
	if(bc_lo_str[2] == "traction" && box.loVect()[2]==domain.loVect()[2])
	{
		// End and triple end cells have already been filled - so no need to iterate over those
		int k = box.loVect()[2] - 1;
		//std::cout << "Neumann BC on back face. k = " << k << std::endl;
		AMREX_D_TERM(	for (int i = box.loVect()[0]; i <= box.hiVect()[0]; i++),
				for (int j = box.loVect()[1]; j <= box.hiVect()[1]; j++),
				)
		{
			if(i == domain.loVect()[0] && bc_lo_str[0] == "traction") continue;
			if(i == domain.hiVect()[0] && bc_hi_str[0] == "traction") continue;
			if(j == domain.loVect()[1] && bc_lo_str[1] == "traction") continue;
			if(j == domain.hiVect()[1] && bc_hi_str[1] == "traction") continue;
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
			if(debug)
			{
				if(m == point_back)
				{
					std::cout << "Looking at point on the back (" << point_back[0] << "," << point_back[1] << "," << point_back[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m+dz, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				//std::cout << "Looking at point on the back (" << point_back[0] << "," << point_back[1] << "," << point_back[2] <<")" <<std::endl; 
				if(m == point_back)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}

			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[5](n);
		}
	}

	/* Step 6: Fill front face of ghost cells */
	if(bc_hi_str[2] == "traction" && box.hiVect()[2]==domain.hiVect()[2])
	{
		int k = box.hiVect()[2] + 1;
		//std::cout << "Neumann BC on the front face. k = " << k << std::endl;
		AMREX_D_TERM(	for (int i = box.loVect()[0] + 1; i <= box.hiVect()[0] - 1; i++),
						for (int j = box.loVect()[1] + 1; j <= box.hiVect()[1] - 1; j++),
					)
		{
			if(i == domain.loVect()[0] && bc_lo_str[0] == "traction") continue;
			if(i == domain.hiVect()[0] && bc_hi_str[0] == "traction") continue;
			if(j == domain.loVect()[1] && bc_lo_str[1] == "traction") continue;
			if(j == domain.hiVect()[1] && bc_hi_str[1] == "traction") continue;
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
			//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
			if(debug)
			{
				if(m == point_front)
				{
					std::cout << "Looking at point on the front (" << point_front[0] << "," << point_front[1] << "," << point_front[2] <<")" <<std::endl; 
					for (int p = 0; p < 7; p++)
						std::cout << "Before state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}
			StencilFill(stencil, traction, points, m-dz, m_amrlev, m_mglev, mfi, debug);
			if(debug)
			{
				//std::cout << "Looking at point on the back (" << point_back[0] << "," << point_back[1] << "," << point_back[2] <<")" <<std::endl; 
				if(m == point_front)
				{
					for (int p = 0; p < 7; p++)
						std::cout << "After state p = " << p+1 << ". (" << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << ")" << std::endl;
				}
			}

			for (int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[6](n);
		}
	}
#endif
#endif
	/*Finally let's do averaging of corner cells - edge corners and triple corners*/
	AMREX_D_TERM(	for(int i = box.loVect()[0]-ngrow; i <= box.hiVect()[0]+ngrow; i++),
					for(int j = box.loVect()[1]-ngrow; j <= box.hiVect()[1]+ngrow; j++),
					for(int k = box.loVect()[2]-ngrow; k <= box.hiVect()[2]+ngrow; k++))
	{
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		int mul1 = 0, mul2 = 0, mul3 = 0;

#if AMREX_SPACEDIM > 1
		if(i == domain.loVect()[0]-1 && j == domain.loVect()[1]-1)
		{
			mul1 = 1; mul2 = 1;
#if AMREX_SPACEDIM > 2
			if(k==domain.loVect()[2]-1) mul3 = 1;
			else if(k==domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
		}
		else if(i == domain.loVect()[0]-1 && j == domain.hiVect()[1]+1)
		{
			mul1 = 1; mul2 = -1;
#if AMREX_SPACEDIM > 2
			if(k==domain.loVect()[2]-1) mul3 = 1;
			else if(k==domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
		}
		else if(i == domain.hiVect()[0]+1 && j == domain.loVect()[1]-1)
		{
			mul1 = -1; mul2 = 1;
#if AMREX_SPACEDIM > 2
			if(k==domain.loVect()[2]-1) mul3 = 1;
			else if(k==domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
		}
		else if(i == domain.hiVect()[0]+1 && j == domain.hiVect()[1]+1)
		{
			mul1 = -1; mul2 = -1;
#if AMREX_SPACEDIM > 2
			if(k==domain.loVect()[2]-1) mul3 = 1;
			else if(k==domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
		}
#if AMREX_SPACEDIM > 2
		else if (i == domain.loVect()[0]-1 && k == domain.loVect()[2]-1){mul1 = 1; mul2 = 0; mul3 = 1;}
		else if (i == domain.loVect()[0]-1 && k == domain.hiVect()[2]+1){mul1 = 1; mul2 = 0; mul3 = -1;}
		else if (i == domain.hiVect()[0]+1 && k == domain.loVect()[2]-1){mul1 = -1; mul2 = 0; mul3 = 1;}
		else if (i == domain.hiVect()[0]+1 && k == domain.hiVect()[2]+1){mul1 = -1; mul2 = 0; mul3 = -1;}
		else if (j == domain.loVect()[1]-1 && k == domain.loVect()[2]-1){mul1 = 0; mul2 = 1; mul3 = 1;}
		else if (j == domain.loVect()[1]-1 && k == domain.hiVect()[2]+1){mul1 = 0; mul2 = 1; mul3 = -1;}
		else if (j == domain.hiVect()[1]+1 && k == domain.loVect()[2]-1){mul1 = 0; mul2 = -1; mul3 = 1;}
		else if (j == domain.hiVect()[1]+1 && k == domain.hiVect()[2]+1){mul1 = 0; mul2 = -1; mul3 = -1;}
#endif
		if (m == amrex::IntVect(-1,3,3))
			std::cout << "mul1 = " << mul1 << ". mul2 = " << mul2 << ". mul3 = " << mul3 << std::endl;

		if(AMREX_D_TERM(std::abs(mul1), + std::abs(mul2), + std::abs(mul3)) == 0) continue;

		for(int n = 0; n<AMREX_SPACEDIM; n++)
				mf_box(m,n) = (AMREX_D_TERM(mf_box(m+mul1*dx,n),+mf_box(m+mul2*dy,n),+mf_box(m+mul3*dz,n)))/(AMREX_D_TERM(std::abs(mul1),+std::abs(mul2),+std::abs(mul3)));
	}
#endif
}

amrex::BCRec
Elastic::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
Elastic::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(bc_lo[0]),BCUtil::IsPeriodic(bc_lo[1]),BCUtil::IsPeriodic(bc_lo[2]))};
}

amrex::Periodicity 
Elastic::Periodicity () const
{
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(bc_lo[0]),
							      m_geom.Domain().length(1) * BCUtil::IsPeriodic(bc_lo[1]),
							      m_geom.Domain().length(2) * BCUtil::IsPeriodic(bc_lo[2]))));
}
amrex::Periodicity 
Elastic::Periodicity (const amrex::Box& b) {
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(bc_lo[0]),
							      b.length(1) * BCUtil::IsPeriodic(bc_lo[1]),
							      b.length(2) * BCUtil::IsPeriodic(bc_lo[2]))));

}

// Stencil Fill routine - takes in a stencil, a list of unknown points and fills the unknown values
// in the stencil.
//#define m_operator->C(i,j,k,l,m,amrlev,mglev,mfi) m_operator->m_operator->C(i,j,k,l,m,amrlev,mglev,mfi)
void 
Elastic::StencilFill(	amrex::Vector<Set::Vector> &stencil,
			const amrex::Vector<Set::Vector> &traction,
			const amrex::Vector<int> &points,
			const amrex::IntVect &m,
			const int amrlev,
			const int mglev,
			const amrex::MFIter &mfi,
			const bool debug)
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

	//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
	if(points.size() > 3 || points.size() < 1)
		Util::Abort("Number of unknown points can not be greater than 3");

	if(points.size() != traction.size())
		Util::Abort("Mismatch between number of unknown points and tractions");
	//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
	const Real* DX = m_geom.CellSize();
	//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
  
#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	if(points.size() == 1)		//Need three equations - non corner cases
	{
		//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
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

			AMREX_D_TERM( 	left(0,0) = m_operator->C(0,0,0,0,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = m_operator->C(0,0,1,0,m,amrlev,mglev,mfi);
					left(1,0) = m_operator->C(1,0,0,0,m,amrlev,mglev,mfi);
					left(1,1) = m_operator->C(1,0,1,0,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = m_operator->C(0,0,2,0,m,amrlev,mglev,mfi);
					left(1,2) = m_operator->C(1,0,2,0,m,amrlev,mglev,mfi);
					left(2,0) = m_operator->C(2,0,0,0,m,amrlev,mglev,mfi);
					left(2,1) = m_operator->C(2,0,1,0,m,amrlev,mglev,mfi);
					left(2,2) = m_operator->C(2,0,2,0,m,amrlev,mglev,mfi););
			AMREX_D_TERM(	right(0) = AMREX_D_TERM(mul*traction[0](0)
								, 
								- m_operator->C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- m_operator->C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								);
					, // 2D
					right(1) = AMREX_D_TERM(mul*traction[0](1)
								, 
								- m_operator->C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- m_operator->C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								);
					, // 3D
					right(2) = AMREX_D_TERM(mul*traction[0](2)
								, 
								- m_operator->C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								,
								- m_operator->C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								););
			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);
					
			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test;
				test(0) = 	m_operator->C(0,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(1) = 	m_operator->C(1,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(2) = 	m_operator->C(2,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);

				if(m == amrex::IntVect(0,7,5))
				{
					std::cout << "mul = " << mul << ". sol = " << sol(0) << ", " << sol(1) << ", " << sol(2) << std::endl;
				}
				if((test-traction[0]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
				}
			}
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
			AMREX_D_TERM( 	left(0,0) = m_operator->C(0,1,0,1,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = m_operator->C(0,1,1,1,m,amrlev,mglev,mfi);
					left(1,0) = m_operator->C(1,1,0,1,m,amrlev,mglev,mfi);
					left(1,1) = m_operator->C(1,1,1,1,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = m_operator->C(0,1,2,1,m,amrlev,mglev,mfi);
					left(1,2) = m_operator->C(1,1,2,1,m,amrlev,mglev,mfi);
					left(2,0) = m_operator->C(2,1,0,1,m,amrlev,mglev,mfi);
					left(2,1) = m_operator->C(2,1,1,1,m,amrlev,mglev,mfi);
					left(2,2) = m_operator->C(2,1,2,1,m,amrlev,mglev,mfi););
			AMREX_D_TERM(	right(0) = AMREX_D_TERM(mul*traction[0](0)
								, 
								- m_operator->C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 2D
					right(1) = AMREX_D_TERM(mul*traction[0](1)
								, 
								- m_operator->C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 3D
					right(2) = AMREX_D_TERM(mul*traction[0](2)
								, 
								- m_operator->C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
								- m_operator->C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
								- m_operator->C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
								- m_operator->C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								););
			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);

			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test;
				test(0) = 	m_operator->C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(0,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(1) = 	m_operator->C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(1,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(2) = 	m_operator->C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(2,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);

				if((test-traction[0]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
				}
			}
			
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
			AMREX_D_TERM( 	left(0,0) = m_operator->C(0,2,0,2,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = m_operator->C(0,2,1,2,m,amrlev,mglev,mfi);
					left(1,0) = m_operator->C(1,2,0,2,m,amrlev,mglev,mfi);
					left(1,1) = m_operator->C(1,2,1,2,m,amrlev,mglev,mfi);
					, // 3D
					left(0,2) = m_operator->C(0,2,2,2,m,amrlev,mglev,mfi);
					left(1,2) = m_operator->C(1,2,2,2,m,amrlev,mglev,mfi);
					left(2,0) = m_operator->C(2,2,0,2,m,amrlev,mglev,mfi);
					left(2,1) = m_operator->C(2,2,1,2,m,amrlev,mglev,mfi);
					left(2,2) = m_operator->C(2,2,2,2,m,amrlev,mglev,mfi););
			AMREX_D_TERM(	right(0) = AMREX_D_TERM(mul*traction[0](0)
								, 
								- m_operator->C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								- m_operator->C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								- m_operator->C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 2D
					right(1) = AMREX_D_TERM(mul*traction[0](1)
								, 
								- m_operator->C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								- m_operator->C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								- m_operator->C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					, // 3D
					right(2) = AMREX_D_TERM(mul*traction[0](2)
								, 
								- m_operator->C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
								- m_operator->C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
								,
								- m_operator->C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
								- m_operator->C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
								- m_operator->C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
								- m_operator->C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
								);
					);
			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);

			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test;
				test(0) = 	m_operator->C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(0,2,0,2,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,2,1,2,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,2,2,2,m,amrlev,mglev,mfi)*sol(2);
				test(1) = 	m_operator->C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(1,2,0,2,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,2,1,2,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,2,2,2,m,amrlev,mglev,mfi)*sol(2);
				test(2) = 	m_operator->C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(2,2,0,2,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,2,1,2,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,2,2,2,m,amrlev,mglev,mfi)*sol(2);

				if((test-traction[0]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
				}
			}
					
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

			AMREX_D_TERM(	left(0,0) = m_operator->C(0,0,0,0,m,amrlev,mglev,mfi);
					, // 2D
					left(0,1) = m_operator->C(0,0,1,0,m,amrlev,mglev,mfi);
					left(0,2) = m_operator->C(0,0,0,1,m,amrlev,mglev,mfi);
					left(0,3) = m_operator->C(0,0,1,1,m,amrlev,mglev,mfi);
					left(1,0) = m_operator->C(1,0,0,0,m,amrlev,mglev,mfi);
					left(1,1) = m_operator->C(1,0,1,0,m,amrlev,mglev,mfi);
					left(1,2) = m_operator->C(1,0,0,1,m,amrlev,mglev,mfi);
					left(1,3) = m_operator->C(1,0,1,1,m,amrlev,mglev,mfi);
					left(2,0) = m_operator->C(0,1,0,0,m,amrlev,mglev,mfi);
					left(2,1) = m_operator->C(0,1,1,0,m,amrlev,mglev,mfi);
					left(2,2) = m_operator->C(0,1,0,1,m,amrlev,mglev,mfi);
					left(2,3) = m_operator->C(0,1,1,1,m,amrlev,mglev,mfi);
					left(3,0) = m_operator->C(1,1,0,0,m,amrlev,mglev,mfi);
					left(3,1) = m_operator->C(1,1,1,0,m,amrlev,mglev,mfi);
					left(3,2) = m_operator->C(1,1,0,1,m,amrlev,mglev,mfi);
					left(3,3) = m_operator->C(1,1,1,1,m,amrlev,mglev,mfi);
					, // 3D
					left(0,4) = m_operator->C(0,0,2,0,m,amrlev,mglev,mfi);
					left(0,5) = m_operator->C(0,0,2,1,m,amrlev,mglev,mfi);
					left(1,4) = m_operator->C(1,0,2,0,m,amrlev,mglev,mfi);
					left(1,5) = m_operator->C(1,0,2,1,m,amrlev,mglev,mfi);
					left(2,4) = m_operator->C(0,1,2,0,m,amrlev,mglev,mfi);
					left(2,5) = m_operator->C(0,1,2,1,m,amrlev,mglev,mfi);
					left(3,4) = m_operator->C(1,1,2,0,m,amrlev,mglev,mfi);
					left(3,5) = m_operator->C(1,1,2,1,m,amrlev,mglev,mfi);
					left(4,0) = m_operator->C(2,0,0,0,m,amrlev,mglev,mfi);
					left(4,1) = m_operator->C(2,0,1,0,m,amrlev,mglev,mfi);
					left(4,2) = m_operator->C(2,0,0,1,m,amrlev,mglev,mfi);
					left(4,3) = m_operator->C(2,0,1,1,m,amrlev,mglev,mfi);
					left(4,4) = m_operator->C(2,0,2,0,m,amrlev,mglev,mfi);
					left(4,5) = m_operator->C(2,0,2,1,m,amrlev,mglev,mfi);
					left(5,0) = m_operator->C(2,1,0,0,m,amrlev,mglev,mfi);
					left(5,1) = m_operator->C(2,1,1,0,m,amrlev,mglev,mfi);
					left(5,2) = m_operator->C(2,1,0,1,m,amrlev,mglev,mfi);
					left(5,3) = m_operator->C(2,1,1,1,m,amrlev,mglev,mfi);
					left(5,4) = m_operator->C(2,1,2,0,m,amrlev,mglev,mfi);
					left(5,5) = m_operator->C(2,1,2,1,m,amrlev,mglev,mfi););

			AMREX_D_TERM(	right(0) = AMREX_D_TERM(	mul1*traction[0](0)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									);
					, // 2D
					right(1) = AMREX_D_TERM(	mul1*traction[0](1)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									);
					right(2) = AMREX_D_TERM(	mul2*traction[1](0)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									);
					right(3) = AMREX_D_TERM(	mul2*traction[1](1)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									);
					, // 3D
					right(4) = AMREX_D_TERM(	mul1*traction[0](2)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									);
					right(5) = AMREX_D_TERM(	mul2*traction[1](2)
									, // 2D
									+ 0.0
									, // 3D
									- m_operator->C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
									- m_operator->C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
									- m_operator->C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2)
									););

			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);

			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test, test2;
				test(0) = 	m_operator->C(0,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,0,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,0,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,0,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,0,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(0,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(0,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(0,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(1) = 	m_operator->C(1,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,0,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,0,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,0,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,0,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(1,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(1,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(1,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test(2) = 	m_operator->C(2,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,0,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,0,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,0,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,0,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(2,0,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(2,0,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(2,0,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test2(0) = 	m_operator->C(0,1,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,1,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,1,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,1,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,1,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,1,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(0,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(0,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(0,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test2(1) = 	m_operator->C(1,1,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,1,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,1,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,1,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,1,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,1,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(1,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(1,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(1,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);
				test2(2) = 	m_operator->C(2,1,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,1,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,1,2,0,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,1,0,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,1,1,1,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,1,2,1,m,amrlev,mglev,mfi)*sol(5)
							+ m_operator->C(2,1,0,2,m,amrlev,mglev,mfi)*gradu_3(0)
							+ m_operator->C(2,1,1,2,m,amrlev,mglev,mfi)*gradu_3(1)
							+ m_operator->C(2,1,2,2,m,amrlev,mglev,mfi)*gradu_3(2);

				if((test-traction[0]).norm() > 1.e-2 || (test2-traction[1]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
					std::cout << "test2 =" << std::endl << test2 << std::endl;
					std::cout << "traction2 = " << std::endl << traction[1] << std::endl;
				}
			}
			AMREX_D_TERM(	stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
					,
					stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
					stencil[points[1]](0) = stencil[0](0) + mul2*DX[1]*sol(2);
					stencil[points[1]](1) = stencil[0](1) + mul2*DX[1]*sol(3);
					,
					stencil[points[0]](2) = stencil[0](2) + mul1*DX[0]*sol(4);
					stencil[points[1]](2) = stencil[0](2) + mul2*DX[1]*sol(5););
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

			left(0,0) = m_operator->C(0,0,0,0,m,amrlev,mglev,mfi);
			left(0,1) = m_operator->C(0,0,1,0,m,amrlev,mglev,mfi);
			left(0,2) = m_operator->C(0,0,2,0,m,amrlev,mglev,mfi);
			left(0,3) = m_operator->C(0,0,0,2,m,amrlev,mglev,mfi);
			left(0,4) = m_operator->C(0,0,1,2,m,amrlev,mglev,mfi);
			left(0,5) = m_operator->C(0,0,2,2,m,amrlev,mglev,mfi);
			left(1,0) = m_operator->C(1,0,0,0,m,amrlev,mglev,mfi);
			left(1,1) = m_operator->C(1,0,1,0,m,amrlev,mglev,mfi);
			left(1,2) = m_operator->C(1,0,2,0,m,amrlev,mglev,mfi);
			left(1,3) = m_operator->C(1,0,0,2,m,amrlev,mglev,mfi);
			left(1,4) = m_operator->C(1,0,1,2,m,amrlev,mglev,mfi);
			left(1,5) = m_operator->C(1,0,2,2,m,amrlev,mglev,mfi);
			left(2,0) = m_operator->C(2,0,0,0,m,amrlev,mglev,mfi);
			left(2,1) = m_operator->C(2,0,1,0,m,amrlev,mglev,mfi);
			left(2,2) = m_operator->C(2,0,2,0,m,amrlev,mglev,mfi);
			left(2,3) = m_operator->C(2,0,0,2,m,amrlev,mglev,mfi);
			left(2,4) = m_operator->C(2,0,1,2,m,amrlev,mglev,mfi);
			left(2,5) = m_operator->C(2,0,2,2,m,amrlev,mglev,mfi);
			left(3,0) = m_operator->C(0,2,0,0,m,amrlev,mglev,mfi);
			left(3,1) = m_operator->C(0,2,1,0,m,amrlev,mglev,mfi);
			left(3,2) = m_operator->C(0,2,2,0,m,amrlev,mglev,mfi);
			left(3,3) = m_operator->C(0,2,0,2,m,amrlev,mglev,mfi);
			left(3,4) = m_operator->C(0,2,1,2,m,amrlev,mglev,mfi);
			left(3,5) = m_operator->C(0,2,2,2,m,amrlev,mglev,mfi);
			left(4,0) = m_operator->C(1,2,0,0,m,amrlev,mglev,mfi);
			left(4,1) = m_operator->C(1,2,1,0,m,amrlev,mglev,mfi);
			left(4,2) = m_operator->C(1,2,2,0,m,amrlev,mglev,mfi);
			left(4,3) = m_operator->C(1,2,0,2,m,amrlev,mglev,mfi);
			left(4,4) = m_operator->C(1,2,1,2,m,amrlev,mglev,mfi);
			left(4,5) = m_operator->C(1,2,2,2,m,amrlev,mglev,mfi);
			left(5,0) = m_operator->C(2,2,0,0,m,amrlev,mglev,mfi);
			left(5,1) = m_operator->C(2,2,1,0,m,amrlev,mglev,mfi);
			left(5,2) = m_operator->C(2,2,2,0,m,amrlev,mglev,mfi);
			left(5,3) = m_operator->C(2,2,0,2,m,amrlev,mglev,mfi);
			left(5,4) = m_operator->C(2,2,1,2,m,amrlev,mglev,mfi);
			left(5,5) = m_operator->C(2,2,2,2,m,amrlev,mglev,mfi);

			right(0) = 	mul1*traction[0](0)
				- m_operator->C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(1) = 	mul1*traction[0](1)
				- m_operator->C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(2) = 	mul1*traction[0](2)
				- m_operator->C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(3) = 	mul2*traction[1](0)
				- m_operator->C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(4) = 	mul2*traction[1](1)
				- m_operator->C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);
			right(5) = 	mul2*traction[1](2)
				- m_operator->C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
				- m_operator->C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
				- m_operator->C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2);

			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);

			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test, test2;
				test(0) = 	m_operator->C(0,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(0,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(0,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(0,0,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,0,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,0,2,2,m,amrlev,mglev,mfi)*sol(5);
				test(1) = 	m_operator->C(1,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(1,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(1,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(1,0,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,0,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,0,2,2,m,amrlev,mglev,mfi)*sol(5);
				test(2) = 	m_operator->C(2,0,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,0,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,0,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,0,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(2,0,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(2,0,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(2,0,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,0,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,0,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(0) = 	m_operator->C(0,2,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,2,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,2,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(0,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(0,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(0,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,2,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(1) = 	m_operator->C(1,2,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,2,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,2,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(1,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(1,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(1,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,2,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(2) = 	m_operator->C(2,2,0,0,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,2,1,0,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,2,2,0,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,2,0,1,m,amrlev,mglev,mfi)*gradu_2(0)
							+ m_operator->C(2,2,1,1,m,amrlev,mglev,mfi)*gradu_2(1)
							+ m_operator->C(2,2,2,1,m,amrlev,mglev,mfi)*gradu_2(2)
							+ m_operator->C(2,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,2,2,2,m,amrlev,mglev,mfi)*sol(5);

				if((test-traction[0]).norm() > 1.e-2 || (test2-traction[1]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
					std::cout << "test2 =" << std::endl << test2 << std::endl;
					std::cout << "traction2 = " << std::endl << traction[1] << std::endl;
				}
			}
			
			stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
			stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
			stencil[points[0]](2) = stencil[0](2) + mul1*DX[0]*sol(2);

			stencil[points[1]](0) = stencil[0](0) + mul2*DX[2]*sol(3);
			stencil[points[1]](1) = stencil[0](1) + mul2*DX[2]*sol(4);
			stencil[points[1]](2) = stencil[0](2) + mul2*DX[2]*sol(5);
		}

		else if ((points[0] == 3 || points[0] == 4) && (points[1] == 5 || points[1] == 6)) 
		{
			mul1 = points[0] == 3 ? -1 : 1;
			mul2 = points[1] == 5 ? -1 : 1;

			//gradu_k,1
			Set::Vector gradu_1;
			gradu_1(0) = (stencil[2](0) - stencil[1](0))/(2.0*DX[0]);
			gradu_1(1) = (stencil[2](1) - stencil[1](1))/(2.0*DX[0]);
			gradu_1(2) = (stencil[2](2) - stencil[1](2))/(2.0*DX[0]);

			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,2*AMREX_SPACEDIM> left;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> right;
			Eigen::Matrix<amrex::Real,2*AMREX_SPACEDIM,1> sol; 

			left(0,0) = m_operator->C(0,1,0,1,m,amrlev,mglev,mfi);
			left(0,1) = m_operator->C(0,1,1,1,m,amrlev,mglev,mfi);
			left(0,2) = m_operator->C(0,1,2,1,m,amrlev,mglev,mfi);
			left(0,3) = m_operator->C(0,1,0,2,m,amrlev,mglev,mfi);
			left(0,4) = m_operator->C(0,1,1,2,m,amrlev,mglev,mfi);
			left(0,5) = m_operator->C(0,1,2,2,m,amrlev,mglev,mfi);
			left(1,0) = m_operator->C(1,1,0,1,m,amrlev,mglev,mfi);
			left(1,1) = m_operator->C(1,1,1,1,m,amrlev,mglev,mfi);
			left(1,2) = m_operator->C(1,1,2,1,m,amrlev,mglev,mfi);
			left(1,3) = m_operator->C(1,1,0,2,m,amrlev,mglev,mfi);
			left(1,4) = m_operator->C(1,1,1,2,m,amrlev,mglev,mfi);
			left(1,5) = m_operator->C(1,1,2,2,m,amrlev,mglev,mfi);
			left(2,0) = m_operator->C(2,1,0,1,m,amrlev,mglev,mfi);
			left(2,1) = m_operator->C(2,1,1,1,m,amrlev,mglev,mfi);
			left(2,2) = m_operator->C(2,1,2,1,m,amrlev,mglev,mfi);
			left(2,3) = m_operator->C(2,1,0,2,m,amrlev,mglev,mfi);
			left(2,4) = m_operator->C(2,1,1,2,m,amrlev,mglev,mfi);
			left(2,5) = m_operator->C(2,1,2,2,m,amrlev,mglev,mfi);
			left(3,0) = m_operator->C(0,2,0,1,m,amrlev,mglev,mfi);
			left(3,1) = m_operator->C(0,2,1,1,m,amrlev,mglev,mfi);
			left(3,2) = m_operator->C(0,2,2,1,m,amrlev,mglev,mfi);
			left(3,3) = m_operator->C(0,2,0,2,m,amrlev,mglev,mfi);
			left(3,4) = m_operator->C(0,2,1,2,m,amrlev,mglev,mfi);
			left(3,5) = m_operator->C(0,2,2,2,m,amrlev,mglev,mfi);
			left(4,0) = m_operator->C(1,2,0,1,m,amrlev,mglev,mfi);
			left(4,1) = m_operator->C(1,2,1,1,m,amrlev,mglev,mfi);
			left(4,2) = m_operator->C(1,2,2,1,m,amrlev,mglev,mfi);
			left(4,3) = m_operator->C(1,2,0,2,m,amrlev,mglev,mfi);
			left(4,4) = m_operator->C(1,2,1,2,m,amrlev,mglev,mfi);
			left(4,5) = m_operator->C(1,2,2,2,m,amrlev,mglev,mfi);
			left(5,0) = m_operator->C(2,2,0,1,m,amrlev,mglev,mfi);
			left(5,1) = m_operator->C(2,2,1,1,m,amrlev,mglev,mfi);
			left(5,2) = m_operator->C(2,2,2,1,m,amrlev,mglev,mfi);
			left(5,3) = m_operator->C(2,2,0,2,m,amrlev,mglev,mfi);
			left(5,4) = m_operator->C(2,2,1,2,m,amrlev,mglev,mfi);
			left(5,5) = m_operator->C(2,2,2,2,m,amrlev,mglev,mfi);

			right(0) = 	mul1*traction[0](0)
				- m_operator->C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(1) = 	mul1*traction[0](1)
				- m_operator->C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(2) = 	mul1*traction[0](2)
				- m_operator->C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(3) = 	mul2*traction[1](0)
				- m_operator->C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(4) = 	mul2*traction[1](1)
				- m_operator->C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);
			right(5) = 	mul2*traction[1](2)
				- m_operator->C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
				- m_operator->C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
				- m_operator->C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2);

			//sol = left.ldlt().solve(right); // we can change this solver as needed
			sol = left.colPivHouseholderQr().solve(right);

			if (debug)
			{
				//std::cout << "DX = " << DX[0] << "," << DX[1] << "," << DX[2] << std::endl;
				//std::cout << "gradu_2 = " << std::endl << gradu_2 << std::endl;
				//std::cout << "gradu_3 = " << std::endl << gradu_3 << std::endl;
				//std::cout << "left = " << std::endl << left << std::endl;
				//std::cout << "right = " << std::endl << right << std::endl;
				//std::cout << "sol = " << std::endl << sol << std::endl;

				Set::Vector test, test2;
				test(0) = 	m_operator->C(0,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(0,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(0,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(0,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,1,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,1,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,1,2,2,m,amrlev,mglev,mfi)*sol(5);
				test(1) = 	m_operator->C(1,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(1,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(1,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(1,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,1,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,1,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,1,2,2,m,amrlev,mglev,mfi)*sol(5);
				test(2) = 	m_operator->C(2,1,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(2,1,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(2,1,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(2,1,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,1,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,1,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,1,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,1,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,1,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(0) = 	m_operator->C(0,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(0,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(0,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(0,2,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(0,2,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(0,2,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(0,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(0,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(0,2,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(1) = 	m_operator->C(1,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(1,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(1,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(1,2,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(1,2,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(1,2,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(1,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(1,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(1,2,2,2,m,amrlev,mglev,mfi)*sol(5);
				test2(2) = 	m_operator->C(2,2,0,0,m,amrlev,mglev,mfi)*gradu_1(0)
							+ m_operator->C(2,2,1,0,m,amrlev,mglev,mfi)*gradu_1(1)
							+ m_operator->C(2,2,2,0,m,amrlev,mglev,mfi)*gradu_1(2)
							+ m_operator->C(2,2,0,1,m,amrlev,mglev,mfi)*sol(0)
							+ m_operator->C(2,2,1,1,m,amrlev,mglev,mfi)*sol(1)
							+ m_operator->C(2,2,2,1,m,amrlev,mglev,mfi)*sol(2)
							+ m_operator->C(2,2,0,2,m,amrlev,mglev,mfi)*sol(3)
							+ m_operator->C(2,2,1,2,m,amrlev,mglev,mfi)*sol(4)
							+ m_operator->C(2,2,2,2,m,amrlev,mglev,mfi)*sol(5);

				if((test-traction[0]).norm() > 1.e-2 || (test2-traction[1]).norm() > 1.e-2)
				{
					std::cout << "test =" << std::endl << test << std::endl;
					std::cout << "traction = " << std::endl << traction[0] << std::endl;
					std::cout << "test2 =" << std::endl << test2 << std::endl;
					std::cout << "traction2 = " << std::endl << traction[1] << std::endl;
				}
			}
			
			stencil[points[0]](0) = stencil[0](0) + mul1*DX[1]*sol(0);
			stencil[points[0]](1) = stencil[0](1) + mul1*DX[1]*sol(1);
			stencil[points[0]](2) = stencil[0](2) + mul1*DX[1]*sol(2);

			stencil[points[1]](0) = stencil[0](0) + mul2*DX[2]*sol(3);
			stencil[points[1]](1) = stencil[0](1) + mul2*DX[2]*sol(4);
			stencil[points[1]](2) = stencil[0](2) + mul2*DX[2]*sol(5);
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

		left(0,0) = m_operator->C(0,0,0,0,m,amrlev,mglev,mfi);
		left(0,1) = m_operator->C(0,0,1,0,m,amrlev,mglev,mfi);
		left(0,2) = m_operator->C(0,0,2,0,m,amrlev,mglev,mfi);
		left(0,3) = m_operator->C(0,0,0,1,m,amrlev,mglev,mfi);
		left(0,4) = m_operator->C(0,0,1,1,m,amrlev,mglev,mfi);
		left(0,5) = m_operator->C(0,0,2,1,m,amrlev,mglev,mfi);
		left(0,6) = m_operator->C(0,0,0,2,m,amrlev,mglev,mfi);
		left(0,7) = m_operator->C(0,0,2,2,m,amrlev,mglev,mfi);
		left(0,8) = m_operator->C(0,0,2,2,m,amrlev,mglev,mfi);
		left(1,0) = m_operator->C(1,0,0,0,m,amrlev,mglev,mfi);
		left(1,1) = m_operator->C(1,0,1,0,m,amrlev,mglev,mfi);
		left(1,2) = m_operator->C(1,0,2,0,m,amrlev,mglev,mfi);
		left(1,3) = m_operator->C(1,0,0,1,m,amrlev,mglev,mfi);
		left(1,4) = m_operator->C(1,0,1,1,m,amrlev,mglev,mfi);
		left(1,5) = m_operator->C(1,0,2,1,m,amrlev,mglev,mfi);
		left(1,6) = m_operator->C(1,0,0,2,m,amrlev,mglev,mfi);
		left(1,7) = m_operator->C(1,0,2,2,m,amrlev,mglev,mfi);
		left(1,8) = m_operator->C(1,0,2,2,m,amrlev,mglev,mfi);
		left(2,0) = m_operator->C(2,0,0,0,m,amrlev,mglev,mfi);
		left(2,1) = m_operator->C(2,0,1,0,m,amrlev,mglev,mfi);
		left(2,2) = m_operator->C(2,0,2,0,m,amrlev,mglev,mfi);
		left(2,3) = m_operator->C(2,0,0,1,m,amrlev,mglev,mfi);
		left(2,4) = m_operator->C(2,0,1,1,m,amrlev,mglev,mfi);
		left(2,5) = m_operator->C(2,0,2,1,m,amrlev,mglev,mfi);
		left(2,6) = m_operator->C(2,0,0,2,m,amrlev,mglev,mfi);
		left(2,7) = m_operator->C(2,0,1,2,m,amrlev,mglev,mfi);
		left(2,8) = m_operator->C(2,0,2,2,m,amrlev,mglev,mfi);
		left(3,0) = m_operator->C(0,1,0,0,m,amrlev,mglev,mfi);
		left(3,1) = m_operator->C(0,1,1,0,m,amrlev,mglev,mfi);
		left(3,2) = m_operator->C(0,1,2,0,m,amrlev,mglev,mfi);
		left(3,3) = m_operator->C(0,1,0,1,m,amrlev,mglev,mfi);
		left(3,4) = m_operator->C(0,1,1,1,m,amrlev,mglev,mfi);
		left(3,5) = m_operator->C(0,1,2,1,m,amrlev,mglev,mfi);
		left(3,6) = m_operator->C(0,1,0,2,m,amrlev,mglev,mfi);
		left(3,7) = m_operator->C(0,1,1,2,m,amrlev,mglev,mfi);
		left(3,8) = m_operator->C(0,1,2,2,m,amrlev,mglev,mfi);
		left(4,0) = m_operator->C(1,1,0,0,m,amrlev,mglev,mfi);
		left(4,1) = m_operator->C(1,1,1,0,m,amrlev,mglev,mfi);
		left(4,2) = m_operator->C(1,1,2,0,m,amrlev,mglev,mfi);
		left(4,3) = m_operator->C(1,1,0,1,m,amrlev,mglev,mfi);
		left(4,4) = m_operator->C(1,1,1,1,m,amrlev,mglev,mfi);
		left(4,5) = m_operator->C(1,1,2,1,m,amrlev,mglev,mfi);
		left(4,6) = m_operator->C(1,1,0,2,m,amrlev,mglev,mfi);
		left(4,7) = m_operator->C(1,1,1,2,m,amrlev,mglev,mfi);
		left(4,8) = m_operator->C(1,1,2,2,m,amrlev,mglev,mfi);
		left(5,0) = m_operator->C(2,1,0,0,m,amrlev,mglev,mfi);
		left(5,1) = m_operator->C(2,1,1,0,m,amrlev,mglev,mfi);
		left(5,2) = m_operator->C(2,1,2,0,m,amrlev,mglev,mfi);
		left(5,3) = m_operator->C(2,1,0,1,m,amrlev,mglev,mfi);
		left(5,4) = m_operator->C(2,1,1,1,m,amrlev,mglev,mfi);
		left(5,5) = m_operator->C(2,1,2,1,m,amrlev,mglev,mfi);
		left(5,6) = m_operator->C(2,1,0,2,m,amrlev,mglev,mfi);
		left(5,7) = m_operator->C(2,1,1,2,m,amrlev,mglev,mfi);
		left(5,8) = m_operator->C(2,1,2,2,m,amrlev,mglev,mfi);
		left(6,0) = m_operator->C(0,2,0,0,m,amrlev,mglev,mfi);
		left(6,1) = m_operator->C(0,2,1,0,m,amrlev,mglev,mfi);
		left(6,2) = m_operator->C(0,2,2,0,m,amrlev,mglev,mfi);
		left(6,3) = m_operator->C(0,2,0,1,m,amrlev,mglev,mfi);
		left(6,4) = m_operator->C(0,2,1,1,m,amrlev,mglev,mfi);
		left(6,5) = m_operator->C(0,2,2,1,m,amrlev,mglev,mfi);
		left(6,6) = m_operator->C(0,2,0,2,m,amrlev,mglev,mfi);
		left(6,7) = m_operator->C(0,2,1,2,m,amrlev,mglev,mfi);
		left(6,8) = m_operator->C(0,2,2,2,m,amrlev,mglev,mfi);
		left(7,0) = m_operator->C(1,2,0,0,m,amrlev,mglev,mfi);
		left(7,1) = m_operator->C(1,2,1,0,m,amrlev,mglev,mfi);
		left(7,2) = m_operator->C(1,2,2,0,m,amrlev,mglev,mfi);
		left(7,3) = m_operator->C(1,2,0,1,m,amrlev,mglev,mfi);
		left(7,4) = m_operator->C(1,2,1,1,m,amrlev,mglev,mfi);
		left(7,5) = m_operator->C(1,2,2,1,m,amrlev,mglev,mfi);
		left(7,6) = m_operator->C(1,2,0,2,m,amrlev,mglev,mfi);
		left(7,7) = m_operator->C(1,2,1,2,m,amrlev,mglev,mfi);
		left(7,8) = m_operator->C(1,2,2,2,m,amrlev,mglev,mfi);
		left(7,0) = m_operator->C(2,2,0,0,m,amrlev,mglev,mfi);
		left(8,1) = m_operator->C(2,2,1,0,m,amrlev,mglev,mfi);
		left(8,2) = m_operator->C(2,2,2,0,m,amrlev,mglev,mfi);
		left(8,3) = m_operator->C(2,2,0,1,m,amrlev,mglev,mfi);
		left(8,4) = m_operator->C(2,2,1,1,m,amrlev,mglev,mfi);
		left(8,5) = m_operator->C(2,2,2,1,m,amrlev,mglev,mfi);
		left(8,6) = m_operator->C(2,2,0,2,m,amrlev,mglev,mfi);
		left(8,7) = m_operator->C(2,2,1,2,m,amrlev,mglev,mfi);
		left(8,8) = m_operator->C(2,2,2,2,m,amrlev,mglev,mfi);

		right(0) = mul1*traction[0](0);
		right(1) = mul1*traction[0](1);
		right(2) = mul1*traction[0](2);
		right(3) = mul2*traction[1](0);
		right(4) = mul2*traction[1](1);
		right(5) = mul2*traction[1](2);
		right(6) = mul3*traction[2](0);
		right(7) = mul3*traction[2](1);
		right(8) = mul3*traction[2](2);

		//sol = left.ldlt().solve(right); // we can change this solver as needed
		sol = left.colPivHouseholderQr().solve(right);

		stencil[points[0]](0) = stencil[0](0) + mul1*DX[0]*sol(0);
		stencil[points[0]](1) = stencil[0](1) + mul1*DX[0]*sol(1);
		stencil[points[0]](2) = stencil[0](2) + mul1*DX[0]*sol(2);

		stencil[points[1]](0) = stencil[0](0) + mul2*DX[1]*sol(3);
		stencil[points[1]](1) = stencil[0](1) + mul2*DX[1]*sol(4);
		stencil[points[1]](2) = stencil[0](2) + mul2*DX[1]*sol(5);

		stencil[points[2]](0) = stencil[0](0) + mul3*DX[2]*sol(6);
		stencil[points[2]](1) = stencil[0](1) + mul3*DX[2]*sol(7);
		stencil[points[2]](2) = stencil[0](2) + mul3*DX[2]*sol(8);
	}
#endif
}
}
