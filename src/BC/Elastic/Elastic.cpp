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

	/* The following steps are for debugging purposes.
		They can be disabled by setting debug = false*/
	bool debug = true;
	amrex::IntVect test_point(-1,3,3);
	if(debug)
	{
		std::cout << "Box loVect = (" << box.loVect()[0] << "," << box.loVect()[1] << "," << box.loVect()[2] << ")" << std::endl;
		std::cout << "Box hiVect = (" << box.hiVect()[0] << "," << box.hiVect()[1] << "," << box.hiVect()[2] << ")" << std::endl;
		std::cout << "Domain loVect = (" << domain.loVect()[0] << "," << domain.loVect()[1] << "," << domain.loVect()[2] << ")" << std::endl;
		std::cout << "Domain hiVect = (" << domain.hiVect()[0] << "," << domain.hiVect()[1] << "," << domain.hiVect()[2] << ")" << std::endl;
		std::cout << "Amrlev = " << m_amrlev << ". Mglev = " << m_mglev << std::endl;
		std::cout << __LINE__ << ": Value at test point = (" << mf_box(test_point,0) << "," << mf_box(test_point,1) << "," << mf_box(test_point,2) << ")" << std::endl;
	}

	//mf.FillBoundary(m_geom.periodicity());

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	
	/* Dirichlet boundaries are first */
	AMREX_D_TERM(	for (int i = box.loVect()[0] - ngrow; i<=box.hiVect()[0] + ngrow; i++),
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


	AMREX_D_TERM(	for (int i = box.loVect()[0] - ngrow; i<=box.hiVect()[0] + ngrow; i++),
			for (int j = box.loVect()[1] - ngrow; j<=box.hiVect()[1] + ngrow; j++),
			for (int k = box.loVect()[2] - ngrow; k<=box.hiVect()[2] + ngrow; k++))
	{
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		amrex::Vector<Set::Vector> stencil;
		amrex::Vector<Set::Vector> traction;
		amrex::Vector<int> points;

		if (i == domain.loVect()[0]-1 && (face == Orientation::xlo || face == Orientation::All) && bc_lo_str[0] == "traction") // Left boundary
		{
			if (j == domain.loVect()[1] - 1) continue;
			if (j == domain.hiVect()[1] + 1) continue;
			if (k == domain.loVect()[2] - 1) continue;
			if (k == domain.hiVect()[2] + 1) continue;
			
			stencil.clear();
			traction.clear();
			points.clear();			
			
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

			traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
			points.push_back(1);
#if AMREX+SAPCEDIM > 1
			if (j==domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
			}
			else if (j==domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
			}
#if AMREX_SPACEDIM > 2
			if (k==domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k==domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
#endif
			StencilFill(stencil, traction, points, m+dx, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[1](n);
		}
		if (i == domain.hiVect()[0]+1 && (face == Orientation::xhi || face == Orientation::All) && bc_hi_str[0] == "traction") // Left boundary
		{
			if (j == domain.loVect()[1] - 1) continue;
			if (j == domain.hiVect()[1] + 1) continue;
			if (k == domain.loVect()[2] - 1) continue;
			if (k == domain.hiVect()[2] + 1) continue;
			stencil.clear();
			traction.clear();
			points.clear();	

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

			traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
			points.push_back(2);
#if AMREX+SAPCEDIM > 1
			if (j==domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
			}
			else if (j==domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
			}
#if AMREX_SPACEDIM > 2
			if (k==domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k==domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
#endif
			StencilFill(stencil, traction, points, m-dx, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[2](n);
		}
#if AMREX_SPACEDIM > 1
		if (j == domain.loVect()[1]-1 && (face == Orientation::ylo || face == Orientation::All) && bc_lo_str[1] == "traction") // Left boundary
		{
			if (i == domain.loVect()[0] - 1) continue;
			if (i == domain.hiVect()[0] + 1) continue;
			if (k == domain.loVect()[2] - 1) continue;
			if (k == domain.hiVect()[2] + 1) continue;
			stencil.clear();
			traction.clear();
			points.clear();	

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
			if (i == domain.loVect()[0] && bc_lo_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
				points.push_back(1);
			}
			else if (i == domain.hiVect()[0] && bc_hi_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));
				points.push_back(2);
			}

			traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
			points.push_back(3);

#if AMREX_SPACEDIM > 2
			if (k==domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k==domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
			StencilFill(stencil, traction, points, m+dy, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[3](n);

		}
		if (j == domain.hiVect()[1]+1 && (face == Orientation::yhi || face == Orientation::All) && bc_hi_str[1] == "traction") // Left boundary
		{
			if (i == domain.loVect()[0] - 1) continue;
			if (i == domain.hiVect()[0] + 1) continue;
			if (k == domain.loVect()[2] - 1) continue;
			if (k == domain.hiVect()[2] + 1) continue;
			stencil.clear();
			traction.clear();
			points.clear();	

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
			if (i == domain.loVect()[0] && bc_lo_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
				points.push_back(1);
			}
			else if (i == domain.hiVect()[0] && bc_hi_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));
				points.push_back(2);
			}

			traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
			points.push_back(4);

#if AMREX_SPACEDIM > 2
			if (k==domain.loVect()[2] && bc_lo_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
				points.push_back(5);
			}
			else if (k==domain.hiVect()[2] && bc_hi_str[2] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
				points.push_back(6);
			}
#endif
			StencilFill(stencil, traction, points, m-dy, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[4](n);
		}
#if AMREX_SPACEDIM > 2
		if (k == domain.loVect()[2]-1 && (face == Orientation::zlo || face == Orientation::All) && bc_lo_str[2] == "traction") // Left boundary
		{
			if (i == domain.loVect()[0] - 1) continue;
			if (i == domain.hiVect()[0] + 1) continue;
			if (j == domain.loVect()[1] - 1) continue;
			if (j == domain.hiVect()[1] + 1) continue;
			stencil.clear();
			traction.clear();
			points.clear();	

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

			if (i == domain.loVect()[0] && bc_lo_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
				points.push_back(1);
			}
			else if (i == domain.hiVect()[0] && bc_hi_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));
				points.push_back(2);
			}
			if (j==domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
			}
			else if (j==domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
			}
			traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[0],bc_lo_3[1],bc_lo_3[2])));
			points.push_back(5);

			StencilFill(stencil, traction, points, m+dz, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[5](n);
		}
		if (k == domain.hiVect()[2]+1 && (face == Orientation::zhi || face == Orientation::All) && bc_hi_str[2] == "traction") // Left boundary
		{
			if (i == domain.loVect()[0] - 1) continue;
			if (i == domain.hiVect()[0] + 1) continue;
			if (j == domain.loVect()[1] - 1) continue;
			if (j == domain.hiVect()[1] + 1) continue;
			stencil.clear();
			traction.clear();
			points.clear();	

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

			if (i == domain.loVect()[0] && bc_lo_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[0],bc_lo_1[1],bc_lo_1[2])));
				points.push_back(1);
			}
			else if (i == domain.hiVect()[0] && bc_hi_str[0] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[0],bc_hi_1[1],bc_hi_1[2])));
				points.push_back(2);
			}
			if (j==domain.loVect()[1] && bc_lo_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[0],bc_lo_2[1],bc_lo_2[2])));
				points.push_back(3);
			}
			else if (j==domain.hiVect()[1] && bc_hi_str[1] == "traction")
			{
				traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[0],bc_hi_2[1],bc_hi_2[2])));
				points.push_back(4);
			}
			traction.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[0],bc_hi_3[1],bc_hi_3[2])));
			points.push_back(6);

			StencilFill(stencil, traction, points, m-dz, m_amrlev, m_mglev, mfi, debug);
			for (int n = 0; n < AMREX_SPACEDIM; n++)
				mf_box(m,n) = stencil[6](n);
		}
#endif
#endif
	}

	/*Finally let's do averaging of corner cells - edge corners and triple corners*/
	AMREX_D_TERM(	for(int i = box.loVect()[0]-ngrow; i <= box.hiVect()[0]+ngrow; i++),
			for(int j = box.loVect()[1]-ngrow; j <= box.hiVect()[1]+ngrow; j++),
			for(int k = box.loVect()[2]-ngrow; k <= box.hiVect()[2]+ngrow; k++))
	{
		amrex::IntVect m(AMREX_D_DECL(i,j,k));
		int mul1 = 0, mul2 = 0, mul3 = 0;

		if (i == domain.loVect()[0]-1 && (face == Orientation::xlo || face == Orientation::All) && bc_lo_str[0] == "traction")
		{
			mul1 = 1;
#if AMREX_SPACEDIM > 1
			if (j == domain.loVect()[1] - 1) mul2 = 1;
			else if (j==domain.hiVect()[1] +1) mul2 = -1;
			else mul2 = 0;
#if AMREX_SPACEDIM > 2
			if (k == domain.loVect()[2] - 1) mul3 = 1;
			else if (k==domain.hiVect()[2] +1) mul3 = -1;
			else mul3 = 0;
#endif
#endif
			if (mul2 == 0 && mul3 == 0) mul1 = 0;
		}
		if (i == domain.hiVect()[0]+1 && (face == Orientation::xhi || face == Orientation::All) && bc_hi_str[0] == "traction")
		{
			mul1 = -1;
#if AMREX_SPACEDIM > 1
			if (j == domain.loVect()[1] - 1) mul2 = 1;
			else if (j==domain.hiVect()[1] +1) mul2 = -1;
			else mul2 = 0;
#if AMREX_SPACEDIM > 2
			if (k == domain.loVect()[2] - 1) mul3 = 1;
			else if (k==domain.hiVect()[2] +1) mul3 = -1;
			else mul3 = 0;
#endif
#endif
			if (mul2 == 0 && mul3 == 0) mul1 = 0;
		}

#if AMREX_SPACEDIM > 1
		if (j == domain.loVect()[1]-1 && (face == Orientation::ylo || face == Orientation::All) && bc_lo_str[1] == "traction")
		{
			mul2 = 1;
			if (i == domain.loVect()[0] -1) mul1 = 1;
			else if (i == domain.hiVect()[0]+1) mul1 = -1;
			else mul1 = 0;
#if AMREX_SPACEDIM > 2
			if (k == domain.loVect()[2]-1) mul3 = 1;
			else if (k == domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
			if (mul1 == 0 && mul3 == 0) mul2 = 0;
		}
		if (j == domain.hiVect()[1]+1 && (face == Orientation::yhi || face == Orientation::All) && bc_hi_str[1] == "traction")
		{
			mul2 = -1;
			if (i == domain.loVect()[0] -1) mul1 = 1;
			else if (i == domain.hiVect()[0]+1) mul1 = -1;
			else mul1 = 0;
#if AMREX_SPACEDIM > 2
			if (k == domain.loVect()[2]-1) mul3 = 1;
			else if (k == domain.hiVect()[2]+1) mul3 = -1;
			else mul3 = 0;
#endif
			if (mul2 == 0 && mul3 == 0) mul1 = 0;
		}
#endif

#if AMREX_SPACEDIM > 2
		if (k == domain.loVect()[2]-1 && (face == Orientation::zlo || face == Orientation::All) && bc_lo_str[2] == "traction")
		{
			mul3 = 1;
			if (i == domain.loVect()[0] -1) mul1 = 1;
			else if (i == domain.hiVect()[0]+1) mul1 = -1;
			else mul1 = 0;

			if (j == domain.loVect()[1] - 1) mul2 = 1;
			else if (j==domain.hiVect()[1] +1) mul2 = -1;
			else mul2 = 0;

			if (mul1 == 0 && mul2 == 0) mul3 = 0;

		}
		if (k == domain.hiVect()[2]+1 && (face == Orientation::zhi || face == Orientation::All) && bc_hi_str[2] == "traction")
		{
			mul3 = -1;
			if (i == domain.loVect()[0] -1) mul1 = 1;
			else if (i == domain.hiVect()[0]+1) mul1 = -1;
			else mul1 = 0;

			if (j == domain.loVect()[1] - 1) mul2 = 1;
			else if (j==domain.hiVect()[1] +1) mul2 = -1;
			else mul2 = 0;

			if (mul1 == 0 && mul2 == 0) mul3 = 0;
		}
#endif
		if(AMREX_D_TERM(std::abs(mul1), + std::abs(mul2), + std::abs(mul3)) == 0) continue;

		for(int n = 0; n<AMREX_SPACEDIM; n++)
			mf_box(m,n) = (AMREX_D_TERM(mf_box(m+mul1*dx,n),+mf_box(m+mul2*dy,n),+mf_box(m+mul3*dz,n)))/(AMREX_D_TERM(1.0*std::abs(mul1),+1.0*std::abs(mul2),+1.0*std::abs(mul3)));
	}

	if(debug)
	{
		std::cout << __LINE__ << ": Value at test point = (" << mf_box(test_point,0) << "," << mf_box(test_point,1) << "," << mf_box(test_point,2) << ")" << std::endl;
	}
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
	if(debug)
	{
		for (int p = 0; p<stencil.size(); p++)
		{
			if(std::isnan(stencil[p](0)) || std::isnan(stencil[p](1)) || std::isnan(stencil[p](2)))
			{
				std::cout << "Error: Nans found in stencil after computation. p = " << p << std::endl;
			}
		}
	}
}
}
