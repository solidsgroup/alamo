#include "Constant.H"

namespace BC
{


Constant::Constant (amrex::Vector<std::string> bc_hi_str,
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
			Util::Abort(INFO, "Invalid BCs cannot be periodic on one side and not the other");
	}
}


//amrex::Mask& m
void
Constant::FillBoundary (amrex::FArrayBox &in,
			const amrex::Box &box,
			int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real /*time*/,
			Orientation face, const amrex::Mask *mask)
{
	// if (mask != nullptr)
	// {
	// 	std::cout << in.loVect()[0] << " " << mask->loVect()[0] << std::endl;
	// 	std::cout << in.hiVect()[0] << " " << mask->hiVect()[0] << std::endl;
	// }

	//	nullptr;
	amrex::Box domain(m_geom.Domain());
	const amrex::Real* DX = m_geom.CellSize();

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	AMREX_D_TERM(for (int i = box.loVect()[0] - ngrow; i<=box.hiVect()[0] + ngrow; i++),
		     for (int j = box.loVect()[1] - ngrow; j<=box.hiVect()[1] + ngrow; j++),
		     for (int k = box.loVect()[2] - ngrow; k<=box.hiVect()[2] + ngrow; k++))
	{
		amrex::IntVect m(AMREX_D_DECL(i,j,k));

		if (mask != nullptr) 
		 	if (AMREX_D_TERM(mask->loVect()[0] <= i && i <= mask->hiVect()[0],
		 			 && mask->loVect()[1] <= j && j <= mask->hiVect()[1],
		 			 && mask->loVect()[2] <= k && k <= mask->hiVect()[2]))
				//std::cout << (*mask)(m) << std::endl;
			if ((*mask)(m) <= 0)
			{
				std::cout << "continuing"<< std::endl;
				continue;
			}

		for (int n = 0; n < in.nComp(); n++)
		{
			if (i == domain.loVect()[0]-1 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
			{
				if (BCUtil::IsDirichlet(bc_lo[0]))
					in(m,n) = bc_lo_1[n];
				else if(BCUtil::IsNeumann(bc_lo[0]))
					in(m,n) = in(m+dx,n) - bc_lo_1[n]*DX[0];
				else if(BCUtil::IsPeriodic(bc_lo[0])) continue;
				//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
			else if (i == domain.hiVect()[0]+1 && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
			{
				//if (mask != nullptr) std::cout << (*mask)(m) << std::endl;//if ((*mask)(m) <= 0) continue;
				if (BCUtil::IsDirichlet(bc_hi[0]))
					in(m,n) = bc_hi_1[n];
				else if(BCUtil::IsNeumann(bc_hi[0]))
					in(m,n) = in(m-dx,n) - bc_hi_1[n]*DX[0];
				else if(BCUtil::IsPeriodic(bc_hi[0])) continue;
					//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
			else if (j == domain.loVect()[1]-1 && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
			{
				if (BCUtil::IsDirichlet(bc_lo[1]))
					in(m,n) = bc_lo_2[n];
				else if (BCUtil::IsNeumann(bc_lo[1]))
					in(m,n) = in(m+dy,n) - bc_lo_2[n]*DX[1];
				else if(BCUtil::IsPeriodic(bc_lo[1])) continue;
				 	//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
			else if (j == domain.hiVect()[1]+1 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
			{
				if (BCUtil::IsDirichlet(bc_hi[1]))
					in(m,n) = bc_hi_2[n];
				else if (BCUtil::IsNeumann(bc_hi[1]))
					in(m,n) = in(m-dy,n) - bc_hi_2[n]*DX[1];
				else if(BCUtil::IsPeriodic(bc_hi[1])) continue;
					//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
#if AMREX_SPACEDIM>2
			else if (k == domain.loVect()[2]-1 && (face == Orientation::zlo || face == Orientation::All))
			{
				if (BCUtil::IsDirichlet(bc_lo[2]))
					in(m,n) = bc_lo_3[n];
				else if (BCUtil::IsNeumann(bc_lo[2]))
					in(m,n) = in(m+dz,n) - bc_lo_3[n]*DX[2];
				else if(BCUtil::IsPeriodic(bc_lo[2])) continue;
					//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
			else if (k == domain.hiVect()[2]+1 && (face == Orientation::zhi || face == Orientation::All))
			{
				if (BCUtil::IsDirichlet(bc_hi[2]))
					in(m,n) = bc_hi_3[n];
				else if(BCUtil::IsNeumann(bc_hi[2]))
					in(m,n) = in(m-dz,n) - bc_hi_3[n]*DX[2];
				else if(BCUtil::IsPeriodic(bc_hi[2])) continue;
				//in(m,n) = in(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
				else
					Util::Abort(INFO, "Incorrect boundary conditions");
			}
#endif

		}
	
	}
}

amrex::BCRec
Constant::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
Constant::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(bc_lo[0]),
			     BCUtil::IsPeriodic(bc_lo[1]),
			     BCUtil::IsPeriodic(bc_lo[2]))};
}
amrex::Periodicity Constant::Periodicity () const
{
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(bc_lo[0]),
							      m_geom.Domain().length(1) * BCUtil::IsPeriodic(bc_lo[1]),
							      m_geom.Domain().length(2) * BCUtil::IsPeriodic(bc_lo[2]))));
}
amrex::Periodicity Constant::Periodicity (const amrex::Box& b) {
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(bc_lo[0]),
							      b.length(1) * BCUtil::IsPeriodic(bc_lo[1]),
							      b.length(2) * BCUtil::IsPeriodic(bc_lo[2]))));

}


}


