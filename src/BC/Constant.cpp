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
Constant::FillBoundary (amrex::FArrayBox &a_in,
			const amrex::Box &a_box,
			int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real /*time*/,
			Orientation face, const amrex::Mask * /*mask*/)
{
	const amrex::Real* DX = m_geom.CellSize();

	amrex::Box box = a_box;
	box.grow(ngrow);
	const amrex::Dim3 lo= amrex::lbound(m_geom.Domain()), hi = amrex::ubound(m_geom.Domain());

	amrex::Array4<amrex::Real> const& in = a_in.array();

	for (int n = 0; n < a_in.nComp(); n++)
	amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
		amrex::IntVect glevel;
		AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
					 glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
					 glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );
		
		if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
		{
			if (BCUtil::IsDirichlet(bc_lo[0]))
				in(i,j,k,n) = bc_lo_1[n];
			else if(BCUtil::IsNeumann(bc_lo[0]))
				in(i,j,k,n) = in(i+1,j,k,n) - (bc_lo_1.size() > 0 ? bc_lo_1[n]*DX[0] : 0);
			else if(BCUtil::IsReflectEven(bc_lo[0]))
				in(i,j,k,n) = in(1-glevel[0],j,k,n);
			else if(BCUtil::IsReflectOdd(bc_lo[0]))
				in(i,j,k,n) = -in(1-glevel[0],j,k,n);
			else if(BCUtil::IsPeriodic(bc_lo[0])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[0]>0)// && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
		{
			if (BCUtil::IsDirichlet(bc_hi[0]))
				in(i,j,k,n) = bc_hi_1[n];
			else if(BCUtil::IsNeumann(bc_hi[0]))
				in(i,j,k,n) = in(i-1,j,k,n) - (bc_hi_1.size() > 0 ? bc_hi_1[n]*DX[0] : 0);
			else if(BCUtil::IsReflectEven(bc_hi[0]))
				in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
			else if(BCUtil::IsReflectOdd(bc_hi[0]))
				in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
			else if(BCUtil::IsPeriodic(bc_hi[0])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		
		else if (glevel[1]<0)// && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
		{
			if (BCUtil::IsDirichlet(bc_lo[1]))
				in(i,j,k,n) = bc_lo_2[n];
			else if (BCUtil::IsNeumann(bc_lo[1]))
				in(i,j,k,n) = in(i,j+1,k,n) - (bc_lo_2.size() > 0 ? bc_lo_2[n]*DX[1] : 0);
			else if (BCUtil::IsReflectEven(bc_lo[1]))
				in(i,j,k,n) = in(i,j-glevel[1],k,n);
			else if (BCUtil::IsReflectOdd(bc_lo[1]))
				in(i,j,k,n) = -in(i,j-glevel[1],k,n);
			else if(BCUtil::IsPeriodic(bc_lo[1])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
		{
			if (BCUtil::IsDirichlet(bc_hi[1]))
				in(i,j,k,n) = bc_hi_2[n];
			else if (BCUtil::IsNeumann(bc_hi[1]))
				in(i,j,k,n) = in(i,j-1,k,n) - (bc_hi_2.size() > 0 ? bc_hi_2[n]*DX[1] : 0);
			else if (BCUtil::IsReflectEven(bc_hi[1]))
				in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
			else if (BCUtil::IsReflectOdd(bc_hi[1]))
				in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
			else if(BCUtil::IsPeriodic(bc_hi[1])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}

#if AMREX_SPACEDIM>2
		else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
		{
			if (BCUtil::IsDirichlet(bc_lo[2]))
				in(i,j,k,n) = bc_lo_3[n];
			else if (BCUtil::IsNeumann(bc_lo[2]))
				in(i,j,k,n) = in(i,j,k+1,n) - (bc_lo_3.size() > 0 ? bc_lo_3[n]*DX[2] : 0);
			else if (BCUtil::IsReflectEven(bc_lo[2]))
				in(i,j,k,n) = in(i,j,1-glevel[2],n);
			else if (BCUtil::IsReflectOdd(bc_lo[2]))
				in(i,j,k,n) = -in(i,j,1-glevel[2],n);
			else if(BCUtil::IsPeriodic(bc_lo[2])) {}
			else Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
		{
			if (BCUtil::IsDirichlet(bc_hi[2]))
				in(i,j,k,n) = bc_hi_3[n];
			else if(BCUtil::IsNeumann(bc_hi[2]))
				in(i,j,k,n) = in(i,j,k-1,n) - (bc_hi_3.size() > 0 ? bc_hi_3[n]*DX[2] : 0);
			else if(BCUtil::IsReflectEven(bc_hi[2]))
				in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
			else if(BCUtil::IsReflectOdd(bc_hi[2]))
				in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
			else if(BCUtil::IsPeriodic(bc_hi[2])) {}
			else Util::Abort(INFO, "Incorrect boundary conditions");
		}
#endif


	});


/*
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

			

		}
	
	}
	*/
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


