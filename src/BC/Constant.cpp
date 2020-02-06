#include "Constant.H"

namespace BC
{


Constant::Constant (int a_ncomp,
			amrex::Vector<std::string> bc_hi_str,
		    amrex::Vector<std::string> bc_lo_str,
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
				 		 amrex::Vector<amrex::Real> _bc_lo_2,
				 		 amrex::Vector<amrex::Real> _bc_lo_3),
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
				 		 amrex::Vector<amrex::Real> _bc_hi_2,
				 		 amrex::Vector<amrex::Real> _bc_hi_3))
	//: 
	//AMREX_D_DECL(bc_lo_1(_bc_lo_1),bc_lo_2(_bc_lo_2),bc_lo_3(_bc_lo_3)),
	//AMREX_D_DECL(bc_hi_1(_bc_hi_1),bc_hi_2(_bc_hi_2),bc_hi_3(_bc_hi_3))
{
	Util::Warning(INFO,"This method is going away. Please use pp.queryclass() instead.");

	m_ncomp = a_ncomp;

	m_bc_type[Face::XLO].resize(m_ncomp,BCUtil::ReadString(bc_lo_str[0]));
	m_bc_type[Face::XHI].resize(m_ncomp,BCUtil::ReadString(bc_hi_str[0]));
	m_bc_type[Face::YLO].resize(m_ncomp,BCUtil::ReadString(bc_lo_str[1]));
	m_bc_type[Face::YHI].resize(m_ncomp,BCUtil::ReadString(bc_hi_str[1]));
	#if AMREX_SPACEDIM == 3
	m_bc_type[Face::ZLO].resize(m_ncomp,BCUtil::ReadString(bc_lo_str[2]));
	m_bc_type[Face::ZHI].resize(m_ncomp,BCUtil::ReadString(bc_hi_str[2]));
	#endif


	m_bc_val[Face::XLO].resize(m_ncomp,NAN);
	m_bc_val[Face::XHI].resize(m_ncomp,NAN);
	m_bc_val[Face::YLO].resize(m_ncomp,NAN);
	m_bc_val[Face::YHI].resize(m_ncomp,NAN);
	#if AMREX_SPACEDIM == 3
	m_bc_val[Face::ZLO].resize(m_ncomp,NAN);
	m_bc_val[Face::ZHI].resize(m_ncomp,NAN);
	#endif

	for (unsigned int i=0;i<m_ncomp;i++)
	{
		if (_bc_lo_1.size() > 0) m_bc_val[Face::XLO][i] = _bc_lo_1[i];
		if (_bc_hi_1.size() > 0) m_bc_val[Face::XHI][i] = _bc_hi_1[i];
		if (_bc_lo_2.size() > 0) m_bc_val[Face::YLO][i] = _bc_lo_2[i];
		if (_bc_hi_2.size() > 0) m_bc_val[Face::YHI][i] = _bc_hi_2[i];
		#if AMREX_SPACEDIM == 3
		if (_bc_lo_3.size() > 0) m_bc_val[Face::ZLO][i] = _bc_lo_3[i];
		if (_bc_hi_3.size() > 0) m_bc_val[Face::ZHI][i] = _bc_hi_3[i];
		#endif
	}
}


//amrex::Mask& m
void
Constant::FillBoundary (amrex::FArrayBox &a_in,
			const amrex::Box &a_box,
			int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real time,
			Orientation face, const amrex::Mask * /*mask*/)
{
	const amrex::Real* DX = m_geom.CellSize();

	Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp));

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
			if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n]))
				in(i,j,k,n) = m_bc_val[Face::XLO][n](time);
			else if(BCUtil::IsNeumann(m_bc_type[Face::XLO][n]))
				in(i,j,k,n) = in(i+1,j,k,n) - (m_bc_val[Face::XLO].size() > 0 ? m_bc_val[Face::XLO][n](time)*DX[0] : 0);
			else if(BCUtil::IsReflectEven(m_bc_type[Face::XLO][n]))
				in(i,j,k,n) = in(1-glevel[0],j,k,n);
			else if(BCUtil::IsReflectOdd(m_bc_type[Face::XLO][n]))
				in(i,j,k,n) = -in(1-glevel[0],j,k,n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::XLO][n])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[0]>0)// && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
		{
			if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n]))
				in(i,j,k,n) = m_bc_val[Face::XHI][n](time);
			else if(BCUtil::IsNeumann(m_bc_type[Face::XHI][n]))
				in(i,j,k,n) = in(i-1,j,k,n) - (m_bc_val[Face::XHI].size() > 0 ? m_bc_val[Face::XHI][n](time)*DX[0] : 0);
			else if(BCUtil::IsReflectEven(m_bc_type[Face::XHI][n]))
				in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
			else if(BCUtil::IsReflectOdd(m_bc_type[Face::XHI][n]))
				in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::XHI][n])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		
		else if (glevel[1]<0)// && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
		{
			if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n]))
				in(i,j,k,n) = m_bc_val[Face::YLO][n](time);
			else if (BCUtil::IsNeumann(m_bc_type[Face::YLO][n]))
				in(i,j,k,n) = in(i,j+1,k,n) - (m_bc_val[Face::YLO].size() > 0 ? m_bc_val[Face::YLO][n](time)*DX[1] : 0);
			else if (BCUtil::IsReflectEven(m_bc_type[Face::YLO][n]))
				in(i,j,k,n) = in(i,j-glevel[1],k,n);
			else if (BCUtil::IsReflectOdd(m_bc_type[Face::YLO][n]))
				in(i,j,k,n) = -in(i,j-glevel[1],k,n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::YLO][n])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
		{
			if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n]))
				in(i,j,k,n) = m_bc_val[Face::YHI][n](time);
			else if (BCUtil::IsNeumann(m_bc_type[Face::YHI][n]))
				in(i,j,k,n) = in(i,j-1,k,n) - (m_bc_val[Face::YHI].size() > 0 ? m_bc_val[Face::YHI][n](time)*DX[1] : 0);
			else if (BCUtil::IsReflectEven(m_bc_type[Face::YHI][n]))
				in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
			else if (BCUtil::IsReflectOdd(m_bc_type[Face::YHI][n]))
				in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::YHI][n])) {}
			else
				Util::Abort(INFO, "Incorrect boundary conditions");
		}

#if AMREX_SPACEDIM>2
		else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
		{
			if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n]))
				in(i,j,k,n) = m_bc_val[Face::ZLO][n](time);
			else if (BCUtil::IsNeumann(m_bc_type[Face::ZLO][n]))
				in(i,j,k,n) = in(i,j,k+1,n) - (m_bc_val[Face::ZLO].size() > 0 ? m_bc_val[Face::ZLO][n](time)*DX[2] : 0);
			else if (BCUtil::IsReflectEven(m_bc_type[Face::ZLO][n]))
				in(i,j,k,n) = in(i,j,1-glevel[2],n);
			else if (BCUtil::IsReflectOdd(m_bc_type[Face::ZLO][n]))
				in(i,j,k,n) = -in(i,j,1-glevel[2],n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::ZLO][n])) {}
			else Util::Abort(INFO, "Incorrect boundary conditions");
		}
		else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
		{
			if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n]))
				in(i,j,k,n) = m_bc_val[Face::ZHI][n](time);
			else if(BCUtil::IsNeumann(m_bc_type[Face::ZHI][n]))
				in(i,j,k,n) = in(i,j,k-1,n) - (m_bc_val[Face::ZHI].size() > 0 ? m_bc_val[Face::ZHI][n](time)*DX[2] : 0);
			else if(BCUtil::IsReflectEven(m_bc_type[Face::ZHI][n]))
				in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
			else if(BCUtil::IsReflectOdd(m_bc_type[Face::ZHI][n]))
				in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
			else if(BCUtil::IsPeriodic(m_bc_type[Face::ZHI][n])) {}
			else Util::Abort(INFO, "Incorrect boundary conditions");
		}
#endif


	});
}

amrex::BCRec
Constant::GetBCRec() 
{
	int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][0],m_bc_type[Face::YLO][0],m_bc_type[Face::XLO][0])};
	int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][0],m_bc_type[Face::YHI][0],m_bc_type[Face::XHI][0])};

	return amrex::BCRec(bc_lo,bc_hi);
}

amrex::Array<int,AMREX_SPACEDIM>
Constant::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
			     BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
			     BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))};
}
amrex::Periodicity Constant::Periodicity () const
{
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
							      			 			  m_geom.Domain().length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
							      			 			  m_geom.Domain().length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}
amrex::Periodicity Constant::Periodicity (const amrex::Box& b) {
	return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
							                              b.length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
							                              b.length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));

}


}


