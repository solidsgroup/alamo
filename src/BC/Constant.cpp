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
Constant::FillBoundary (amrex::BaseFab<Set::Scalar> &a_in,
            const amrex::Box &a_box,
            int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real time,
            Orientation face, const amrex::Mask * /*mask*/)
{
    const auto DX = m_geom.CellSizeArray();

    Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp));

    amrex::Box box = a_box;
    box.grow(ngrow);
    amrex::Box domain = m_geom.Domain();
    const amrex::IndexType type = a_in.box().ixType();
    domain.convert(type);
    const bool nodal = type == amrex::IndexType::TheNodeType();
    const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

    amrex::Array4<amrex::Real> const& in = a_in.array();

    for (int n = 0; n < a_in.nComp(); n++)
    {
        const auto bc_type_xlo = m_bc_type[Face::XLO][n];
        const auto bc_type_xhi = m_bc_type[Face::XHI][n];
        const auto bc_type_ylo = m_bc_type[Face::YLO][n];
        const auto bc_type_yhi = m_bc_type[Face::YHI][n];
        const Set::Scalar bc_val_xlo =
            m_bc_val[Face::XLO].empty() ? 0.0 : m_bc_val[Face::XLO][n](time);
        const Set::Scalar bc_val_xhi =
            m_bc_val[Face::XHI].empty() ? 0.0 : m_bc_val[Face::XHI][n](time);
        const Set::Scalar bc_val_ylo =
            m_bc_val[Face::YLO].empty() ? 0.0 : m_bc_val[Face::YLO][n](time);
        const Set::Scalar bc_val_yhi =
            m_bc_val[Face::YHI].empty() ? 0.0 : m_bc_val[Face::YHI][n](time);
#if AMREX_SPACEDIM > 2
        const auto bc_type_zlo = m_bc_type[Face::ZLO][n];
        const auto bc_type_zhi = m_bc_type[Face::ZHI][n];
        const Set::Scalar bc_val_zlo =
            m_bc_val[Face::ZLO].empty() ? 0.0 : m_bc_val[Face::ZLO][n](time);
        const Set::Scalar bc_val_zhi =
            m_bc_val[Face::ZHI].empty() ? 0.0 : m_bc_val[Face::ZHI][n](time);
#endif
        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            amrex::IntVect glevel;
            AMREX_D_TERM(   glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
                            glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
                            glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );

            if ((glevel[0]<0 || (nodal && i == lo.x)) && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xlo))
                {
                    if (nodal) in(i,j,k,n) = bc_val_xlo;
                    else if (glevel[0]<0) in(i,j,k,n) = 2.0 * bc_val_xlo - in(2 * lo.x - i - 1,j,k,n);
                }
                else if(glevel[0]<0 && BCUtil::IsNeumann(bc_type_xlo))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_val_xlo*DX[0];
                else if(glevel[0]<0 && BCUtil::IsReflectEven(bc_type_xlo))
                    in(i,j,k,n) = in(1-glevel[0],j,k,n);
                else if(glevel[0]<0 && BCUtil::IsReflectOdd(bc_type_xlo))
                    in(i,j,k,n) = -in(1-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xlo) || (nodal && i == lo.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[0]>0 || (nodal && i == hi.x)) && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xhi))
                {
                    if (nodal) in(i,j,k,n) = bc_val_xhi;
                    else if (glevel[0]>0) in(i,j,k,n) = 2.0 * bc_val_xhi - in(2 * hi.x - i + 1,j,k,n);
                }
                else if(glevel[0]>0 && BCUtil::IsNeumann(bc_type_xhi))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_val_xhi*DX[0];
                else if(glevel[0]>0 && BCUtil::IsReflectEven(bc_type_xhi))
                    in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
                else if(glevel[0]>0 && BCUtil::IsReflectOdd(bc_type_xhi))
                    in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xhi) || (nodal && i == hi.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }

            else if ((glevel[1]<0 || (nodal && j == lo.y)) && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
            {
                if (BCUtil::IsDirichlet(bc_type_ylo))
                {
                    if (nodal) in(i,j,k,n) = bc_val_ylo;
                    else if (glevel[1]<0) in(i,j,k,n) = 2.0 * bc_val_ylo - in(i,2 * lo.y - j - 1,k,n);
                }
                else if (glevel[1]<0 && BCUtil::IsNeumann(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_val_ylo*DX[1];
                else if (glevel[1]<0 && BCUtil::IsReflectEven(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n);
                else if (glevel[1]<0 && BCUtil::IsReflectOdd(bc_type_ylo))
                    in(i,j,k,n) = -in(i,j-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_ylo) || (nodal && j == lo.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[1]>0 || (nodal && j == hi.y)) && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
            {
                if (BCUtil::IsDirichlet(bc_type_yhi))
                {
                    if (nodal) in(i,j,k,n) = bc_val_yhi;
                    else if (glevel[1]>0) in(i,j,k,n) = 2.0 * bc_val_yhi - in(i,2 * hi.y - j + 1,k,n);
                }
                else if (glevel[1]>0 && BCUtil::IsNeumann(bc_type_yhi))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_val_yhi*DX[1];
                else if (glevel[1]>0 && BCUtil::IsReflectEven(bc_type_yhi))
                    in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
                else if (glevel[1]>0 && BCUtil::IsReflectOdd(bc_type_yhi))
                    in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_yhi) || (nodal && j == hi.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }

#if AMREX_SPACEDIM>2
            else if ((glevel[2]<0 || (nodal && k == lo.z)) && (face == Orientation::zlo || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zlo))
                {
                    if (nodal) in(i,j,k,n) = bc_val_zlo;
                    else if (glevel[2]<0) in(i,j,k,n) = 2.0 * bc_val_zlo - in(i,j,2 * lo.z - k - 1,n);
                }
                else if (glevel[2]<0 && BCUtil::IsNeumann(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_val_zlo*DX[2];
                else if (glevel[2]<0 && BCUtil::IsReflectEven(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,1-glevel[2],n);
                else if (glevel[2]<0 && BCUtil::IsReflectOdd(bc_type_zlo))
                    in(i,j,k,n) = -in(i,j,1-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zlo) || (nodal && k == lo.z)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[2]>0 || (nodal && k == hi.z)) && (face == Orientation::zhi || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zhi))
                {
                    if (nodal) in(i,j,k,n) = bc_val_zhi;
                    else if (glevel[2]>0) in(i,j,k,n) = 2.0 * bc_val_zhi - in(i,j,2 * hi.z - k + 1,n);
                }
                else if(glevel[2]>0 && BCUtil::IsNeumann(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_val_zhi*DX[2];
                else if(glevel[2]>0 && BCUtil::IsReflectEven(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
                else if(glevel[2]>0 && BCUtil::IsReflectOdd(bc_type_zhi))
                    in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zhi) || (nodal && k == hi.z)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
#endif
        });

        // Average the value of corner cells based on the neighboring ghost cells.
        // This fixes NAN issues that can arise from neumann conditions calculated in ghost cells
        // Fixes this issue in 2D only for now.
        //
        // TODO: more general fix for 3D cell-based fields
        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (i < lo.x && j < lo.y) in(i,j,k,n) = 0.5*( in(i+1,j,k,n) + in(i,j+1,k,n) );
            if (i < lo.x && j > hi.y) in(i,j,k,n) = 0.5*( in(i+1,j,k,n) + in(i,j-1,k,n) );
            if (i > hi.x && j < lo.y) in(i,j,k,n) = 0.5*( in(i-1,j,k,n) + in(i,j+1,k,n) );
            if (i > hi.x && j > hi.y) in(i,j,k,n) = 0.5*( in(i-1,j,k,n) + in(i,j-1,k,n) );
        });

    }
}

amrex::BCRec
Constant::GetBCRec(int component)
{
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][component],m_bc_type[Face::YLO][component],m_bc_type[Face::ZLO][component])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][component],m_bc_type[Face::YHI][component],m_bc_type[Face::ZHI][component])};

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
