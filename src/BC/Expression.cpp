#include "Expression.H"

namespace BC
{

void
Expression::FillBoundary (amrex::BaseFab<Set::Scalar> &a_in,
                        const amrex::Box &a_box,
                        int ngrow, int /*dcomp*/, int /*ncomp*/, Set::Scalar time,
                        Orientation face, const amrex::Mask * /*mask*/)
{
    using Base = BC<Set::Scalar>;

    const amrex::Real* DX = Base::m_geom.CellSize();

    Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp)," - ", m_prefix);

    amrex::Box box = a_box;
    box.grow(ngrow);
    const amrex::Dim3 lo= amrex::lbound(Base::m_geom.Domain()), hi = amrex::ubound(Base::m_geom.Domain());

    amrex::Array4<amrex::Real> const& in = a_in.array();

    amrex::IndexType type = amrex::IndexType::TheCellType();

    for (int n = 0; n < a_in.nComp(); n++)
    amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Set::Vector pos = Set::Position(i, j, k, Base::m_geom, type);
        Set::Scalar x = 0.0, y=0.0, z=0.0, t=time;
        x = pos(0);
        #if AMREX_SPACEDIM > 1
        y = pos(1);
        #if AMREX_SPACEDIM > 2
        z = pos(2);
        #endif
        #endif


        amrex::IntVect glevel;
        AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
                    glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
                    glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );
        
        if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n]))
                in(i,j,k,n) = m_bc_func[Face::XLO][n](x,y,z,t);
            else if(BCUtil::IsNeumann(m_bc_type[Face::XLO][n]))
                in(i,j,k,n) = in(i-glevel[0],j,k,n) - (m_bc_func[Face::XLO].size() > 0 ? m_bc_func[Face::XLO][n](x,y,z,t)*DX[0] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::XLO][n]))
                in(i,j,k,n) = in(1-glevel[0],j,k,n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::XLO][n]))
                in(i,j,k,n) = -in(1-glevel[0],j,k,n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::XLO][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }
        else if (glevel[0]>0 && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n]))
                in(i,j,k,n) = m_bc_func[Face::XHI][n](x,y,z,t);
            else if(BCUtil::IsNeumann(m_bc_type[Face::XHI][n]))
                in(i,j,k,n) = in(i-glevel[0],j,k,n) - (m_bc_func[Face::XHI].size() > 0 ? m_bc_func[Face::XHI][n](x,y,z,t)*DX[0] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::XHI][n]))
                in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::XHI][n]))
                in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::XHI][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }
        
        else if (glevel[1]<0 && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n]))
                in(i,j,k,n) = m_bc_func[Face::YLO][n](x,y,z,t);
            else if (BCUtil::IsNeumann(m_bc_type[Face::YLO][n]))
                in(i,j,k,n) = in(i,j-glevel[1],k,n) - (m_bc_func[Face::YLO].size() > 0 ? m_bc_func[Face::YLO][n](x,y,z,t)*DX[1] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::YLO][n]))
                in(i,j,k,n) = in(i,j-glevel[1],k,n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::YLO][n]))
                in(i,j,k,n) = -in(i,j-glevel[1],k,n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::YLO][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }
        else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n]))
                in(i,j,k,n) = m_bc_func[Face::YHI][n](x,y,z,t);
            else if (BCUtil::IsNeumann(m_bc_type[Face::YHI][n]))
                in(i,j,k,n) = in(i,j-glevel[1],k,n) - (m_bc_func[Face::YHI].size() > 0 ? m_bc_func[Face::YHI][n](x,y,z,t)*DX[1] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::YHI][n]))
                in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::YHI][n]))
                in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::YHI][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }

#if AMREX_SPACEDIM>2
        else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n]))
                in(i,j,k,n) = m_bc_func[Face::ZLO][n](x,y,z,t);
            else if (BCUtil::IsNeumann(m_bc_type[Face::ZLO][n]))
                in(i,j,k,n) = in(i,j,k-glevel[2],n) - (m_bc_func[Face::ZLO].size() > 0 ? m_bc_func[Face::ZLO][n](x,y,z,t)*DX[2] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::ZLO][n]))
                in(i,j,k,n) = in(i,j,1-glevel[2],n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::ZLO][n]))
                in(i,j,k,n) = -in(i,j,1-glevel[2],n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::ZLO][n])) {}
            else Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }
        else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n]))
                in(i,j,k,n) = m_bc_func[Face::ZHI][n](x,y,z,t);
            else if(BCUtil::IsNeumann(m_bc_type[Face::ZHI][n]))
                in(i,j,k,n) = in(i,j,k-glevel[2],n) - (m_bc_func[Face::ZHI].size() > 0 ? m_bc_func[Face::ZHI][n](x,y,z,t)*DX[2] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::ZHI][n]))
                in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::ZHI][n]))
                in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::ZHI][n])) {}
            else Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
        }
#endif
    });
}


void
Expression::FillBoundary (amrex::BaseFab<Set::Vector> &a_in,
                        const amrex::Box &a_box,
                        int ngrow, int /*dcomp*/, int /*ncomp*/, Set::Scalar time,
                        Orientation face, const amrex::Mask * /*mask*/)
{
    using Base = BC<Set::Vector>;

    const amrex::Real* DX = Base::m_geom.CellSize();

    Util::Assert(INFO,TEST(a_in.nComp() == 1)," - ", m_prefix);
    Util::Assert(INFO,TEST(AMREX_SPACEDIM == (int)m_ncomp), " - ", m_prefix);

    amrex::Box box = a_box;
    box.grow(ngrow);
    const amrex::Dim3 lo= amrex::lbound(Base::m_geom.Domain()), hi = amrex::ubound(Base::m_geom.Domain());

    amrex::Array4<Set::Vector> const& in = a_in.array();

    amrex::IndexType type = amrex::IndexType::TheCellType();

    amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Set::Vector pos = Set::Position(i, j, k, Base::m_geom, type);
        Set::Scalar x = 0.0, y=0.0, z=0.0, t=time;
        x = pos(0);
        #if AMREX_SPACEDIM > 1
        y = pos(1);
        #if AMREX_SPACEDIM > 2
        z = pos(2);
        #endif
        #endif


        amrex::IntVect glevel;
        AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
                    glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
                    glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );
        
        for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
            if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n]))
                    in(i,j,k)(n) = m_bc_func[Face::XLO][n](x,y,z,t);
                else if(BCUtil::IsNeumann(m_bc_type[Face::XLO][n]))
                    in(i,j,k)(n) = in(i-glevel[0],j,k)(n) - (m_bc_func[Face::XLO].size() > 0 ? m_bc_func[Face::XLO][n](x,y,z,t)*DX[0] : 0);
                else if(BCUtil::IsReflectEven(m_bc_type[Face::XLO][n]))
                    in(i,j,k)(n) = in(1-glevel[0],j,k)(n);
                else if(BCUtil::IsReflectOdd(m_bc_type[Face::XLO][n]))
                    in(i,j,k)(n) = -in(1-glevel[0],j,k)(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::XLO][n])) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }
            else if (glevel[0]>0 && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n]))
                    in(i,j,k)(n) = m_bc_func[Face::XHI][n](x,y,z,t);
                else if(BCUtil::IsNeumann(m_bc_type[Face::XHI][n]))
                    in(i,j,k)(n) = in(i-glevel[0],j,k)(n) - (m_bc_func[Face::XHI].size() > 0 ? m_bc_func[Face::XHI][n](x,y,z,t)*DX[0] : 0);
                else if(BCUtil::IsReflectEven(m_bc_type[Face::XHI][n]))
                    in(i,j,k)(n) = in(hi.x-glevel[0],j,k)(n);
                else if(BCUtil::IsReflectOdd(m_bc_type[Face::XHI][n]))
                    in(i,j,k)(n) = -in(hi.x-glevel[0],j,k)(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::XHI][n])) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }
        
            else if (glevel[1]<0 && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n]))
                    in(i,j,k)(n) = m_bc_func[Face::YLO][n](x,y,z,t);
                else if (BCUtil::IsNeumann(m_bc_type[Face::YLO][n]))
                    in(i,j,k)(n) = in(i,j-glevel[1],k)(n) - (m_bc_func[Face::YLO].size() > 0 ? m_bc_func[Face::YLO][n](x,y,z,t)*DX[1] : 0);
                else if (BCUtil::IsReflectEven(m_bc_type[Face::YLO][n]))
                    in(i,j,k)(n) = in(i,j-glevel[1],k)(n);
                else if (BCUtil::IsReflectOdd(m_bc_type[Face::YLO][n]))
                    in(i,j,k)(n) = -in(i,j-glevel[1],k)(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::YLO][n])) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }
            else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n]))
                    in(i,j,k)(n) = m_bc_func[Face::YHI][n](x,y,z,t);
                else if (BCUtil::IsNeumann(m_bc_type[Face::YHI][n]))
                    in(i,j,k)(n) = in(i,j-glevel[1],k)(n) - (m_bc_func[Face::YHI].size() > 0 ? m_bc_func[Face::YHI][n](x,y,z,t)*DX[1] : 0);
                else if (BCUtil::IsReflectEven(m_bc_type[Face::YHI][n]))
                    in(i,j,k)(n) = in(i,hi.y-glevel[1],k)(n);
                else if (BCUtil::IsReflectOdd(m_bc_type[Face::YHI][n]))
                    in(i,j,k)(n) = -in(i,hi.y-glevel[1],k)(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::YHI][n])) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }

#if AMREX_SPACEDIM>2
            else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n]))
                    in(i,j,k)(n) = m_bc_func[Face::ZLO][n](x,y,z,t);
                else if (BCUtil::IsNeumann(m_bc_type[Face::ZLO][n]))
                    in(i,j,k)(n) = in(i,j,k-glevel[2])(n) - (m_bc_func[Face::ZLO].size() > 0 ? m_bc_func[Face::ZLO][n](x,y,z,t)*DX[2] : 0);
                else if (BCUtil::IsReflectEven(m_bc_type[Face::ZLO][n]))
                    in(i,j,k)(n) = in(i,j,1-glevel[2])(n);
                else if (BCUtil::IsReflectOdd(m_bc_type[Face::ZLO][n]))
                    in(i,j,k)(n) = -in(i,j,1-glevel[2])(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::ZLO][n])) {}
                else Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }
            else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n]))
                    in(i,j,k)(n) = m_bc_func[Face::ZHI][n](x,y,z,t);
                else if(BCUtil::IsNeumann(m_bc_type[Face::ZHI][n]))
                    in(i,j,k)(n) = in(i,j,k-glevel[2])(n) - (m_bc_func[Face::ZHI].size() > 0 ? m_bc_func[Face::ZHI][n](x,y,z,t)*DX[2] : 0);
                else if(BCUtil::IsReflectEven(m_bc_type[Face::ZHI][n]))
                    in(i,j,k)(n) = in(i,j,hi.z-glevel[2])(n);
                else if(BCUtil::IsReflectOdd(m_bc_type[Face::ZHI][n]))
                    in(i,j,k)(n) = -in(i,j,hi.z-glevel[2])(n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::ZHI][n])) {}
                else Util::Abort(INFO, "Incorrect boundary conditions", m_prefix);
            }
#endif
        }

    });
}



amrex::BCRec
Expression::GetBCRec() 
{
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][0],m_bc_type[Face::YLO][0],m_bc_type[Face::XLO][0])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][0],m_bc_type[Face::YHI][0],m_bc_type[Face::XHI][0])};

    return amrex::BCRec(bc_lo,bc_hi);
}

amrex::Array<int,AMREX_SPACEDIM>
Expression::IsPeriodic()
{
    return {AMREX_D_DECL(BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))};
}
amrex::Periodicity Expression::Periodicity () const
{
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(BC<Set::Scalar>::m_geom.Domain().length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                          BC<Set::Scalar>::m_geom.Domain().length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                          BC<Set::Scalar>::m_geom.Domain().length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}
amrex::Periodicity Expression::Periodicity (const amrex::Box& b) {
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                          b.length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                          b.length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));

}


}
