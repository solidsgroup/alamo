#include "BC/IntBC.H"
#include "IO/ParmParse.H" // Including IO for ParmParse

namespace BC {

/// \brief FillBoundary implementation for NothingInt
void NothingInt::FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                              const amrex::Box& box,
                              int ngrow, int dcomp, int ncomp, amrex::Real time,
                              Orientation face, const amrex::Mask* mask) {
    // No action is taken; this is a "do-nothing" boundary condition.
}

/// \brief Constructor for ConstantInt with default values
ConstantInt::ConstantInt(int a_ncomp) : m_ncomp(a_ncomp) {}

/// \brief Constructor for ConstantInt using ParmParse
ConstantInt::ConstantInt(int a_ncomp, IO::ParmParse& pp, std::string name) : m_ncomp(a_ncomp) {
    pp_queryclass(name, *this);
}

/// \brief Constructor for ConstantInt with boundary condition values
ConstantInt::ConstantInt(int ncomp, amrex::Vector<std::string> bc_hi_str,
                         amrex::Vector<std::string> bc_lo_str,
                         AMREX_D_DECL(amrex::Vector<int> _bc_lo_1,
                                      amrex::Vector<int> _bc_lo_2,
                                      amrex::Vector<int> _bc_lo_3),
                         AMREX_D_DECL(amrex::Vector<int> _bc_hi_1,
                                      amrex::Vector<int> _bc_hi_2,
                                      amrex::Vector<int> _bc_hi_3))
    : m_ncomp(ncomp) 
{
    m_bc_type[Face::XLO].resize(m_ncomp, BCUtil::ReadString(bc_lo_str[0]));
    m_bc_type[Face::XHI].resize(m_ncomp, BCUtil::ReadString(bc_hi_str[0]));
    m_bc_type[Face::YLO].resize(m_ncomp, BCUtil::ReadString(bc_lo_str[1]));
    m_bc_type[Face::YHI].resize(m_ncomp, BCUtil::ReadString(bc_hi_str[1]));

#if AMREX_SPACEDIM == 3
    m_bc_type[Face::ZLO].resize(m_ncomp, BCUtil::ReadString(bc_lo_str[2]));
    m_bc_type[Face::ZHI].resize(m_ncomp, BCUtil::ReadString(bc_hi_str[2]));
#endif

    m_bc_val[Face::XLO].resize(m_ncomp, 0);
    m_bc_val[Face::XHI].resize(m_ncomp, 0);
    m_bc_val[Face::YLO].resize(m_ncomp, 0);
    m_bc_val[Face::YHI].resize(m_ncomp, 0);

#if AMREX_SPACEDIM == 3
    m_bc_val[Face::ZLO].resize(m_ncomp, 0);
    m_bc_val[Face::ZHI].resize(m_ncomp, 0);
#endif

    for (unsigned int i = 0; i < m_ncomp; i++) {
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

/// \brief FillBoundary implementation for ConstantInt
void ConstantInt::FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                               const amrex::Box& box,
                               int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real time,
                               Orientation face, const amrex::Mask* /*mask*/) {
    const amrex::Real* DX = m_geom.CellSize();

    Util::Assert(INFO, TEST(in.nComp() == (int)m_ncomp));

    amrex::Box grown_box = box;
    grown_box.grow(ngrow);
    const amrex::Dim3 lo = amrex::lbound(m_geom.Domain()), hi = amrex::ubound(m_geom.Domain());

    amrex::Array4<int> const& array_in = in.array();

    for (int n = 0; n < in.nComp(); n++) {
        amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::IntVect glevel;
            AMREX_D_TERM(glevel[0] = std::max(std::min(0, i - lo.x), i - hi.x);,
                         glevel[1] = std::max(std::min(0, j - lo.y), j - hi.y);,
                         glevel[2] = std::max(std::min(0, k - lo.z), k - hi.z););

            if (glevel[0] < 0 && (face == Orientation::xlo || face == Orientation::All)) {  // Left boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::XLO][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::XLO][n])) {
                    array_in(i, j, k, n) = array_in(i - glevel[0], j, k, n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            } else if (glevel[0] > 0 && (face == Orientation::xhi || face == Orientation::All)) {  // Right boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::XHI][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::XHI][n])) {
                    array_in(i, j, k, n) = array_in(i - glevel[0], j, k, n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            } else if (glevel[1] < 0 && (face == Orientation::ylo || face == Orientation::All)) {  // Bottom boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::YLO][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::YLO][n])) {
                    array_in(i, j, k, n) = array_in(i, j - glevel[1], k, n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            } else if (glevel[1] > 0 && (face == Orientation::yhi || face == Orientation::All)) {  // Top boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::YHI][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::YHI][n])) {
                    array_in(i, j, k, n) = array_in(i, j - glevel[1], k, n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            }

#if AMREX_SPACEDIM > 2
            else if (glevel[2] < 0 && (face == Orientation::zlo || face == Orientation::All)) {  // Front boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::ZLO][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::ZLO][n])) {
                    array_in(i, j, k, n) = array_in(i, j, k - glevel[2], n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            } else if (glevel[2] > 0 && (face == Orientation::zhi || face == Orientation::All)) {  // Back boundary
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n])) {
                    array_in(i, j, k, n) = m_bc_val[Face::ZHI][n];
                } else if (BCUtil::IsNeumann(m_bc_type[Face::ZHI][n])) {
                    array_in(i, j, k, n) = array_in(i, j, k - glevel[2], n);
                } else {
                    Util::Abort(INFO, "Incorrect boundary conditions");
                }
            }
#endif
        });
    }
}

/// \brief GetBCRec for ConstantInt
amrex::BCRec ConstantInt::GetBCRec() {
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][0], m_bc_type[Face::YLO][0], m_bc_type[Face::ZLO][0])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][0], m_bc_type[Face::YHI][0], m_bc_type[Face::ZHI][0])};
    return amrex::BCRec(bc_lo, bc_hi);
}

/// \brief IsPeriodic for ConstantInt
amrex::Array<int, AMREX_SPACEDIM> ConstantInt::IsPeriodic() {
    return {AMREX_D_DECL(BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                         BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                         BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))};
}

/// \brief Periodicity for ConstantInt
amrex::Periodicity ConstantInt::Periodicity() const {
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                          m_geom.Domain().length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                          m_geom.Domain().length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}

/// \brief Periodicity with Box input for ConstantInt
amrex::Periodicity ConstantInt::Periodicity(const amrex::Box& b) {
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                         b.length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                         b.length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}

} // namespace BC

