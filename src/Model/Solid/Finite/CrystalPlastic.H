//
// A finite-deformation crystal-plasticity model with configurable slip systems
// and orthorhombic elasticity. FCC and orthorhombic ammonium-perchlorate slip
// sets are provided; custom slip systems can also be supplied in the input.
//
// Inherits from :ref:`Model::Solid::Finite::PseudoLinear::Orthotropic` to provide elastic response.
//
// Plastic flow modeled as:
//
// .. math::
//
//    \mathbf{F}^p\mathbf{F}^{p-1} =
//    \sum \dot{\gamma}_n \mathbf{a}_n\otimes\mathbf{N}_n
//
// where :math:`\mathbf{F}=\mathbf{F}^e\mathbf{F}^p` and :math:`\gamma_n` are
// slips on system :math:`n`.
// 
// Power law rate hardening is used in the integration of slips:
//
// .. math::
//
//    \dot{\gamma}_n = \dot{\gamma}_0 \Big(\frac{\tau}{\tau_{crss}}\Big)^m \operatorname{sign}(\tau)
//
//

#ifndef MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_


#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinear/Orthotropic.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
class CrystalPlastic : public PseudoLinear::Orthotropic
{
public:
    struct PlasticRate
    {
        Set::Matrix3d Fp = Set::Matrix3d::Zero();
        Eigen::Matrix<Set::Scalar,12,1> gamma =
            Eigen::Matrix<Set::Scalar,12,1>::Zero();
    };

    CrystalPlastic() = default;

    ~CrystalPlastic() = default;

    Set::Scalar W(const Set::Matrix & F) const 
    {
        return PseudoLinear::Orthotropic::W(F * Set::reduce(Fp).inverse());
    }
    Set::Matrix DW(const Set::Matrix & F) const 
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        return PseudoLinear::Orthotropic::DW(F * Fpinv) * Fpinv.transpose();
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const 
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        auto DDW = PseudoLinear::Orthotropic::DDW(F * Fpinv);
        auto ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            for (int J = 0; J < AMREX_SPACEDIM; J++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    for (int L = 0; L < AMREX_SPACEDIM; L++)
                    {
                        ret(i,J,k,L) = 0.0;
                        for (int Q = 0; Q < AMREX_SPACEDIM; Q++)
                            for (int R = 0; R < AMREX_SPACEDIM; R++)
                                ret(i,J,k,L) += DDW(i,Q,k,R)*Fpinv(J,Q)*Fpinv(L,R);
                    }
        return ret;
    }
    
    PlasticRate Rate(const Set::Matrix & F,
                     const Set::Matrix3d & a_Fp,
                     Set::Scalar time,
                     Set::Scalar temperature) const
    {
        PlasticRate rate;
        if (time < tstart) return rate;

        // With F = Fe Fp and Fp_dot = Lp Fp, the stress power that is
        // work-conjugate to Lp is the Mandel stress M = Fe^T Pe, where
        // Pe = dW/dFe. Slip systems are defined in this intermediate
        // configuration, so their resolved shear is s dot M m.
        const Set::Matrix Fpinv = Set::reduce(a_Fp).inverse();
        const Set::Matrix Fe = F * Fpinv;
        const Set::Matrix Pe = PseudoLinear::Orthotropic::DW(Fe);
        const Set::Matrix3d mandel = Set::expand(Fe.transpose() * Pe);

        Set::Quaternion normalized_q = q;
        normalized_q.normalize();
        Set::Matrix3d R = normalized_q.toRotationMatrix();
        std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ss = slipSystems(R);

        Set::Matrix3d L = Set::Matrix3d::Zero();

        for (int n = 0; n < 12; n++)
        {
            Set::Vector3d A = ss[n].first;
            Set::Vector3d N = ss[n].second;
            Set::Scalar tau = A.dot(mandel*N);
            Set::Scalar sign = tau>0 ? 1.0 : -1.0;

            // Apply local thermal softening to the initial CRSS. Mechanical
            // self/latent hardening remains additive.
            const Set::Scalar temperature_increment =
                Util::Max(temperature - temperature_reference, 0.0);
            const Set::Scalar thermal_factor = Util::Max(
                minimum_crss_fraction,
                std::exp(-thermal_softening_coefficient *
                         temperature_increment));
            Set::Scalar resistance = tau_crss[n] * thermal_factor +
                hardening_self * Util::Max(gamma[n], 0.0);
            for (int beta = 0; beta < 12; ++beta)
                if (beta != n)
                    resistance += hardening_latent *
                        Util::Max(gamma[beta], 0.0);
            rate.gamma[n] = gammadot0 *
                std::pow(fabs(tau/resistance),m_rate_inv) * sign;
            
            L += rate.gamma[n] * A * N.transpose();
        }

        rate.Fp = L * a_Fp;
        return rate;
    }

    void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix /*P*/, Set::Scalar time)
    {
        PlasticRate rate = Rate(F, Fp, time, temperature_reference);
        Fp += rate.Fp * dt;
        for (int n = 0; n < 12; ++n)
            gamma(n) += Util::Abs(rate.gamma(n)) * dt;
    }


public:
    //Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;
    //static constexpr KinematicVariable kinvar = KinematicVariable::F;
    //Set::Scalar C11=1.68, C12=1.21, C44=0.75;
    //mutable Set::Quaternion q = Set::Quaternion(1.0, 0.0, 0.0, 0.0);

    // Plastic state variables
    Set::Matrix3d Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,12,1> gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Eigen::Matrix<Set::Scalar,3,12> slip_directions = FCCDirections();
    Eigen::Matrix<Set::Scalar,3,12> slip_normals = FCCNormals();

    // Parameters
    Set::Scalar m_rate_inv = NAN;
    //Set::Scalar tau_crss = NAN;
    Eigen::Matrix<Set::Scalar,12,1> tau_crss = Eigen::Matrix<Set::Scalar,12,1>::Zero();
    Set::Scalar gammadot0 = NAN;
    Set::Scalar hardening_self = 0.0;
    Set::Scalar hardening_latent = 0.0;
    Set::Scalar temperature_reference = 300.0;
    Set::Scalar thermal_softening_coefficient = 0.0;
    Set::Scalar minimum_crss_fraction = 1.0;

    // Control variables
    Set::Scalar tstart = NAN;

private:
    using SlipMatrix = Eigen::Matrix<Set::Scalar,3,12>;

    static SlipMatrix FCCDirections()
    {
        SlipMatrix a;
        a.col(0)=Set::Vector3d( 1, 0,-1); a.col(1)=Set::Vector3d( 0,-1, 1);
        a.col(2)=Set::Vector3d( 1,-1, 0); a.col(3)=Set::Vector3d( 1, 0,-1);
        a.col(4)=Set::Vector3d( 1, 1, 0); a.col(5)=Set::Vector3d( 0, 1, 1);
        a.col(6)=Set::Vector3d( 1, 1, 0); a.col(7)=Set::Vector3d( 0,-1, 1);
        a.col(8)=Set::Vector3d( 1, 0, 1); a.col(9)=Set::Vector3d( 0, 1, 1);
        a.col(10)=Set::Vector3d(1, 0, 1); a.col(11)=Set::Vector3d(1,-1, 0);
        for (int n=0; n<12; ++n) a.col(n).normalize();
        return a;
    }
    static SlipMatrix FCCNormals()
    {
        SlipMatrix n;
        for (int i=0;i<3;++i) n.col(i)=Set::Vector3d( 1, 1, 1);
        for (int i=3;i<6;++i) n.col(i)=Set::Vector3d(-1, 1,-1);
        for (int i=6;i<9;++i) n.col(i)=Set::Vector3d( 1,-1,-1);
        for (int i=9;i<12;++i) n.col(i)=Set::Vector3d(-1,-1, 1);
        for (int i=0; i<12; ++i) n.col(i).normalize();
        return n;
    }
    static void SetAPSystems(CrystalPlastic& value,
                             const std::vector<Set::Scalar>& lattice)
    {
        // Orthorhombic AP systems assembled from reported slip-plane families.
        // Directions are crystallographic directions lying in each plane.
        const int planes[12][3] = {
            {0,1,0},{0,1,0},{0,0,1},{0,0,1},{1,0,0},{1,0,0},
            {1,1,0},{-1,1,0},{2,0,1},{2,0,-1},{4,0,1},{4,0,-1}};
        const int dirs[12][3] = {
            {0,0,1},{1,0,0},{1,0,0},{0,1,0},{0,0,1},{0,1,0},
            {0,0,1},{0,0,1},{0,1,0},{0,1,0},{0,1,0},{0,1,0}};
        for (int n=0; n<12; ++n)
        {
            value.slip_normals.col(n) = Set::Vector3d(
                planes[n][0]/lattice[0], planes[n][1]/lattice[1],
                planes[n][2]/lattice[2]);
            value.slip_directions.col(n) = Set::Vector3d(
                dirs[n][0]*lattice[0], dirs[n][1]*lattice[1],
                dirs[n][2]*lattice[2]);
            value.slip_normals.col(n).normalize();
            value.slip_directions.col(n).normalize();
        }
    }
    std::array<std::pair<Set::Vector3d,Set::Vector3d>,12>
    slipSystems(const Set::Matrix3d& R) const
    {
        std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ret;
        for (int n=0; n<12; ++n)
        {
            ret[n].first = R*slip_directions.col(n);
            ret[n].second = R*slip_normals.col(n);
        }
        return ret;
    }

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.C11=ret.C22=ret.C33=0.0;
        ret.C12=ret.C13=ret.C23=0.0;
        ret.C44=ret.C55=ret.C66=0.0;
        ret.q = Set::Quaternion(0.0,0.0,0.0,0.0);
        ret.Fp = Set::Matrix3d::Zero();
        ret.gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
        ret.slip_directions = SlipMatrix::Zero();
        ret.slip_normals = SlipMatrix::Zero();
        ret.m_rate_inv = 0.0;
        ret.tau_crss = Eigen::Matrix<Set::Scalar,12,1>::Zero();
        ret.gammadot0 = 0.0;
        ret.hardening_self = 0.0;
        ret.hardening_latent = 0.0;
        ret.temperature_reference = 300.0;
        ret.thermal_softening_coefficient = 0.0;
        ret.minimum_crss_fraction = 1.0;
        ret.tstart = 0.0;
        return ret;
    }
    static CrystalPlastic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static CrystalPlastic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        CrystalPlastic ret = CrystalPlastic::Zero();
        ret.C11 = ret.C22 = ret.C33 = C11;
        ret.C12 = ret.C13 = ret.C23 = C12;
        ret.C44 = ret.C55 = ret.C66 = C44;
        ret.q = Set::Quaternion::UnitRandom();
        ret.Fp = Set::Matrix3d::Identity();
        ret.slip_directions = FCCDirections();
        ret.slip_normals = FCCNormals();
        return ret;
    }
    static void Parse(CrystalPlastic & value, IO::ParmParse & pp)
    {
        PseudoLinear::Orthotropic::Parse(value,pp);

        std::string slip_systems;
        pp_query_default("slip_systems",slip_systems,std::string("fcc"));
        if (slip_systems == "fcc")
        {
            value.slip_directions = FCCDirections();
            value.slip_normals = FCCNormals();
        }
        else if (slip_systems == "ap_orthorhombic")
        {
            std::vector<Set::Scalar> lattice;
            pp_queryarr_default("lattice_parameters",lattice,
                                std::vector<Set::Scalar>{9.231,5.813,7.453});
            if (!IO::ParmParse::InTraversalMode() &&
                (lattice.size()!=3 || lattice[0]<=0.0 ||
                 lattice[1]<=0.0 || lattice[2]<=0.0))
                Util::Exception(INFO,"lattice_parameters requires three positive values");
            if (lattice.size()==3) SetAPSystems(value,lattice);
        }
        else if (slip_systems == "custom")
        {
            std::vector<Set::Scalar> directions, normals;
            pp_queryarr("slip_directions",directions);
            pp_queryarr("slip_normals",normals);
            if (!IO::ParmParse::InTraversalMode() &&
                (directions.size()!=36 || normals.size()!=36))
                Util::Exception(INFO,
                    "custom slip_directions and slip_normals require 36 values each");
            if (directions.size()==36 && normals.size()==36)
                for (int n=0;n<12;++n)
                    for (int d=0;d<3;++d)
                    {
                        value.slip_directions(d,n)=directions[3*n+d];
                        value.slip_normals(d,n)=normals[3*n+d];
                    }
        }
        else
            Util::Exception(INFO,"Unknown crystal slip_systems selection: ",slip_systems);

        if (!IO::ParmParse::InTraversalMode())
            for (int n=0;n<12;++n)
            {
                const Set::Scalar anorm=value.slip_directions.col(n).norm();
                const Set::Scalar nnorm=value.slip_normals.col(n).norm();
                if (!(anorm>0.0) || !(nnorm>0.0))
                    Util::Exception(INFO,"Slip directions and normals must be nonzero");
                value.slip_directions.col(n) /= anorm;
                value.slip_normals.col(n) /= nnorm;
                if (std::abs(value.slip_directions.col(n).dot(
                                 value.slip_normals.col(n))) > 1.0e-10)
                    Util::Exception(INFO,"Slip direction is not contained in its slip plane");
            }

        std::vector<Set::Scalar> tau_crss;
        // Critical resolved shear stress :math:`\tau_{crss}`
        pp_queryarr("tau_crss",tau_crss);

        // Rate hardening coefficient :math:`\dot{\gamma}_0`
        pp_query_default("gammadot0",value.gammadot0,1.0);
        pp_query_default("hardening_self",value.hardening_self,0.0);
        pp_query_default("hardening_latent",value.hardening_latent,0.0);
        pp_query_default("temperature_reference",value.temperature_reference,
                         "300.0_K",Unit::Temperature());
        pp_query_default("thermal_softening_coefficient",
                         value.thermal_softening_coefficient,0.0);
        pp_query_default("minimum_crss_fraction",
                         value.minimum_crss_fraction,1.0);

        // Inverse of the hardening exponent :math:`\frac{1}{m}`
        pp_query_default("m_rate_inv",value.m_rate_inv,0.5);
        
        // Time to activate plastic slip
        pp_query_default("tstart",value.tstart,0.0);    

        if (IO::ParmParse::InTraversalMode()) return;

        if (tau_crss.size() == 1)
            for (int n = 0; n < 12; n++) value.tau_crss[n] = tau_crss[0];
        else if (tau_crss.size() == 12)
            for (int n = 0; n < 12; n++) value.tau_crss[n] = tau_crss[n];
        else
            Util::Exception(INFO,"Expected either 1 or 12 values for tau_crss but got ",tau_crss.size());

        if (!(value.gammadot0 >= 0.0) || !std::isfinite(value.gammadot0) ||
            !(value.m_rate_inv > 0.0) || !std::isfinite(value.m_rate_inv) ||
            !(value.hardening_self >= 0.0) || !std::isfinite(value.hardening_self) ||
            !(value.hardening_latent >= 0.0) || !std::isfinite(value.hardening_latent) ||
            !(value.temperature_reference > 0.0) ||
            !std::isfinite(value.temperature_reference) ||
            !(value.thermal_softening_coefficient >= 0.0) ||
            !std::isfinite(value.thermal_softening_coefficient) ||
            !(value.minimum_crss_fraction > 0.0) ||
            !(value.minimum_crss_fraction <= 1.0) ||
            !std::isfinite(value.minimum_crss_fraction) ||
            !std::isfinite(value.tstart))
            Util::Exception(INFO, "Invalid crystal plasticity rate parameters");
        for (int n = 0; n < 12; ++n)
            if (!(value.tau_crss[n] > 0.0) || !std::isfinite(value.tau_crss[n]))
                Util::Exception(INFO, "tau_crss values must be finite and positive");

    }

#define OP_CLASS CrystalPlastic
#define OP_VARS X(C11) X(C22) X(C33) X(C12) X(C13) X(C23) X(C44) X(C55) X(C66) X(q) X(Fp) X(gamma) X(slip_directions) X(slip_normals) X(m_rate_inv) X(tau_crss) X(gammadot0) X(hardening_self) X(hardening_latent) X(temperature_reference) X(thermal_softening_coefficient) X(minimum_crss_fraction) X(tstart)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"





}
}
}


template<>
inline int Set::Field<Model::Solid::Finite::CrystalPlastic>::NComp() const 
{
    return 12;
}

template<>
inline std::string Set::Field<Model::Solid::Finite::CrystalPlastic>::Name(int i) const 
{
    if (i==0) return name + "gamma1";
    if (i==1) return name + "gamma2";
    if (i==2) return name + "gamma3";
    if (i==3) return name + "gamma4";
    if (i==4) return name + "gamma5";
    if (i==5) return name + "gamma6";
    if (i==6) return name + "gamma7";
    if (i==7) return name + "gamma8";
    if (i==8) return name + "gamma9";
    if (i==9) return name + "gamma10";
    if (i==10) return name + "gamma11";
    if (i==11) return name + "gamma12";
    return name;
}

template<>
inline void Set::Field<Model::Solid::Finite::CrystalPlastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Finite::CrystalPlastic> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < 12; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k).gamma(n);
                });
            }
        }
    }    
}



#endif
