#ifndef MODEL_SOLID_LINEAR_CUBIC_H_
#define MODEL_SOLID_LINEAR_CUBIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Cubic : public Solid<Set::Sym::MajorMinor>
{
public:

    Cubic() {};
    virtual ~Cubic() {};

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (DDW(gradu)*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return DDW(gradu)*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & /*gradu*/) const override
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> ddw;
    

        Eigen::Matrix3d R = q.normalized().toRotationMatrix();

        for(int p = 0; p < AMREX_SPACEDIM; p++) 
            for(int q = 0; q < AMREX_SPACEDIM; q++) 
                for(int s = 0; s < AMREX_SPACEDIM; s++) 
                    for(int t = 0; t < AMREX_SPACEDIM; t++)
                    {
                        ddw(p,q,s,t) = 0.0;
                        for(int i = 0; i < 3; i++) 
                            for(int j = 0; j < 3; j++) 
                                for(int k = 0; k < 3; k++) 
                                    for(int l = 0; l < 3; l++) 
                                    {
                                        if(i == j && j == k && k == l) ddw(p,q,s,t) += R(p,i)*R(s,k)*C11*R(q,j)*R(t,l);
                                        else if (i==k && j==l)         ddw(p,q,s,t) += R(p,i)*R(s,k)*C44*R(q,j)*R(t,l);
                                        else if (i==j && k==l)         ddw(p,q,s,t) += R(p,i)*R(s,k)*C12*R(q,j)*R(t,l);
                                        else continue;
                                    }
                    }
        return ddw;
    }

    virtual void Print(std::ostream &out) const //override 
    {
        out << "C11=" << C11 << ", C12=" << C12 << ", C44=" << C44 << " R=" << q.toRotationMatrix() << std::endl;
        // out << "C110=" << C110 << ", C120=" << C120 << ", C440=" << C440 << std::endl;
    }

    void DegradeModulus(const Set::Scalar eta)
    {
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
        // if(std::isnan(C110)) Util::Abort(INFO, "Nans in C110");
        // if(std::isnan(C120)) Util::Abort(INFO, "Nans in C120");
        // if(std::isnan(C440)) Util::Abort(INFO, "Nans in C440");

        // C11 = C110*(1.-eta);
        // C12 = C120*(1.-eta);
        // C44 = C440*(1.-eta);
    }
    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
    {
        if(eta.size() == 0) return;
        else if (eta.size() == 1) DegradeModulus(eta[0]);
        else if(eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            // C11 = C110*(1.-eta[0]);
            // C12 = C120*(1.-eta[1]);
            // C44 = C440*(1.-eta[1]);
        }
        else if (eta.size() == 3)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0 || eta[2] < 0.0 || eta[2] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            // C11 = C110*(1.-eta[0]);
            // C12 = C120*(1.-eta[1]);
            // C44 = C440*(1.-eta[2]);
        }
        else
            return;
    }

public:

    Set::Scalar C11=NAN, C12=NAN, C44=NAN;
    //Set::Scalar C110=NAN, C120=NAN, C440=NAN;
    Eigen::Quaterniond q;    

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;

    AMREX_FORCE_INLINE
    static Cubic Combine(const std::vector<Cubic> &models, const std::vector<Set::Scalar> &eta)
    {
        Cubic ret = Cubic::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret += models[n] * (eta[n] / etasum);
        }
        return ret;
    }

    static Cubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static Cubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        Cubic ret;
        ret.C11 = C11; ret.C12 = C12; ret.C44 = C44;
        // ret.C110 = C11; ret.C120 = C12; ret.C440 = C44;
        ret.q = Eigen::Quaterniond::UnitRandom();
        return ret;
    }
    static Cubic Zero()
    {
        Cubic ret;
        ret.q.coeffs() *= 0.0; 
        ret.C11 = 0.0;
        ret.C12 = 0.0;
        ret.C44 = 0.0;
        // ret.C110 = 0.0;
        // ret.C120 = 0.0;
        // ret.C440 = 0.0;
        return ret;
    }

    static void Parse(Cubic & value, IO::ParmParse & pp)
    {
        //Set::Scalar C11 = 1.68, C12 = 1.21, C44 = 0.75;
        pp.query("C11",value.C11);
        pp.query("C12",value.C12);
        pp.query("C44",value.C44);
        // value.C110 = value.C11;
        // value.C120 = value.C12;
        // value.C440 = value.C44;

        Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
        pp.query("phi1",phi1);
        pp.query("Phi",Phi);
        pp.query("phi2",phi2);

        Eigen::Matrix3d R;
        R = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
            Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
        value.q = Eigen::Quaterniond(R);

        if (pp.contains("random"))
        {
            value.q = Eigen::Quaterniond::UnitRandom();
        }

        if (pp.contains("phi1") || pp.contains("Phi") || pp.contains("phi2"))
        {
            pp.query("phi1",phi1);
            pp.query("Phi",Phi);
            pp.query("phi2",phi2);

            Eigen::Matrix3d R;
            R = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
                Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
            value.q = Eigen::Quaterniond(R);
        }    
        // std::cout << "q = " << value.q.w() << " + " << value.q.vec() << std::endl;
        // std::cout << "R = " << std::endl << value.q.normalized().toRotationMatrix() << std::endl;
    }    

    #define OP_CLASS Cubic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) /*X(C110) X(C120) X(C440)*/ 
    #include "Model/Solid/InClassOperators.H"
    
    
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
