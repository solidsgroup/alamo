#ifndef MODEL_PROPELLANT_CONSTANT_H
#define MODEL_PROPELLANT_CONSTANT_H

#include "Propellant.H"

namespace Model
{
namespace Propellant
{
class Constant
{
public:
    static constexpr const char* name = "constant";
    
    Constant() = default;
    Constant(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }


    template<int I = 0>
    AMREX_GPU_HOST
    Set::Scalar get (std::string /*property*/)
    {
        Util::Abort(INFO);
        return NAN;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void set_pressure(const Set::Scalar /*P*/) 
    {}
    
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_K (const Set::Scalar /*phi*/) 
    {return this->K;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_rho (const Set::Scalar /*phi*/) 
    {return this->rho;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_cp (const Set::Scalar /*phi*/) 
    {return this->cp;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_qdot (const Set::Scalar /*mdot*/, const Set::Scalar /*phi*/) 
    {return this->qdot;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_L (Set::Scalar /*phi*/, Set::Scalar /*T*/) 
    {return L;}

    static void Parse(Constant& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
        pp_query_required("K", value.K);   
        pp_query_required("rho", value.rho);
        pp_query_required("cp", value.cp);
        pp_query_required("qdot", value.qdot);
        pp_query_required("L", value.L);
    }

    Set::Scalar K = NAN;
    Set::Scalar rho = NAN;
    Set::Scalar cp = NAN;
    Set::Scalar qdot = NAN;
    Set::Scalar L = NAN;
};
}
}

#endif
