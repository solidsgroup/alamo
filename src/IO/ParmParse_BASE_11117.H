#include "Set/Set.H"
#include "AMReX_ParmParse.H"

namespace IO
{
///
/// \brief This is a thin wrapper to the amrex::ParmParse class
///
/// This class exists to add some additional parsing capability, 
/// e.g. parsing Set::Matrix and Set::Vector data types.
/// Otherwise, it is used in exactly the same way as amrex::ParmParse
///
class ParmParse : public amrex::ParmParse
{
public:
    ParmParse(std::string arg) : amrex::ParmParse::ParmParse(arg) {} ;
    void queryarr(std::string name, Set::Vector & value)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() < AMREX_SPACEDIM) 
        {
            Util::Abort(INFO,name," requires at least ", AMREX_SPACEDIM, " arguments, got ",vals.size());
        }
        for (int i = 0; i < AMREX_SPACEDIM; i++) value(i) = vals[i];
    }
    void queryarr(std::string name, Set::Matrix & value)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() < AMREX_SPACEDIM*AMREX_SPACEDIM) 
        {
            Util::Abort(INFO,name," requires at least ", AMREX_SPACEDIM*AMREX_SPACEDIM, " arguments, got ",vals.size());
        }    
        for (int i = 0; i < AMREX_SPACEDIM*AMREX_SPACEDIM; i++) 
        {
            value(i/AMREX_SPACEDIM, i%AMREX_SPACEDIM) = vals[i];
        }
    }

};
}