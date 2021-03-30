#ifndef IO_PARMPARSE
#define IO_PARMPARSE
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
    ParmParse() : amrex::ParmParse::ParmParse() {} ;
    int queryarr(std::string name, Set::Vector & value)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() < AMREX_SPACEDIM) 
        {
            Util::Abort(INFO,name," requires at least ", AMREX_SPACEDIM, " arguments, got ",vals.size());
        }
        for (int i = 0; i < AMREX_SPACEDIM; i++) value(i) = vals[i];
        return 0;
    }
    int queryarr(std::string name, Set::Matrix & value)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() == 9)
        {
            #if AMREX_SPACEDIM==2
            Util::Warning(INFO,"Reading a 3D matrix (",name,")into a 2D code - some values will be ignored.");
            value(0,0) = vals[0];  value(0,1)= vals[1];
            value(1,0) = vals[3];  value(1,1)= vals[4];
            #endif
            #if AMREX_SPACEDIM==3
            value(0,0) = vals[0];  value(0,1)= vals[1];  value(0,2)= vals[2];
            value(1,0) = vals[3];  value(1,1)= vals[4];  value(1,2)= vals[5];
            value(2,0) = vals[6];  value(2,1)= vals[7];  value(2,2)= vals[8];
            #endif
        }
        else if (vals.size() == 4)
        {
            #if AMREX_SPACEDIM==2
            value(0,0) = vals[0];  value(0,1)= vals[1];
            value(1,0) = vals[2];  value(1,1)= vals[3];
            #endif
            #if AMREX_SPACEDIM==3
            Util::Warning(INFO,"Reading a 2D matrix (",name,")into a 3D code - remaining values will be set to zero.");
            value(0,0) = vals[0];  value(0,1)= vals[1];  value(0,2)= 0.0;
            value(1,0) = vals[2];  value(1,1)= vals[3];  value(1,2)= 0.0;
            value(2,0) = 0.0;      value(2,1)= 0.0;      value(2,2)= 0.0;
            #endif
        }
        else
        {
            Util::Abort(INFO,name," needs either 4 or 9 components, but got ",vals.size());
        }    
        return 0;
    }

    int AnyUnusedInputs()
    {
        int cnt = 0;
        for (std::list<amrex::ParmParse::PP_entry>::const_iterator li = m_table.begin(), 
            End = m_table.end(); 
            li != End; ++li)
        {
            if (!li->m_queried && li->m_name.find(getPrefix()) != std::string::npos)
            {
                Util::Warning(INFO,li->m_name);
                cnt++;
            }
        }
        return cnt;
    }

    using amrex::ParmParse::queryarr;
    template<class T>
    void queryclass(std::string name, T * value)
    {
        this->pushPrefix(name);
        T::Parse(*value, *this);
        if (AnyUnusedInputs()) 
            Util::Abort(INFO,"DANGER - the above inputs were specified but not used!");
        this->popPrefix();
    }
    template<class T>
    void queryclass(std::string name, T & value)
    {
        this->pushPrefix(name);
        T::Parse(value, *this);
        if (AnyUnusedInputs()) 
            Util::Abort(INFO,"DANGER - the above inputs were specified but not used!");
        this->popPrefix();
    }
    template<class T>
    void queryclass(T * value)
    {
        T::Parse(*value, *this);
        if (AnyUnusedInputs()) 
            Util::Abort(INFO,"DANGER - the above inputs were specified but not used!");
    }
    template<class T>
    void queryclass(T & value)
    {
        T::Parse(value, *this);
        if (AnyUnusedInputs()) 
            Util::Abort(INFO,"DANGER - the above inputs were specified but not used!");
    }
};
}
#endif
