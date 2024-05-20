//
// This is a thin wrapper to the amrex::ParmParse class
// This class exists to add some additional parsing capability, 
// e.g. parsing Set::Matrix and Set::Vector data types.
//
// :ref:`IO::ParmParse` uses static :code:`Parse()` functions to 
// perform cascading class-based input parsing.
// See the :ref:`autodoc` section for instructions on adding documentation.
// 
// :bdg-warning-line:`This is standard infrastructure code; make sure you know what you ard doing before you change it.`
//

#ifndef IO_PARMPARSE
#define IO_PARMPARSE
#include "Set/Set.H"
#include "AMReX_ParmParse.H"

namespace IO
{
class ParmParse : public amrex::ParmParse
{
public:
    ParmParse(std::string arg) : amrex::ParmParse::ParmParse(arg) {} ;
    ParmParse() : amrex::ParmParse::ParmParse() {} ;
    std::string getPrefix() const {return amrex::ParmParse::getPrefix();};
    void ignore(std::string name)
    {
        (void)amrex::ParmParse::contains(name.c_str());
    }
    template<typename T>
    int query_required(std::string name, T & value)
    {
        if (!contains(name.c_str())) Util::Abort(INFO,name," must be specified");
        return query(name.c_str(),value);
    }
    template<typename T>
    int query_default(std::string name, T & value, T defaultvalue)
    {
        if (!contains(name.c_str()))
        {
            add(name.c_str(),defaultvalue);
        }
        return query(name.c_str(),value);
    } 
    int query_validate(std::string name, std::string & value, std::vector<const char *> possiblecharvals, bool firstbydefault = true)
    {
        
        // if not using default values, then the input must be specified
        if (!firstbydefault)
        {
            if (!contains(name.c_str())) Util::Abort(INFO,name," must be specified");
        }

        // set default value
        value = std::string(possiblecharvals[0]);

        // get the read value (if it exists)
        int retval = query(name.c_str(),value);

        // check to make sure the read value matches one of the inpus
        bool ok = false;
        for (unsigned int i = 0; i < possiblecharvals.size(); i++)
        {
            if (value == std::string(possiblecharvals[i])) ok = true;
        }
        
        if (ok) return retval;

        std::stringstream ss;
        ss << possiblecharvals[0];
        for (unsigned int i = 1; i < possiblecharvals.size(); i++)
            ss << "," << possiblecharvals[i];

        Util::Abort(INFO,"Reading '",name,"' expected [", ss.str(), "] but got ", value);

        return -1;
    }
    
    // special case for strings
    int query_default(std::string name, std::string & value, const char *defaultvalue)
    {
        return query_default(name, value, std::string(defaultvalue));
    }
    // special case for bools
    int query_default(std::string name, int & value, bool defaultvalue)
    {
        int defaultint = 0;
        if (defaultvalue) defaultint = 1;
        return query_default(name, value, defaultint);
    }

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
            if (!li->m_queried && li->m_name.rfind(getPrefix()+".",0) != std::string::npos)
            {
                Util::Warning(INFO,li->m_name);
                cnt++;
            }
        }
        return cnt;
    }

    static int AllUnusedInputs()
    {
        ParmParse pp;
        int cnt = 0;
        for (std::list<amrex::ParmParse::PP_entry>::const_iterator li = pp.m_table.begin(), 
            End = pp.m_table.end(); 
            li != End; ++li)
        {
            if (!li->m_queried)
            {
                Util::Warning(INFO,li->m_name);
                cnt++;
            }
        }
        return cnt;
    }
    std::string prefix ()
    {
        return getPrefix();
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
