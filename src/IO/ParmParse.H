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
// .. _query-directives:
//
// **Query directives**
//
// Alamo uses different query directives to standardize reading inputs and automatic documentation.
// The type of query that is used (query vs query_required, etc) triggers different handlers and different
// automatic documentation procedures.
// For instance, the use of a :code:`query_default` causes a default value to be set and added to the metadata
// file, and it also serves as a flag for the autodoc system to document that a default value is available.
// The following table is a reference for the different kinds of query directives.
//
//
// .. table::
//    :widths: 1 99
//
//    +---------------------------------+-----------------------------------------------------------------+
//    | BC Type                         | Description                                                     |
//    +=================================+=================================================================+
//    | :bdg-warning:`query`            | Standard IO for bool, string, Set::Scalarthat does not enforce  |
//    |                                 | defaults or required values. Not recommended for general use.   |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_required`   | Similar to query, but will abort if no value is specified.      |
//    |                                 | Required values are indicated by :bdg-danger-line:`required`.   |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_default`    | Similar to query, but will fill with default value if no value  |
//    |                                 | is provided. Will also add default value to metadata.           |
//    |                                 | Default values are indicated by green badge values, e.g.        |
//    |                                 | :bdg-success:`0.0`.                                             |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_validate`   | For strings, read in a value and enforce that the value is one  |
//    |                                 | of a supplied number of values. Optionally make required, or    |
//    |                                 | set the default value to the first supplied value.              | 
//    |                                 | Acceptable options are indicated by blue badge values, e.g.     | 
//    |                                 | :bdg-primary:`1.0`, :bdg-primary:`2.0`. If a default value is   | 
//    |                                 | available, it is indicated by a green badge, e.g.               | 
//    |                                 | :bdg-success:`1.0`.                                             | 
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_file`       | Read in a string that defines a file name.                      |
//    |                                 | Check to make sure that the file exists and is a regular file,  |
//    |                                 | and print an informative error message if not (this can be      |
//    |                                 | disabled).                                                      |
//    |                                 | Also, copy the file to the output directory, with the full path |
//    |                                 | preserved by replacing / with _.                                |
//    |                                 | (This can also be disabled, but doing so is discouraged.)       |
//    |                                 | Default values are not allowed.                                 |
//    |                                 | File paths are indicated by :bdg-secondary-line:`file path`.    |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-primary:`queryarr`         | Read in an array of numbers into either a standard vector or    |
//    |                                 | into a :code:`Set::Vector` or :code:`Set::Scalar`.              |
//    |                                 | No defaults or existence checking is performed.                 |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-primary:`queryclass`       | Read a class object with a specified prefix.                    |
//    |                                 | How that class is read in is determined by its :code:`Parse`    |
//    |                                 | function.                                                       |
//    +---------------------------------+-----------------------------------------------------------------+
// 
// .. _query_locator_macros:
//
// **Query macros**
//
// A set of preprocessor macros are defined so that you can call a query function using :code:`pp_` instead
// of :code:`pp.`.
// For instance, the following two can be used interchangeably:
//
// .. code-block:: cpp
//
//     pp.query_required("myvar",myvar); /* function version */
//     pp_query_required("myvar",myvar); /* preprocessor macro version - preferred*/
// 
// Using the preprocessor macros enables code location information to be passed to the parser, so that
// more informative error messages will be printed out.
// Note that **the ParmParse object must always be called** :code:`pp` **for this to work**.
// 
//

#ifndef IO_PARMPARSE
#define IO_PARMPARSE

#include <filesystem>
#include <exception>

#include "Util/Util.H"
#include "Set/Set.H"
#include "AMReX_ParmParse.H"


#define pp_query_required(...) pp.query_required(__VA_ARGS__,INFO)
#define pp_query_default(...) pp.query_default(__VA_ARGS__,INFO)
#define pp_query_validate(...) pp.query_validate(__VA_ARGS__,INFO)
#define pp_query_file(...) pp.query_file(__VA_ARGS__,INFO)
#define pp_queryarr(...) pp.queryarr(__VA_ARGS__,INFO)
#define pp_query(...) pp.query(__VA_ARGS__)
#define pp_queryclass(...) pp.queryclass(__VA_ARGS__,INFO)


namespace IO
{
class ParmParse : public amrex::ParmParse
{
private:
    void Define() 
    {
        if (checked_for_input_files) return;
        int k = 0;
        std::string inputfile = "";
        while (this->querykth("input",k,inputfile)) 
        {
            Util::Message(INFO,"Including inputs from "+inputfile);
            this->addfile(inputfile);
            k++;
        }
        checked_for_input_files = true;
    }
    static bool checked_for_input_files;
    
public:
    ParmParse(std::string arg) : amrex::ParmParse::ParmParse(arg) {Define();} ;
    ParmParse() : amrex::ParmParse::ParmParse() {Define();} ;
    std::string getPrefix() const {return amrex::ParmParse::getPrefix();};
    void ignore(std::string name)
    {
        (void)amrex::ParmParse::contains(name.c_str());
    }

    template<typename T> 
    int query_required( std::string name, T & value,
                        std::string file = "", std::string func = "", int line = -1)
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,file,func,line,full(name),"required value for ",full(name)," missing");
        }
        return query(name.c_str(),value);
    }

    template<typename T>
    int query_default(  std::string name, T & value, T defaultvalue,
                        std::string = "", std::string = "", int = -1)
    {
        if (!contains(name.c_str()))
        {
            add(name.c_str(),defaultvalue);
        }
        return query(name.c_str(),value);
    } 

    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals, bool firstbydefault,
                        std::string file = "", std::string func = "", int line = -1)
    {        
        // if not using default values, then the input must be specified
        if (!firstbydefault)
        {
            if (!contains(name.c_str()))
            {
                Util::ParmParseException(INFO,file,func,line,full(name),"required value for ",full(name)," missing");
            }
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

        Util::ParmParseException(INFO,file,func,line,full(name),"' expected [", ss.str(), "] but got ", value);

        return -1;
    }

    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals,
                        std::string file = "", std::string func = "", int line = -1)
    {
        return query_validate(name,value,possiblecharvals,true,file,func,line);
    }

    
    // special case for strings
    int query_default(  std::string name, std::string & value, const char *defaultvalue,
                        std::string file = "", std::string func = "", int line = -1)
    {
        return query_default(name, value, std::string(defaultvalue), file, func, line);
    }
    // special case for bools
    int query_default(  std::string name, int & value, bool defaultvalue,
                        std::string file = "", std::string func = "query_default", int line = -1)
    {
        int defaultint = 0;
        if (defaultvalue) defaultint = 1;
        return query_default(name, value, defaultint, file, func, line);
    }



    // validate filenames
    int query_file( std::string name, std::string & value, bool copyfile, bool checkfile,
                    std::string file = "", std::string func = "query_file", int line = -1)
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,file,func,line,full(name),full(name)," must be specified");
        }

        int retval = query(name.c_str(),value);

        if ( checkfile && ! std::filesystem::exists(value))
        {
            Util::ParmParseException(INFO,file,func,line,full(name),full(name)," does not exist");
        }

        if ( checkfile && !std::filesystem::is_regular_file(value))
        {
            Util::ParmParseException(INFO,file,func,line,full(name),full(name)," is not a regular file");
        }

        if ( copyfile )
        {
            Util::CopyFileToOutputDir(value, true);
        }

        return retval;
    }
    int query_file( std::string name, std::string & value, bool copyfile, 
                    std::string file = "", std::string func = "query_file", int line = -1)
    {
        return query_file(name,value,copyfile,true,file,func,line);
    }
    int query_file( std::string name, std::string & value, 
                    std::string file = "", std::string func = "query_file", int line = -1)
    {
        return query_file(name,value,true,true,file,func,line);
    }


    template<typename T>
    int queryarr(   std::string name, std::vector<T> & value,
                    std::string /*file*/, std::string /*func*/, int /*line*/)
    {
        return amrex::ParmParse::queryarr(name.c_str(),value);
    }
    int queryarr(   std::string name, Set::Vector & value,
                    std::string file = "", std::string func = "queryarr", int line = -1)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() < AMREX_SPACEDIM) 
        {
            Util::ParmParseException(INFO,file,func,line,full(name),full(name)," requires at least ", AMREX_SPACEDIM, " arguments, got ",vals.size());
        }
        for (int i = 0; i < AMREX_SPACEDIM; i++) value(i) = vals[i];
        return 0;
    }
    int queryarr(   std::string name, Set::Matrix & value,
                    std::string file = "", std::string func = "queryarr", int line = -1)
    {
        std::vector<Set::Scalar> vals;
        amrex::ParmParse::queryarr(name.c_str(), vals);
        if (vals.size() == 9)
        {
#if AMREX_SPACEDIM==2
            Util::Warning(file,func,line,"Reading a 3D matrix (",full(name),")into a 2D code - some values will be ignored.");
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
            Util::Warning(file,func,line,"Reading a 2D matrix (",full(name),")into a 3D code - remaining values will be set to zero.");
            value(0,0) = vals[0];  value(0,1)= vals[1];  value(0,2)= 0.0;
            value(1,0) = vals[2];  value(1,1)= vals[3];  value(1,2)= 0.0;
            value(2,0) = 0.0;      value(2,1)= 0.0;      value(2,2)= 0.0;
#endif
        }
        else
        {
            Util::ParmParseException(INFO,file,func,line,full(name),full(name)," needs either 4 or 9 components, but got ",vals.size());
        }    
        return 0;
    }

    int AnyUnusedInputs()
    {
        int cnt = 0;
        for (std::list<amrex::ParmParse::PP_entry>::const_iterator li = m_table->begin(), End = m_table->end(); li != End; ++li)
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
        for (std::list<amrex::ParmParse::PP_entry>::const_iterator li = pp.m_table->begin(), End = pp.m_table->end(); li != End; ++li)
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
    std::string full (std::string name)
    {
        std::string prefix = getPrefix();
        if (prefix != "") return getPrefix() + "." + name;
        else return name;
    }


    using amrex::ParmParse::queryarr;
    template<class T>
    void queryclass(std::string name, T * value,
                    std::string file = "", std::string func = "", int line = -1)
    {
        try
        {
            this->pushPrefix(name);
            T::Parse(*value, *this);
            if (AnyUnusedInputs())
            {
                Util::ParmParseException(INFO,file,func,line,full(name),"Above inputs specified but not used");
            }
        }
        catch (...)
        {
            this->popPrefix();
            Util::ParmParseException(INFO,file,func,line,full(name));
        }
        this->popPrefix();
    }
    template<class T>
    void queryclass(std::string name, T & value,
                    std::string file = "", std::string func = "", int line = __LINE__)
    {
        try
        {
            this->pushPrefix(name);
            T::Parse(value, *this);
            if (AnyUnusedInputs())
            {
                Util::ParmParseException(INFO,file,func,line,full(name),"Above inputs specified but not used");
            }
        }
        catch (...)
        {
            this->popPrefix();
            Util::ParmParseException(INFO,file,func,line,full(name));
        }
        this->popPrefix();
    }
    template<class T>
    void queryclass(T * value,
                    std::string file = "", std::string func = "", int line = __LINE__)
    {
        try
        {
            T::Parse(*value, *this);
            if (AnyUnusedInputs())
            {
                Util::ParmParseException(INFO,file,func,line,"Above inputs specified but not used");
            }
        }
        catch (...)
        {
            Util::ParmParseException(INFO,file,func,line);
        }
    }
    template<class T>
    void queryclass(T & value,
                    std::string file = "", std::string func = "", int line = __LINE__)
    {
        try
        {
            T::Parse(value, *this);
            if (AnyUnusedInputs())
            {
                Util::ParmParseException(INFO,file,func,line,"Above inputs specified but not used");
            }
        }
        catch (...)
        {
            Util::ParmParseException(INFO,file,func,line);
        }
    }
};
}
#endif
