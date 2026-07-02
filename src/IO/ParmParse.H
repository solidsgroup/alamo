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

#include <array>
#include <filesystem>
#include <exception>
#include <list>
#include <stdexcept>

#include "Util/Util.H"
#include "Unit/Unit.H"
#include "Set/Base.H"
#include "AMReX_ParmParse.H"


#define pp_query_required(...) pp.query_required(__VA_ARGS__)
#define pp_query_default(...) pp.query_default(__VA_ARGS__)
#define pp_query_validate(...) pp.query_validate(__VA_ARGS__)
#define pp_query_file(...) pp.query_file(__VA_ARGS__)
#define pp_queryarr(...) pp.queryarr(__VA_ARGS__)
#define pp_queryarr_required(...) pp.queryarr_required(__VA_ARGS__)
#define pp_queryarr_default(...) pp.queryarr_default(__VA_ARGS__)
#define pp_query(...) pp.query(__VA_ARGS__)
#define pp_queryclass(...) pp.queryclass(__VA_ARGS__)
#define pp_forbid(...) pp.forbid(__VA_ARGS__)


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
    std::string getPrefix() const {return m_prefix;};
    void ignore(std::string name)
    {
        (void)amrex::ParmParse::contains(name.c_str());
    }
    
    void pushPrefix(const std::string prefix)
    {
        if (m_prefix.length())
            m_prefix = m_prefix + "." + prefix;
        else
            m_prefix = prefix;
    }

    void popPrefix()
    {
        size_t pos = m_prefix.rfind('.');
        if (pos == std::string::npos)
            m_prefix = "";
        else
            m_prefix = m_prefix.substr(0, pos);
    }


    void forbid(std::string name, std::string explanation)
    {
        if (amrex::ParmParse::contains(full(name).c_str()))
        {
            Util::ParmParseException(INFO,full(name),full(name)," forbidden: ", explanation);
        }
        std::set<std::string> subs = amrex::ParmParse::getEntries(full(name));
        if (subs.size())
        {
            Util::ParmParseException(INFO,full(name),full(name)," forbidden: ", explanation);
        }
    }

    bool contains(std::string name)
    {
        if (amrex::ParmParse::contains(name.c_str()))
            return true;
        if (amrex::ParmParse::contains(full(name).c_str()))
            return true;
        {
            std::set<std::string> subs = amrex::ParmParse::getEntries(name.c_str());
            if (subs.size())
                return true;
        }
        {
            std::set<std::string> subs = amrex::ParmParse::getEntries(full(name).c_str());
            if (subs.size())
                return true;
        }
        return false;
    }

    int queryunit (std::string name, Unit &value)
    {
        try
        {
            std::string strvalue;
            int retval = amrex::ParmParse::query(name.c_str(),strvalue);
            value = Unit::Parse(strvalue);
            Util::DebugMessage(INFO,full(name),": ",strvalue," ==> ",value);
            return retval;
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryunit (std::string name, Unit &value, const Unit type)
    {
        try
        {
            int retval = queryunit(name,value);
            if (!value.isType(type) && !value.isType(Unit::Less()))
            {
                Util::ParmParseException(INFO, full(name),
                                "value requiested had wrong units:", value.normalized_unitstring(),
                                ", units requested are of SI type ",type.normalized_unitstring());
            }
            return retval;
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryunit (std::string name, Set::Scalar &value, const Unit type)
    {
        try
        {
            Unit read;
            int retval = queryunit(name,read,type);
            value = read.normalized_value();
            return retval;
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        return -1;
    }

    template<typename T> 
    int query_required( std::string name, T & value)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
            }
            return query(name.c_str(),value);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T> 
    int query_required( std::string name, T & value, const Unit type)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
            }
            return queryunit(name.c_str(), value, type);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T>
    int query_default(  std::string name, T & value, T defaultvalue)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                add(name.c_str(),defaultvalue);
            }
            return query(name.c_str(),value);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    } 

    template<typename T>
    int query_default(  std::string name, T & value, std::string defaultvalue, const Unit type)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                add(name.c_str(),defaultvalue);
            }
            return queryunit(name.c_str(),value,type);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    } 

    int query_validate( std::string name, int & value, std::vector<int> possibleintvals)
    {
        try
        {
            // First value is accepted by default...

            // set default value
            if (!contains(name.c_str()))
            {
                add(name.c_str(), possibleintvals[0]);
            }

            // get the read value (if it exists)
            int retval = query(name.c_str(),value);

            // check to make sure the read value matches one of the inpus
            bool ok = false;
            for (unsigned int i = 0; i < possibleintvals.size(); i++)
            {
                if (value == possibleintvals[i]) ok = true;
            }

            if (ok) return retval;

            std::stringstream ss;
            ss << possibleintvals[0];
            for (unsigned int i = 1; i < possibleintvals.size(); i++)
                ss << "," << possibleintvals[i];

            Util::ParmParseException(INFO,full(name),"' expected [", ss.str(), "] but got ", value);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals, bool firstbydefault)
    {        
        try
        {
            // if not using default values, then the input must be specified
            if (!firstbydefault)
            {
                if (!amrex::ParmParse::contains(name.c_str()))
                {
                    Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
                }
            }

            // set default value
            if (!amrex::ParmParse::contains(name.c_str()))
            {
                add(name.c_str(), std::string(possiblecharvals[0]));
            }

            // get the read value (if it exists)
            int retval = amrex::ParmParse::query(name.c_str(),value);

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

            Util::ParmParseException(INFO,full(name),"' expected [", ss.str(), "] but got ", value);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }

        return -1;
    }

    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals)
    {
        try
        {
            return query_validate(name,value,possiblecharvals,true);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    
    // special case for strings
    int query_default(  std::string name, std::string & value, const char *defaultvalue)
    {
        try
        {
            return query_default(name, value, std::string(defaultvalue));
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    // special case for bools
    int query_default(  std::string name, int & value, bool defaultvalue)
    {
        try
        {
            int defaultint = 0;
            if (defaultvalue) defaultint = 1;
            return query_default(name, value, defaultint);
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    // validate filenames
    int query_file( std::string name, std::string & value, bool copyfile, bool checkfile)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                Util::ParmParseException(INFO,full(name),full(name)," must be specified");
            }

            int retval = query(name.c_str(),value);

            if (amrex::ParallelDescriptor::IOProcessor())
            {
                if ( checkfile && ! std::filesystem::exists(value))
                {
                    Util::ParmParseException(INFO,full(name),full(name)," does not exist");
                }
                if ( checkfile && !std::filesystem::is_regular_file(value))
                {
                    Util::ParmParseException(INFO,full(name),full(name)," is not a regular file");
                }
                if ( copyfile )
                {
                    Util::CopyFileToOutputDir(value, true, full(name));
                }
            }
            return retval;
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int query_file( std::string name, std::string & value, bool copyfile)
    {
        return query_file(name,value,copyfile,true);
    }
    int query_file( std::string name, std::string & value)
    {
        return query_file(name,value,true,true);
    }


    template<typename T>
    int queryarr(   std::string name, std::vector<T> & value)
    {
        try
        {
            return amrex::ParmParse::queryarr(name.c_str(),value);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, std::vector<Set::Scalar> & value, Unit unit = Unit::Less())
    {
        try
        {
            if (unit.isType(Unit::Less()))
            {
                return amrex::ParmParse::queryarr(name.c_str(),value);
            }
            else
            {
                value.clear();
                std::vector<std::string> valstrings;
                int retval = amrex::ParmParse::queryarr(name.c_str(), valstrings);
                for (unsigned int i = 0; i < valstrings.size(); i++)
                {
                    Unit unitvalue = Unit::Parse(valstrings[i]);
                    Util::DebugMessage(INFO,full(name),": ",valstrings[i]," ==> ",unitvalue);

                    if (!unitvalue.isType(unit) && !unitvalue.isType(Unit::Less()))
                    {
                        Util::Exception(INFO,"value specified had wrong units:", valstrings[i]);
                    }
                    value.push_back(unitvalue.normalized_value());
                }
                return retval;
            }
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, Set::Vector & value, Unit unit = Unit::Less())
    {
        try
        {
            std::vector<Set::Scalar> vals;
            queryarr(name.c_str(), vals, unit);
            if (vals.size() < AMREX_SPACEDIM) 
            {
                Util::ParmParseException(   INFO,full(name),full(name),
                                            " requires at least ", AMREX_SPACEDIM, 
                                            " arguments, got ",vals.size());
            }
            for (int i = 0; i < AMREX_SPACEDIM; i++) value(i) = vals[i];
            return 0;
        }
        catch(...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, Set::Matrix & value, Unit unit = Unit::Less())
    {
        std::vector<Set::Scalar> vals;
        queryarr(name.c_str(), vals,unit);
        if (vals.size() == 9)
        {
#if AMREX_SPACEDIM==2
            Util::Warning(INFO, "Reading a 3D matrix (",full(name),")into a 2D code - some values will be ignored.");
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
            Util::Warning(INFO,"Reading a 2D matrix (",full(name),")into a 3D code - remaining values will be set to zero.");
            value(0,0) = vals[0];  value(0,1)= vals[1];  value(0,2)= 0.0;
            value(1,0) = vals[2];  value(1,1)= vals[3];  value(1,2)= 0.0;
            value(2,0) = 0.0;      value(2,1)= 0.0;      value(2,2)= 0.0;
#endif
        }
        else
        {
            Util::ParmParseException(INFO,full(name),full(name)," needs either 4 or 9 components, but got ",vals.size());
        }    
        return 0;
    }
    int queryarr(   std::string name, std::vector<Set::Scalar> & value, int n, Unit unit = Unit::Less())
    {
        try
        {
            if (unit.isType(Unit::Less()))
            {
                return amrex::ParmParse::queryarr(name.c_str(),value,0,n);
            }
            else
            {
                value.clear();
                std::vector<std::string> valstrings;
                int retval = amrex::ParmParse::queryarr(name.c_str(), valstrings, 0, n);
                for (unsigned int i = 0; i < valstrings.size(); i++)
                {
                    Unit unitvalue = Unit::Parse(valstrings[i]);
                    Util::DebugMessage(INFO,full(name),": ",valstrings[i]," ==> ",unitvalue);

                    if (!unitvalue.isType(unit) && !unitvalue.isType(Unit::Less()))
                    {
                        Util::Exception(INFO,"value specified had wrong units:", valstrings[i]);
                    }
                    value.push_back(unitvalue.normalized_value());
                }
                return retval;
            }
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T, std::size_t n>
    int queryarr(   std::string name, std::array<T,n> & value)
    {
        try
        {
            std::vector<T> vals;
            int retval = amrex::ParmParse::queryarr(name.c_str(), vals, 0, static_cast<int>(n));
            if (retval)
            {
                for (std::size_t i = 0; i < n; i++) value[i] = vals[i];
            }
            return retval;
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<std::size_t n>
    int queryarr(   std::string name, std::array<Set::Scalar,n> & value, Unit unit = Unit::Less())
    {
        try
        {
            std::vector<Set::Scalar> vals;
            int retval = queryarr(name, vals, static_cast<int>(n), unit);
            if (retval)
            {
                for (std::size_t i = 0; i < n; i++) value[i] = vals[i];
            }
            return retval;
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T>
    int queryarr_required(  std::string name, std::vector<T> & value)
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr<T>(name,value);
    }
    int queryarr_required(  std::string name, std::vector<Set::Scalar> & value, Unit unit = Unit::Less())
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr(name,value,unit);
    }

    template<typename T, std::size_t n>
    int queryarr_required(  std::string name, std::array<T,n> & value)
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr(name,value);
    }

    template<std::size_t n>
    int queryarr_required(  std::string name, std::array<Set::Scalar,n> & value, Unit unit = Unit::Less())
    {
        if (!contains(name.c_str()))
        {
            Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr(name,value,unit);
    }

    int queryarr_default(   std::string name, std::vector<std::string> & value, std::vector<std::string> defaultvalue)
    {
        if (!contains(name.c_str()))
        {
            addarr(name.c_str(),defaultvalue);
        }
        return amrex::ParmParse::queryarr(name.c_str(),value);
    }

    int queryarr_default(   std::string name, Set::Vector & value, std::string defaultvalue, Unit unit)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                this->addarr(name.c_str(),Util::String::Split(defaultvalue));
            }
            return queryarr(name.c_str(),value,unit);
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr_default(   std::string name, Set::Vector & value, Set::Vector defaultvalue)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                std::vector<Set::Scalar> def_data(AMREX_SPACEDIM);
                for (unsigned int i = 0; i < def_data.size(); i++)
                    def_data[i] = defaultvalue[i];
                
                add(name.c_str(),Util::String::Join(def_data));
                value = defaultvalue;
                return 0;
            }
            return queryarr(name.c_str(),value);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    
    int queryarr_default(   std::string name, Set::Matrix & value, Set::Matrix def)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                std::vector<Set::Scalar> def_data(AMREX_SPACEDIM*AMREX_SPACEDIM);
                for (unsigned int i = 0; i < def_data.size(); i++)
                    def_data[i] = def.data()[i];

                add(name.c_str(),Util::String::Join(def_data));
                value = def;
                return 0;
            }
            return queryarr(name.c_str(),value);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryarr_default(   std::string name, std::vector<double> & value, std::vector<double> defaultvalue)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                addarr(name.c_str(),defaultvalue);
            }
            return queryarr(name.c_str(),value);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryarr_default(   std::string name, std::vector<double> & value, std::vector<std::string> defaultvalue, Unit unit)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                addarr(name.c_str(),defaultvalue);
            }
            return queryarr(name,value, unit);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    template<typename T, std::size_t n>
    int queryarr_default(   std::string name, std::array<T,n> & value, std::array<T,n> defaultvalue)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                std::vector<T> def_data(defaultvalue.begin(), defaultvalue.end());
                addarr(name.c_str(),def_data);
            }
            return queryarr(name.c_str(),value);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<std::size_t n>
    int queryarr_default(   std::string name, std::array<Set::Scalar,n> & value,
                            std::array<std::string,n> defaultvalue, Unit unit)
    {
        try
        {
            if (!contains(name.c_str()))
            {
                std::vector<std::string> def_data(defaultvalue.begin(), defaultvalue.end());
                addarr(name.c_str(),def_data);
            }
            return queryarr(name,value, unit);
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    template <typename T>
    int
    queryclass_enumerate(std::string a_name, std::vector<T> &value, int number = 1)
    {
        value.clear();

        //
        // If only one is present with no subscript, then read
        // it only and return.
        //
        std::string name = a_name;
        if (this->contains(name))
        {
            for (int n = 0; n < number; n++)
            {
                T tmp;
                this->queryclass<T>(name, tmp);
                value.push_back(tmp);
            }
            return 0;
        }

        //
        // Some logic to determine whether we are starting with zero
        // (model0, model1, model2, ...)
        // or one
        // (model1, model2, model3, ...)
        // since both are supported
        //
        int start = -1;
        std::string name0 = a_name + std::to_string(0);
        std::string name1 = a_name + std::to_string(1);
        if (this->contains(name0.c_str()))
        {
            start = 0;
            name = name0;
        }
        else if (this->contains(name1.c_str()))
        {
            start = 1;
            name = name1;
        }
        else
        {
            Util::ParmParseException(INFO,full(name), "Enumerations must begin with 0 or 1");
        }

        //
        // Iterate over items called (model0), model1, model2, etc
        // until no more are found then exit. 
        //
        for (int cntr = start; this->contains(name.c_str()); cntr++)
        {
            if (this->contains(name.c_str()))
            {
                T tmp;
                this->queryclass<T>(name, tmp);
                value.push_back(tmp);
            }
            name = a_name + std::to_string(cntr+1);
        }

        return 0;
    }


    template <typename T>
    int
    query_enumerate(std::string a_name, std::vector<T> &value, int number = 1)
    {
        value.clear();

        //
        // If only one is present with no subscript, then read
        // it only and return.
        //
        std::string name = a_name;
        if (this->contains(name))
        {
            for (int n = 0; n < number; n++)
            {
                T tmp;
                this->query_required(name, tmp);
                value.push_back(tmp);
            }
            return 0;
        }

        //
        // Some logic to determine whether we are starting with zero
        // (model0, model1, model2, ...)
        // or one
        // (model1, model2, model3, ...)
        // since both are supported
        //
        int start = -1;
        std::string name0 = a_name + std::to_string(0);
        std::string name1 = a_name + std::to_string(1);
        if (this->contains(name0.c_str()))
        {
            start = 0;
            name = name0;
        }
        else if (this->contains(name1.c_str()))
        {
            start = 1;
            name = name1;
        }
        else
        {
            Util::ParmParseException(INFO,full(name), "Enumerations must begin with 0 or 1");
        }

        //
        // Iterate over items called (model0), model1, model2, etc
        // until no more are found then exit. 
        //
        for (int cntr = start; this->contains(name.c_str()); cntr++)
        {
            if (this->contains(name.c_str()))
            {
                T tmp;
                this->query_required(name, tmp);
                value.push_back(tmp);
            }
            name = a_name + std::to_string(cntr+1);
        }

        return 0;
    }



    int AnyUnusedInputs(bool inscopeonly = true, bool verbose = false)
    {
        int cnt = 0;
        for (auto li = m_table->begin(), End = m_table->end(); li != End; ++li)
        {
                if (!li->second.m_count)
                {
                    if (inscopeonly && getPrefix() != "")
                    {
                        if (li->first.rfind(getPrefix()+".",0) != std::string::npos)
                        {
                            if (verbose) Util::Warning(INFO,li->first);
                            cnt++;
                        }
                    }
                    else
                    {
                        if (verbose) Util::Warning(INFO,li->first);
                        cnt++;
                    }
                }
        }
        return cnt;
    }

    std::vector<std::string> GetUnusedInputs()
    {
        std::vector<std::string> ret;
        for (auto li = m_table->begin(), End = m_table->end(); li != End; ++li)
        {
            if (!li->second.m_count && li->first.rfind(getPrefix()+".",0) != std::string::npos)
            {
                ret.push_back(li->first);
            }
        }
        return ret;
    }

    static int AllUnusedInputs()
    {
        ParmParse pp;
        int cnt = 0;
        for (auto li = pp.m_table->begin(), End = pp.m_table->end(); li != End; ++li)
        {
            if (!li->second.m_count)
            {
                Util::Warning(INFO,li->first);
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


    template<class T>
    void queryclass(std::string name, T * value)
    {
        auto old_prefix = m_prefix;
        try
        {
            if (old_prefix.empty()) m_prefix = name;
            else m_prefix.append(".").append(name);
            T::Parse(*value, *this);
            std::vector<std::string> unused_inputs = GetUnusedInputs();
            if (unused_inputs.size())
            {
                std::stringstream ss;
                for (unsigned int i=0; i < unused_inputs.size(); i++)
                    ss << "\n\t" << unused_inputs[i];
                Util::ParmParseException(INFO,name,"The following inputs were specified but not used",ss.str());
            }
        }
        catch (...)
        {
            m_prefix = old_prefix;
            Util::ParmParseException(INFO,full(name));
        }
        m_prefix = old_prefix;
    }
    template<class T>
    void queryclass(std::string name, T & value)
    {
        auto old_prefix = m_prefix;
        try
        {
            if (old_prefix.empty()) m_prefix = name;
            else m_prefix.append(".").append(name);
            T::Parse(value, *this);
            std::vector<std::string> unused_inputs = GetUnusedInputs();
            if (unused_inputs.size())
            {
                std::stringstream ss;
                for (unsigned int i=0; i < unused_inputs.size(); i++)
                    ss << "\n\t" << unused_inputs[i];
                Util::ParmParseException(INFO,name,"The following inputs were specified but not used",ss.str());
            }
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            m_prefix = old_prefix;
            Util::ParmParseException(INFO,full(name));
        }
        m_prefix = old_prefix;
    }

    template<class T>
    void queryclass(T * value)
    {
        try
        {
            T::Parse(*value, *this);
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,getPrefix(), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,getPrefix());
        }
    }
    template<class T>
    void queryclass(T & value)
    {
        try
        {
            T::Parse(value, *this);
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,getPrefix(), e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO, getPrefix());
        }
    }

    //
    // Variadic template parsing operator to assign a pointer to
    // one of a set of possible class objects, then call that method's
    // Parse function.
    //
    // If there is more than one instantiating class, then type must be set.
    // Otherwise, no type is necessary.
    //
    template<typename... IC, typename... Args, typename PTRTYPE>
    void select (std::string name, PTRTYPE *& ic_eta, Args&&... args)
    {
        try
        {
            // if there is only one IC arg provided, we don't need to check the type - assume that
            // it is the default.
            if constexpr (sizeof...(IC) == 0)
            {
                using first_IC = std::tuple_element_t<0, std::tuple<IC...>>;
                ic_eta = new first_IC(std::forward<Args>(args)..., (*this), name + "." + std::string(first_IC::name));
            }
            // otherwise, check the type.
            else
            {
                std::string type = "";
                this->query_required(name + ".type", type);
                bool matched = ((   type == IC::name 
                                    ? (ic_eta = new IC(std::forward<Args>(args)..., (*this), name + "." + std::string(IC::name))),
                                    true
                                    : false) || ...);
                if (!matched)
                    Util::ParmParseException(INFO,full(name), type, " not a valid type for ", name);
            }
        }
        catch (std::runtime_error &e)
        {
            Util::ParmParseException(INFO,getPrefix(),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,getPrefix());
        }
    }

    //
    // Identical to the above, except sets the type automatically to the
    // firs specified 
    //
    template<typename FirstIC, typename... IC, typename... Args, typename PTRTYPE>
    void select_default (std::string name, PTRTYPE *& ic_eta, Args&&... args)
    {
        std::string type = "";

        this->query_default(name + ".type", type, FirstIC::name);

        bool matched =
            ((  type == FirstIC::name 
                ? (ic_eta = new FirstIC(std::forward<Args>(args)..., (*this), name + "." + std::string(FirstIC::name))),
                true : false))
            || 
            ((  type == IC::name 
                ? (ic_eta = new IC(std::forward<Args>(args)..., (*this), name + "." + std::string(IC::name))),
                true : false) || ...);


        if (!matched)
            Util::ParmParseException(INFO,full(name), type," not a valid type for ",name);
    }


    //
    //
    template<typename... IC, typename... Args, typename PTRTYPE>
    void select_enumerate (std::string a_name, std::vector<PTRTYPE*> & value, Args&&... args)
    {

        value.clear();

        //
        // If only one is present with no subscript, then read
        // it only and return.
        //
        std::string name = a_name;
        if (this->contains(name))
        {
            PTRTYPE *tmp;
            this->select<IC...>(a_name, tmp, args...);
            value.push_back(tmp);
            return;
        }

        //
        // Some logic to determine whether we are starting with zero
        // (model0, model1, model2, ...)
        // or one
        // (model1, model2, model3, ...)
        // since both are supported
        //
        int start = -1;
        std::string name0 = a_name + std::to_string(0);
        std::string name1 = a_name + std::to_string(1);
        if (this->contains(name0.c_str()))
        {
            start = 0;
            name = name0;
        }
        else if (this->contains(name1.c_str()))
        {
            start = 1;
            name = name1;
        }
        else
        {
            Util::ParmParseException(INFO,full(name), "Enumerations must begin with 0 or 1");
        }

        //
        // Iterate over items called (model0), model1, model2, etc
        // until no more are found then exit. 
        //
        for (int cntr = start; this->contains(name.c_str()); cntr++)
        {
            if (this->contains(name.c_str()))
            {
                PTRTYPE *tmp;
                this->select<IC...>(name, tmp, args...);
                value.push_back(tmp);
            }
            name = a_name + std::to_string(cntr+1);
        }
    }



    //
    // Similar to select but specialized for main functions
    //
    template<typename... INTEGRATOR, typename... Args, typename PTRTYPE>
    void select_main (PTRTYPE *& ic_eta, Args&&... args)
    {
        std::string type = "";

        this->query_required("alamo.program", type);
        
        bool matched = ((type == INTEGRATOR::name 
                        ? (ic_eta = new INTEGRATOR(std::forward<Args>(args)..., (*this))),
                        true
                        : false) || ...);
        if (!matched)
            Util::ParmParseException(INFO,getPrefix(),type," not a valid type for ",type);
    }

    //
    // Similar to select_main but works for one function only
    // and doesn't require a type specifier.
    //

    // with variadic arguments
    template<typename INTEGRATOR, typename Args, typename PTRTYPE>
    void select_only (PTRTYPE *& ic_eta, Args&& args)
    {
        ic_eta = new INTEGRATOR(std::forward<Args>(args), (*this));
    }
    // without variadic arguments
    template<typename INTEGRATOR, typename PTRTYPE>
    void select_only (PTRTYPE *& ic_eta)
    {
        ic_eta = new INTEGRATOR((*this));
    }

    
    //
    // functionally simlar to other kinds of select, execpt works on
    // template-based static dispatch.
    //

    template <typename... OBJ, typename CLASS, typename... Args>
    void select(std::string name, CLASS& value, Args&&... args)
    {
        pushPrefix(name);
        static_polymorphism_parser<CLASS, 0, OBJ...>(value, std::forward<Args>(args)...);
        popPrefix();
    }

private:
    //
    // This function is a virtual "swtich / case" block for static
    // dispatch. It works on a faux-virtual class "CLASS" that must 
    // contain a tuple of "obj" with each "obj" being a derived type.
    // Recursion is used to unroll the list, identify the correct class
    // based on "name/type", set the "selected" index, then call the
    // appropriate class' Parse function.
    // 
    // This is only for use by the "select" functions.
    // 
    template <typename CLASS, int N, typename... OBJ, typename... Args>
    void static_polymorphism_parser(CLASS& value, Args&&... args)
    {
        if constexpr (N == 0)
        {
            std::string type;
            query_default("type",type,value.names[0]);

            for (unsigned int i = 0; i < sizeof...(OBJ); i++)
                if (type == value.names[i])
                    value.selected = i;

            if (value.selected < 0)
                Util::ParmParseException(INFO,getPrefix(),
                                "Error reading " + getPrefix() + 
                                ", invalid type " + type);

            
            pushPrefix(type);
        }

        if constexpr (N < sizeof...(OBJ))
        {
            if (value.selected == N)
            {
                std::get<N>(value.obj).Parse(
                    std::get<N>(value.obj), *this, std::forward<Args>(args)...);
                popPrefix();
                return;
            }
            else 
                return static_polymorphism_parser<CLASS, N+1, OBJ...>(
                    value, std::forward<Args>(args)...);
        }
        else Util::Abort(INFO);
    }

public:

    template <int N>
    void query_exactly( std::vector<std::string> names, std::pair<std::string, Set::Scalar> values[N],
                        std::vector<Unit> units = std::vector<Unit>())
    {
        try
        {
            Util::AssertException(  INFO,TEST(units.size() == 0 || units.size() == names.size()),
                                    "# of units must be 0, 1, or ", names.size(), " but got ", units.size());

            int cnt = 0;
            std::vector<std::string> read;
            for (unsigned int n = 0; n < names.size(); n++)
            {
                if (amrex::ParmParse::contains(names[n].c_str()))
                {
                    read.push_back(names[n]);
                    cnt++;
                }
            }
            Util::AssertException(  INFO, TEST(cnt == N),
                                    " Incorrect number of values specified: only ", N,
                                    " values are allowed, but received ",
                                    Util::String::Join(read,", "));


            cnt = 0;
            for (unsigned int n = 0; n < names.size(); n++)
            {
                if (amrex::ParmParse::contains(names[n].c_str()))
                {
                    values[cnt].first = names[n];

                    if (units.size() == 0)
                        query_required(names[n],values[cnt].second);
                    else if (units.size() == 1)
                        query_required(names[n],values[cnt].second,units[0]);
                    else
                        query_required(names[n],values[cnt].second,units[n]);
                    cnt++;
                }
            }
        }
        catch(...)
        {
            std::string names_str = "[";
            for (unsigned int i = 0; i < names.size(); i++)
            {
                names_str += full(names[i]);
                if (i < names.size()-1) names_str += ", ";
            }
            names_str += "]";
            Util::ParmParseException(INFO,names_str);
        }
    }
};
}
#endif
