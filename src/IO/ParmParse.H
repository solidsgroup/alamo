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
//    | :bdg-success:`query_switch`     | Validate a string against allowed values and run the branch     |
//    |                                 | lambda associated with the selected value.                     |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_if`         | Read a Boolean value and run the supplied lambda when enabled.  |
//    +---------------------------------+-----------------------------------------------------------------+
//    | :bdg-success:`query_if_else`    | Read a Boolean value and run the corresponding supplied lambda. |
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
#include <functional>
#include <iostream>
#include <list>
#include <optional>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "Util/Util.H"
#include "Unit/Unit.H"
#include "Set/Base.H"
#include "IO/InputScraper.H"
#include "AMReX_ParmParse.H"


#define pp_query_required(...) pp.query_required(__VA_ARGS__)
#define pp_query_default(...) pp.query_default(__VA_ARGS__)
#define pp_query_validate(...) pp.query_validate(__VA_ARGS__)
#define pp_query_switch(...) pp.query_switch(__VA_ARGS__)
#define pp_query_if(...) pp.query_if(__VA_ARGS__)
#define pp_query_if_else(...) pp.query_if_else(__VA_ARGS__)
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
public:
    using InputNode = InputScraper::InputNode;
    template<typename... Args>
    struct ForwardArgs
    {
        std::tuple<Args...> values;

        explicit ForwardArgs(Args... args) : values(std::forward<Args>(args)...) {}
    };

    template<typename Arg>
    using ForwardArgStorage = std::conditional_t<
        std::is_lvalue_reference_v<Arg>, Arg, std::remove_cvref_t<Arg>>;

    template<typename... Args>
    static ForwardArgs<ForwardArgStorage<Args>...> forward_args(Args&&... args)
    {
        return ForwardArgs<ForwardArgStorage<Args>...>(std::forward<Args>(args)...);
    }

private:
    friend class InputScraper;
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
    void RecordInput(   std::string name,
                        std::string directive,
                        const std::source_location &location,
                        std::vector<std::string> options = {},
                        std::optional<std::string> default_value = std::nullopt,
                        bool has_unnamed_default = false);
    void RecordConstraint(  std::string kind,
                            int count,
                            std::vector<std::string> members,
                            std::vector<std::string> units,
                            const std::source_location &location);
    void PrintTraversalBranch(std::string name, const std::string &value);
    void RecordInputAndContinue(
        std::string name,
        const std::source_location &location,
        std::vector<std::string> options = {},
        std::optional<std::string> default_value = std::nullopt,
        bool has_unnamed_default = false,
        const std::source_location &handler_location = std::source_location::current())
    {
        if (InTraversalMode())
            RecordInput(name, InputScraper::DirectiveName(handler_location), location,
                        std::move(options), std::move(default_value),
                        has_unnamed_default);
    }
    using TraversalBranchScope = InputScraper::TraversalBranchScope;

    template<typename T, typename... Args>
    static T *ConstructSelected(ForwardArgs<Args...> &args, ParmParse &pp, const std::string &name)
    {
        return std::apply(
            [&](auto &...constructor_args) -> T *
            {
                return new T(constructor_args..., pp, name);
            },
            args.values);
    }

    template<typename T, typename = void>
    struct HasSelectName : std::false_type {};

    template<typename T>
    struct HasSelectName<T, std::void_t<decltype(T::name)>> : std::true_type {};

    template<typename T, typename... Args>
    static constexpr bool HasSelectedConstructor =
        std::is_constructible_v<T, std::add_lvalue_reference_t<Args>...,
                                ParmParse &, std::string>;

    template<typename T, typename... Args>
    static T *ConstructSelectDefault(ForwardArgs<Args...> &args)
    {
        return std::apply(
            [&](auto &...constructor_args) -> T *
            {
                static_assert(  std::is_constructible_v<T, decltype(constructor_args)...>,
                                "The unnamed select default cannot be constructed from forward_args");
                return new T(constructor_args...);
            },
            args.values);
    }

    template<typename T, typename... Args>
    static void AppendSelectName(std::vector<std::string> &names, ForwardArgs<Args...> &)
    {
        if constexpr (HasSelectName<T>::value && HasSelectedConstructor<T, Args...>)
            names.emplace_back(T::name);
    }

    template<typename T, typename... Args, typename PTRTYPE>
    static bool TrySelectNamed( const std::string &type, const std::string &name,
                                PTRTYPE *&value, ForwardArgs<Args...> &args,
                                ParmParse &pp)
    {
        if constexpr (HasSelectName<T>::value && HasSelectedConstructor<T, Args...>)
        {
            if (type == T::name)
            {
                value = ConstructSelected<T>(args, pp, name + "." + std::string(T::name));
                return true;
            }
        }
        return false;
    }
    
public:
    static void SetTraversalMode(bool enabled);
    static bool InTraversalMode();
    static bool ShouldExecute();
    static const InputNode &InputTree();
    static void ClearInputTree();
    static void SetTraversalOutputFile(std::string path);
    static const std::string &TraversalOutputFile();
    static void WriteInputTreeJson(std::ostream &os);
    static void WriteInputTreeJsonFile(const std::string &path);
    static bool IgnoreInTraversalMode(
        std::string note = "This input parser is not yet compliant with traversal mode; its inputs are not fully documented.",
        const std::source_location &location = std::source_location::current());

    ParmParse(std::string arg) : amrex::ParmParse::ParmParse(arg) {Define();} ;
    ParmParse() : amrex::ParmParse::ParmParse() {Define();} ;
    std::string getPrefix() const {return m_prefix;};
    void ignore(std::string name,
                const std::source_location & location = std::source_location::current())
    {
        (void)location;
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


    void forbid(std::string name, std::string explanation,
                const std::source_location & location = std::source_location::current())
    {
        (void)location;
        if (amrex::ParmParse::contains(full(name).c_str()))
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),full(name)," forbidden: ", explanation);
        }
        std::set<std::string> subs = amrex::ParmParse::getEntries(full(name));
        if (subs.size())
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),full(name)," forbidden: ", explanation);
        }
    }

    bool contains(  std::string name,
                    const std::source_location & location = std::source_location::current())
    {
        (void)location;
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

    template<typename T>
    int query(  std::string name, T & value,
                const std::source_location & location = std::source_location::current())
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, bool>)
            RecordInputAndContinue(name, location, {"0", "1"});
        else
            RecordInputAndContinue(name, location);
        if (InTraversalMode() && !contains(name.c_str())) return 0;
        return amrex::ParmParse::query(name.c_str(), value);
    }

    int queryunit ( std::string name, Unit &value,
                    const std::source_location & location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode() && !contains(name.c_str())) return 0;
        try
        {
            std::string strvalue;
            int retval = amrex::ParmParse::query(name.c_str(),strvalue);
            value = Unit::Parse(strvalue);
            return retval;
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryunit ( std::string name, Unit &value, const Unit type,
                    const std::source_location & location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode() && !contains(name.c_str())) return 0;
        try
        {
            int retval = queryunit(name,value, location);
            if (!value.isType(type) && !value.isType(Unit::Less()))
            {
                if (!InTraversalMode())
                    Util::ParmParseException(   INFO, full(name),
                                                "value requiested had wrong units:", value.normalized_unitstring(),
                                                ", units requested are of SI type ",type.normalized_unitstring());
            }
            return retval;
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryunit ( std::string name, Set::Scalar &value, const Unit type,
                    const std::source_location & location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode() && !contains(name.c_str())) return 0;
        try
        {
            Unit read;
            int retval = queryunit(name,read,type, location);
            value = read.normalized_value();
            return retval;
        }
        catch (std::runtime_error &e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name), e.what());
        }
        return -1;
    }

    template<typename T> 
    int query_required( std::string name, T & value,
                        const std::source_location & location = std::source_location::current())
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, bool>)
            RecordInputAndContinue(name, location, {"0", "1"});
        else
            RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (!contains(name.c_str())) return 0;
            return amrex::ParmParse::query(name.c_str(), value);
        }
        try
        {
            if (!contains(name.c_str()))
            {
                if (!InTraversalMode())
                    Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
            }
            return query(name.c_str(),value);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T> 
    int query_required( std::string name, T & value, const Unit type,
                        const std::source_location & location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (!contains(name.c_str())) return 0;
            std::string strvalue;
            int retval = amrex::ParmParse::query(name.c_str(), strvalue);
            Unit read = Unit::Parse(strvalue);
            value = read.normalized_value();
            return retval;
        }
        try
        {
            if (!contains(name.c_str()))
            {
                if (!InTraversalMode())
                    Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
            }
            return queryunit(name.c_str(), value, type, location);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    template<typename T>
    int query_default(  std::string name, T & value, T defaultvalue,
                        const std::source_location & location = std::source_location::current())
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, bool>)
            RecordInputAndContinue( name, location, {"0", "1"},
                                    InputScraper::InputValueString(defaultvalue));
        else
            RecordInputAndContinue( name, location, {},
                                    InputScraper::InputValueString(defaultvalue));
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::query(name.c_str(), value);
            value = defaultvalue;
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    } 

    template<typename T>
    int query_default(  std::string name, T & value, std::string defaultvalue, const Unit type,
                        const std::source_location & location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, defaultvalue);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::string strvalue;
                int retval = amrex::ParmParse::query(name.c_str(), strvalue);
                Unit read = Unit::Parse(strvalue);
                value = read.normalized_value();
                return retval;
            }
            Unit read = Unit::Parse(defaultvalue);
            value = read.normalized_value();
            return 0;
        }
        try
        {
            if (!contains(name.c_str()))
            {
                add(name.c_str(),defaultvalue);
            }
            return queryunit(name.c_str(),value,type, location);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    } 

    int query_validate( std::string name, int & value, std::vector<int> possibleintvals,
                        const std::source_location & location = std::source_location::current())
    {
        std::vector<std::string> options;
        for (const auto &possibleval : possibleintvals)
            options.push_back(std::to_string(possibleval));
        RecordInputAndContinue(name, location, options);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::query(name.c_str(), value);
            value = possibleintvals[0];
            return 0;
        }

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

            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),"' expected [", ss.str(), "] but got ", value);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals, bool firstbydefault,
                        const std::source_location &location = std::source_location::current())
    {        
        std::vector<std::string> options;
        for (const auto &possibleval : possiblecharvals)
            options.push_back(std::string(possibleval));
        RecordInputAndContinue(name, location, options);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::query(name.c_str(), value);
            value = std::string(possiblecharvals[0]);
            return 0;
        }

        try
        {
            // if not using default values, then the input must be specified
            if (!firstbydefault)
            {
                if (!amrex::ParmParse::contains(name.c_str()))
                {
                    if (!InTraversalMode())
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }

        return -1;
    }

    int query_validate( std::string name, std::string & value, std::vector<const char *> possiblecharvals,
                        const std::source_location &location = std::source_location::current())
    {
        try
        {
            return query_validate(name,value,possiblecharvals,true, location);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int query_switch(   std::string name,
                        std::initializer_list<std::pair<std::string, std::function<void()>>> cases,
                        const std::source_location &location = std::source_location::current())
    {
        std::string value;
        std::vector<std::string> possiblevals;
        possiblevals.reserve(cases.size());
        for (const auto &entry : cases)
            possiblevals.push_back(entry.first);
        {
            RecordInputAndContinue(name, location, possiblevals);
        }

        auto set_switch_value = [&](const std::string &switch_value)
        {
            this->remove(name.c_str());
            this->add(name.c_str(), switch_value);
        };

        auto print_switch_branch = [&](const std::string &switch_value)
        {
            PrintTraversalBranch(name, switch_value);
        };

        auto run_case = [&](const std::string &switch_value, auto &&action)
        {
            set_switch_value(switch_value);
            print_switch_branch(switch_value);
            TraversalBranchScope branch(*this, name, switch_value);
            action();
        };

        if (InTraversalMode())
        {
            std::string original_value;
            bool had_original_value = amrex::ParmParse::contains(name.c_str());
            if (had_original_value)
                amrex::ParmParse::query(name.c_str(), original_value);

            for (const auto &entry : cases)
                run_case(entry.first, entry.second);

            if (had_original_value)
                set_switch_value(original_value);
            else
                this->remove(name.c_str());
            return 0;
        }

        int retval = -1;
        try
        {
            if (!amrex::ParmParse::contains(name.c_str()))
            {
                add(name.c_str(), possiblevals[0]);
            }

            retval = amrex::ParmParse::query(name.c_str(),value);

            bool ok = false;
            for (const auto & possibleval : possiblevals)
            {
                if (value == possibleval) ok = true;
            }

            if (!ok)
            {
                std::stringstream ss;
                ss << possiblevals[0];
                for (unsigned int i = 1; i < possiblevals.size(); i++)
                    ss << "," << possiblevals[i];

                Util::ParmParseException(INFO,full(name),"' expected [", ss.str(), "] but got ", value);
            }
        }
        catch (std::runtime_error & e)
        {
            Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            Util::ParmParseException(INFO,full(name));
        }

        bool matched = false;
        for (const auto &entry : cases)
        {
            if (!matched && value == entry.first)
            {
                entry.second();
                matched = true;
            }
        }

        if (!matched)
        {
            Util::ParmParseException(INFO,full(name),value," not a valid value for ",full(name));
        }

        return retval;
    }

    template <typename Action>
    int query_if(   std::string name, Action &&action,
                    const std::source_location &location = std::source_location::current())
    {
        static_assert(  std::is_invocable_v<Action &>,
                        "query_if action must be callable without arguments");

        RecordInputAndContinue(name, location, {"0", "1"}, "0");

        auto set_value = [&](const std::string &value)
        {
            this->remove(name.c_str());
            this->add(name.c_str(), value);
        };

        if (InTraversalMode())
        {
            std::string original_value;
            bool had_original_value = amrex::ParmParse::contains(name.c_str());
            if (had_original_value)
                amrex::ParmParse::query(name.c_str(), original_value);

            set_value("1");
            PrintTraversalBranch(name, "1");
            {
                TraversalBranchScope branch(*this, name, "1");
                action();
            }

            if (had_original_value)
                set_value(original_value);
            else
                this->remove(name.c_str());
            return 0;
        }

        int enabled = 0;
        int retval = amrex::ParmParse::query(name.c_str(), enabled);
        if (retval && enabled != 0 && enabled != 1)
            Util::ParmParseException(INFO, full(name), "expected 0 or 1 but got ", enabled);
        if (enabled)
            action();
        return retval;
    }

    template <typename TrueAction, typename FalseAction>
    int query_if_else(  std::string name,
                        TrueAction &&true_action,
                        FalseAction &&false_action,
                        const std::source_location &location = std::source_location::current())
    {
        static_assert(  std::is_invocable_v<TrueAction &>,
                        "query_if_else true action must be callable without arguments");
        static_assert(  std::is_invocable_v<FalseAction &>,
                        "query_if_else false action must be callable without arguments");

        RecordInputAndContinue(name, location, {"0", "1"}, "0");

        auto set_value = [&](const std::string &value)
        {
            this->remove(name.c_str());
            this->add(name.c_str(), value);
        };

        auto run_action = [&](const std::string &value, auto &&action)
        {
            set_value(value);
            PrintTraversalBranch(name, value);
            TraversalBranchScope branch(*this, name, value);
            action();
        };

        if (InTraversalMode())
        {
            std::string original_value;
            bool had_original_value = amrex::ParmParse::contains(name.c_str());
            if (had_original_value)
                amrex::ParmParse::query(name.c_str(), original_value);

            run_action("0", false_action);
            run_action("1", true_action);

            if (had_original_value)
                set_value(original_value);
            else
                this->remove(name.c_str());
            return 0;
        }

        int enabled = 0;
        int retval = amrex::ParmParse::query(name.c_str(), enabled);
        if (retval && enabled != 0 && enabled != 1)
            Util::ParmParseException(INFO, full(name), "expected 0 or 1 but got ", enabled);
        if (enabled)
            true_action();
        else
            false_action();
        return retval;
    }

    // special case for strings
    int query_default(  std::string name, std::string & value, const char *defaultvalue,
                        const std::source_location &location = std::source_location::current())
    {
        try
        {
            return query_default(name, value, std::string(defaultvalue), location);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    // special case for bools
    int query_default(  std::string name, int & value, bool defaultvalue,
                        const std::source_location &location = std::source_location::current())
    {
        try
        {
            int defaultint = 0;
            if (defaultvalue) defaultint = 1;
            return query_default(name, value, defaultint, location);
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    // validate filenames
    int query_file( std::string name, std::string & value, bool copyfile, bool checkfile,
                    const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (!contains(name.c_str())) return 0;
            return amrex::ParmParse::query(name.c_str(), value);
        }
        try
        {
            if (!contains(name.c_str()))
            {
                if (!InTraversalMode())
                    Util::ParmParseException(INFO,full(name),full(name)," must be specified");
            }

            int retval = query(name.c_str(),value);

            if (ShouldExecute() && amrex::ParallelDescriptor::IOProcessor())
            {
                if ( checkfile && ! std::filesystem::exists(value))
                {
                    if (!InTraversalMode())
                        Util::ParmParseException(INFO,full(name),full(name)," does not exist");
                }
                if ( checkfile && !std::filesystem::is_regular_file(value))
                {
                    if (!InTraversalMode())
                        Util::ParmParseException(INFO,full(name),full(name)," is not a regular file");
                }
                if ( copyfile )
                {
                    if (!InTraversalMode())
                        Util::CopyFileToOutputDir(value, true, full(name));
                }
            }
            return retval;
        }
        catch (std::runtime_error & e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int query_file( std::string name, std::string & value, bool copyfile,
                    const std::source_location & location = std::source_location::current())
    {
        return query_file(name,value,copyfile,true, location);
    }
    int query_file( std::string name, std::string & value,std::source_location loc = std::source_location::current())
    {
        return query_file(name,value,true,true, loc);
    }


    template<typename T>
    int queryarr(   std::string name, std::vector<T> & value,
                    const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::queryarr(name.c_str(), value);
            value.clear();
            return 0;
        }
        try
        {
            return amrex::ParmParse::queryarr(name.c_str(),value);
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, std::vector<Set::Scalar> & value, Unit unit = Unit::Less(),
                    const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                if (unit.isType(Unit::Less()))
                    return amrex::ParmParse::queryarr(name.c_str(), value);

                value.clear();
                std::vector<std::string> valstrings;
                int retval = amrex::ParmParse::queryarr(name.c_str(), valstrings);
                for (unsigned int i = 0; i < valstrings.size(); i++)
                {
                    Unit unitvalue = Unit::Parse(valstrings[i]);
                    value.push_back(unitvalue.normalized_value());
                }
                return retval;
            }
            value.clear();
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, Set::Vector & value, Unit unit = Unit::Less(),
                    const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::vector<Set::Scalar> vals;
                if (unit.isType(Unit::Less()))
                    amrex::ParmParse::queryarr(name.c_str(), vals);
                else
                {
                    std::vector<std::string> valstrings;
                    amrex::ParmParse::queryarr(name.c_str(), valstrings);
                    for (unsigned int i = 0; i < valstrings.size(); i++)
                    {
                        Unit unitvalue = Unit::Parse(valstrings[i]);
                        vals.push_back(unitvalue.normalized_value());
                    }
                }
                value = Set::Vector::Zero();
                for (int i = 0; i < AMREX_SPACEDIM && i < static_cast<int>(vals.size()); i++)
                    value(i) = vals[i];
                return 0;
            }
            value = Set::Vector::Zero();
            return 0;
        }
        try
        {
            std::vector<Set::Scalar> vals;
            queryarr(name.c_str(), vals, unit);
            if (vals.size() < AMREX_SPACEDIM) 
            {
                if (!InTraversalMode())
                    Util::ParmParseException(   INFO,full(name),full(name),
                                                " requires at least ", AMREX_SPACEDIM,
                                                " arguments, got ",vals.size());
            }
            for (int i = 0; i < AMREX_SPACEDIM; i++) value(i) = vals[i];
            return 0;
        }
        catch(...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr(   std::string name, Set::Matrix & value, Unit unit = Unit::Less(),
                    const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::vector<Set::Scalar> vals;
                if (unit.isType(Unit::Less()))
                    amrex::ParmParse::queryarr(name.c_str(), vals);
                else
                {
                    std::vector<std::string> valstrings;
                    amrex::ParmParse::queryarr(name.c_str(), valstrings);
                    for (unsigned int i = 0; i < valstrings.size(); i++)
                    {
                        Unit unitvalue = Unit::Parse(valstrings[i]);
                        vals.push_back(unitvalue.normalized_value());
                    }
                }
                value = Set::Matrix::Zero();
                for (int i = 0; i < AMREX_SPACEDIM * AMREX_SPACEDIM && i < static_cast<int>(vals.size()); i++)
                    value.data()[i] = vals[i];
                return 0;
            }
            value = Set::Matrix::Zero();
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),full(name)," needs either 4 or 9 components, but got ",vals.size());
        }    
        return 0;
    }

    template<typename T>
    int queryarr_required(  std::string name, std::vector<T> & value,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::queryarr(name.c_str(), value);
            value.clear();
            return 0;
        }
        if (!contains(name.c_str()))
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr<T>(name,value);
    }

    int queryarr_required(  std::string name, std::vector<Set::Scalar> & value, Unit unit = Unit::Less(),
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                if (unit.isType(Unit::Less()))
                    return amrex::ParmParse::queryarr(name.c_str(), value);

                value.clear();
                std::vector<std::string> valstrings;
                int retval = amrex::ParmParse::queryarr(name.c_str(), valstrings);
                for (unsigned int i = 0; i < valstrings.size(); i++)
                {
                    Unit unitvalue = Unit::Parse(valstrings[i]);
                    value.push_back(unitvalue.normalized_value());
                }
                return retval;
            }
            value.clear();
            return 0;
        }
        if (!contains(name.c_str()))
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        }        
        return queryarr(name,value,unit);
    }
    int queryarr_required(  std::string name, Set::Vector & value, Unit unit = Unit::Less(),
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location);
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::vector<Set::Scalar> vals;
                if (unit.isType(Unit::Less()))
                    amrex::ParmParse::queryarr(name.c_str(), vals);
                else
                {
                    std::vector<std::string> valstrings;
                    amrex::ParmParse::queryarr(name.c_str(), valstrings);
                    for (const auto &valstring : valstrings)
                        vals.push_back(Unit::Parse(valstring).normalized_value());
                }
                value = Set::Vector::Zero();
                for (int i = 0; i < AMREX_SPACEDIM && i < static_cast<int>(vals.size()); i++)
                    value(i) = vals[i];
                return 0;
            }
            value = Set::Vector::Zero();
            return 0;
        }
        if (!contains(name.c_str()))
            Util::ParmParseException(INFO,full(name),"required value for ",full(name)," missing");
        return queryarr(name,value,unit);
    }
    int queryarr_default(   std::string name, std::vector<std::string> & value, std::vector<std::string> defaultvalue,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, InputScraper::InputArrayString(defaultvalue));
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::queryarr(name.c_str(), value);
            value = defaultvalue;
            return 0;
        }
        if (!contains(name.c_str()))
        {
            addarr(name.c_str(),defaultvalue);
        }
        return amrex::ParmParse::queryarr(name.c_str(),value);
    }

    int queryarr_default(   std::string name, Set::Vector & value, std::string defaultvalue, Unit unit,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, defaultvalue);
        if (InTraversalMode())
        {
            std::vector<std::string> data = Util::String::Split(defaultvalue);
            if (contains(name.c_str()))
                amrex::ParmParse::queryarr(name.c_str(), data);

            value = Set::Vector::Zero();
            for (int i = 0; i < AMREX_SPACEDIM && i < static_cast<int>(data.size()); i++)
            {
                Unit unitvalue = Unit::Parse(data[i]);
                value(i) = unitvalue.normalized_value();
            }
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name), e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    int queryarr_default(   std::string name, Set::Vector & value, Set::Vector defaultvalue,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, InputScraper::InputArrayString(defaultvalue));
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::vector<Set::Scalar> vals;
                amrex::ParmParse::queryarr(name.c_str(), vals);
                value = Set::Vector::Zero();
                for (int i = 0; i < AMREX_SPACEDIM && i < static_cast<int>(vals.size()); i++)
                    value(i) = vals[i];
                return 0;
            }
            value = defaultvalue;
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }
    
    int queryarr_default(   std::string name, Set::Matrix & value, Set::Matrix def,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, InputScraper::InputArrayString(def));
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
            {
                std::vector<Set::Scalar> vals;
                amrex::ParmParse::queryarr(name.c_str(), vals);
                value = Set::Matrix::Zero();
                for (int i = 0; i < AMREX_SPACEDIM * AMREX_SPACEDIM && i < static_cast<int>(vals.size()); i++)
                    value.data()[i] = vals[i];
                return 0;
            }
            value = def;
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryarr_default(   std::string name, std::vector<double> & value, std::vector<double> defaultvalue,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, InputScraper::InputArrayString(defaultvalue));
        if (InTraversalMode())
        {
            if (contains(name.c_str()))
                return amrex::ParmParse::queryarr(name.c_str(), value);
            value = defaultvalue;
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }

    int queryarr_default(   std::string name, std::vector<double> & value, std::vector<std::string> defaultvalue, Unit unit,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(name, location, {}, InputScraper::InputArrayString(defaultvalue));
        if (InTraversalMode())
        {
            std::vector<std::string> data = defaultvalue;
            if (contains(name.c_str()))
                amrex::ParmParse::queryarr(name.c_str(), data);

            value.clear();
            for (const auto &entry : data)
            {
                Unit unitvalue = Unit::Parse(entry);
                value.push_back(unitvalue.normalized_value());
            }
            return 0;
        }
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        return -1;
    }


    template <typename T>
    int
    queryclass_enumerate(   std::string a_name, std::vector<T> &value, int number = 1,
                            const std::source_location &location = std::source_location::current())
    {
        {
            RecordInputAndContinue(a_name, location);
        }
        value.clear();
        if (InTraversalMode())
        {
            const std::string template_name = a_name + "0";
            T tmp;
            this->queryclass<T>(template_name, tmp);
            InputScraper::CaptureSequenceTemplate(*this, a_name, template_name);
            return 0;
        }

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
            if (!InTraversalMode())
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
    query_enumerate(std::string a_name, std::vector<T> &value, int number = 1,
                    const std::source_location &location = std::source_location::current())
    {
        {
            RecordInputAndContinue(a_name, location);
        }
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
            if (!InTraversalMode())
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

    template <typename T>
    int
    queryarr_enumerate( std::string a_name, std::vector<std::vector<T>> &value, int number = 1,
                        const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(a_name, location);
        value.clear();

        std::string name = a_name;
        if (this->contains(name))
        {
            for (int n = 0; n < number; n++)
            {
                std::vector<T> tmp;
                this->queryarr_required(name, tmp);
                value.push_back(std::move(tmp));
            }
            return 0;
        }

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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name), "Enumerations must begin with 0 or 1");
        }

        for (int cntr = start; this->contains(name.c_str()); cntr++)
        {
            std::vector<T> tmp;
            this->queryarr_required(name, tmp);
            value.push_back(std::move(tmp));
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
    void queryclass(std::string name, T * value,
                    const std::source_location &location = std::source_location::current())
    {
        {
            RecordInputAndContinue(name, location);
        }
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
                if (!InTraversalMode())
                    Util::ParmParseException(INFO,name,"The following inputs were specified but not used",ss.str());
            }
        }
        catch (...)
        {
            m_prefix = old_prefix;
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        m_prefix = old_prefix;
    }
    template<class T>
    void queryclass(std::string name, T & value,
                    const std::source_location &location = std::source_location::current())
    {
        {
            RecordInputAndContinue(name, location);
        }
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
                if (!InTraversalMode())
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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,full(name));
        }
        m_prefix = old_prefix;
    }

    template<class T>
    void queryclass(T * value,
                    const std::source_location &location = std::source_location::current())
    {
        (void)location;
        try
        {
            T::Parse(*value, *this);
        }
        catch (std::runtime_error &e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,getPrefix(), e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,getPrefix());
        }
    }
    template<class T>
    void queryclass(T & value,
                    const std::source_location &location = std::source_location::current())
    {
        (void)location;
        try
        {
            T::Parse(value, *this);
        }
        catch (std::runtime_error &e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,getPrefix(), e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
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
    template<typename... IC, typename PTRTYPE>
    void select (   std::string name, PTRTYPE *& ic_eta,
                    const std::source_location &location = std::source_location::current())
    {
        auto args = forward_args();
        select<IC...>(std::move(name), ic_eta, args, location);
    }

    template<typename... IC, typename... Args, typename PTRTYPE>
    void select (   std::string name, PTRTYPE *& ic_eta, ForwardArgs<Args...> args,
                    const std::source_location &location = std::source_location::current())
    {
        try
        {
            std::string type_name = name + ".type";
            std::vector<std::string> possiblevals = {std::string(IC::name)...};

            if (InTraversalMode())
            {
                {
                    RecordInputAndContinue(type_name, location, possiblevals);
                }

                std::string original_type;
                bool had_original_type = this->contains(type_name.c_str());
                if (had_original_type)
                    amrex::ParmParse::query(type_name.c_str(), original_type);

                auto set_select_value = [&](const std::string &select_value)
                {
                    this->remove(type_name.c_str());
                    this->add(type_name.c_str(), select_value);
                };

                auto restore_select_value = [&]()
                {
                    this->remove(type_name.c_str());
                    if (had_original_type)
                        this->add(type_name.c_str(), original_type);
                };

                auto traverse_one = [&]<typename T>()
                {
                    std::string select_value = T::name;
                    set_select_value(select_value);
                    PrintTraversalBranch(type_name, select_value);
                    TraversalBranchScope branch(*this, type_name, select_value);
                    PTRTYPE *tmp = ConstructSelected<T>(args, *this, name + "." + select_value);
                    delete tmp;
                };

                (traverse_one.template operator()<IC>(), ...);
                restore_select_value();
                ic_eta = nullptr;
                return;
            }

            // if there is only one IC arg provided, we don't need to check the type - assume that
            // it is the default.
            if constexpr (sizeof...(IC) == 0)
            {
                using first_IC = std::tuple_element_t<0, std::tuple<IC...>>;
                ic_eta = ConstructSelected<first_IC>(args, *this, name + "." + std::string(first_IC::name));
            }
            // otherwise, check the type.
            else
            {
                std::string type = "";
                this->query_required(type_name, type);
                bool matched = ((   type == IC::name 
                                    ? (ic_eta = ConstructSelected<IC>(args, *this, name + "." + std::string(IC::name))),
                                    true
                                    : false) || ...);
                if (!matched)
                    if (!InTraversalMode())
                        Util::ParmParseException(INFO,full(name), type, " not a valid type for ", name);
            }
        }
        catch (std::runtime_error &e)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,getPrefix(),e.what());
        }
        catch (...)
        {
            if (!InTraversalMode())
                Util::ParmParseException(INFO,getPrefix());
        }
    }

    //
    // Identical to the above, except the first class is the default. A first
    // class without a ParmParse constructor is constructed directly when
    // no type is specified.
    //
    template<typename FirstIC, typename... IC, typename PTRTYPE>
    void select_default (   std::string name, PTRTYPE *& ic_eta,
                            const std::source_location &location = std::source_location::current())
    {
        auto args = forward_args();
        select_default<FirstIC, IC...>(std::move(name), ic_eta, args, location);
    }

    template<typename FirstIC, typename... IC, typename... Args, typename PTRTYPE>
    void select_default (   std::string name, PTRTYPE *& ic_eta, ForwardArgs<Args...> args,
                            const std::source_location &location = std::source_location::current())
    {
        std::string type = "";
        std::string type_name = name + ".type";
        constexpr bool has_unnamed_default =
            !(  HasSelectName<FirstIC>::value &&
                HasSelectedConstructor<FirstIC, Args...>);
        std::vector<std::string> possiblevals;
        AppendSelectName<FirstIC>(possiblevals, args);
        (AppendSelectName<IC>(possiblevals, args), ...);

        if (InTraversalMode())
        {
            {
                RecordInputAndContinue( type_name, location, possiblevals,
                                        std::nullopt, has_unnamed_default);
            }

            std::string original_type;
            bool had_original_type = this->contains(type_name.c_str());
            if (had_original_type)
                amrex::ParmParse::query(type_name.c_str(), original_type);

            auto set_select_value = [&](const std::string &select_value)
            {
                this->remove(type_name.c_str());
                this->add(type_name.c_str(), select_value);
            };

            auto restore_select_value = [&]()
            {
                this->remove(type_name.c_str());
                if (had_original_type)
                    this->add(type_name.c_str(), original_type);
            };

            auto traverse_one = [&]<typename T>()
            {
                if constexpr (HasSelectName<T>::value && HasSelectedConstructor<T, Args...>)
                {
                    std::string select_value = T::name;
                    set_select_value(select_value);
                    PrintTraversalBranch(type_name, select_value);
                    TraversalBranchScope branch(*this, type_name, select_value);
                    PTRTYPE *tmp = ConstructSelected<T>(args, *this, name + "." + select_value);
                    delete tmp;
                }
            };

            traverse_one.template operator()<FirstIC>();
            (traverse_one.template operator()<IC>(), ...);
            restore_select_value();
            ic_eta = nullptr;
            return;
        }

        if constexpr (  HasSelectName<FirstIC>::value &&
                        HasSelectedConstructor<FirstIC, Args...>)
        {
            this->query_default(type_name, type, std::string(FirstIC::name));
        }
        else
        {
            if (!this->contains(type_name.c_str()))
            {
                ic_eta = ConstructSelectDefault<FirstIC>(args);
                return;
            }
            this->query_required(type_name, type);
        }

        bool matched = TrySelectNamed<FirstIC>(type, name, ic_eta, args, *this);
        matched = matched || (TrySelectNamed<IC>(type, name, ic_eta, args, *this) || ...);

        if (!matched && !InTraversalMode())
            Util::ParmParseException(INFO,full(name), type," not a valid type for ",name);
    }


    //
    //
    template<typename... IC, typename PTRTYPE>
    void select_enumerate ( std::string a_name, std::vector<PTRTYPE*> & value,
                            const std::source_location &location = std::source_location::current())
    {
        auto args = forward_args();
        select_enumerate<IC...>(std::move(a_name), value, args, location);
    }

    template<typename... IC, typename... Args, typename PTRTYPE>
    void select_enumerate ( std::string a_name, std::vector<PTRTYPE*> & value, ForwardArgs<Args...> args,
                            const std::source_location &location = std::source_location::current())
    {
        RecordInputAndContinue(a_name, location);
        value.clear();

        if (InTraversalMode())
        {
            PTRTYPE *tmp = nullptr;
            const std::string template_name = a_name + "0";
            this->select<IC...>(template_name, tmp, args, location);
            InputScraper::CaptureSequenceTemplate(*this, a_name, template_name);
            return;
        }

        //
        // If only one is present with no subscript, then read
        // it only and return.
        //
        std::string name = a_name;
        if (this->contains(name))
        {
            PTRTYPE *tmp;
            this->select<IC...>(a_name, tmp, args, location);
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
            if (!InTraversalMode())
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
                this->select<IC...>(name, tmp, args, location);
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
            if (!InTraversalMode())
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

    template <typename... OBJ, typename CLASS>
        requires (!std::is_pointer_v<std::remove_reference_t<CLASS>>)
    void select(std::string name, CLASS& value,
                const std::source_location &location = std::source_location::current())
    {
#if !AMREX_DEVICE_COMPILE
        auto args = forward_args();
        select<OBJ...>(std::move(name), value, args, location);
#else
        Util::IgnoreUnused(name, value, location);
#endif
    }

    template <typename... OBJ, typename CLASS, typename... Args>
        requires (!std::is_pointer_v<std::remove_reference_t<CLASS>>)
    void select(std::string name, CLASS& value, ForwardArgs<Args...> args,
                const std::source_location &location = std::source_location::current())
    {
#if !AMREX_DEVICE_COMPILE
        pushPrefix(name);
        if (InTraversalMode())
        {
            std::vector<std::string> possiblevals = {std::string(OBJ::name)...};
            {
                RecordInputAndContinue("type", location, possiblevals);
            }

            std::string original_type;
            bool had_original_type = this->contains("type");
            if (had_original_type)
                amrex::ParmParse::query("type", original_type);
            int original_selected = value.selected;

            static_polymorphism_traverser<CLASS, 0, OBJ...>(value, args);

            this->remove("type");
            if (had_original_type)
                this->add("type", original_type);
            value.selected = original_selected;
        }
        else
            static_polymorphism_parser<CLASS, 0, OBJ...>(value, args);
        popPrefix();
#else
        Util::IgnoreUnused(name, value, args, location);
#endif
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
    void static_polymorphism_traverser(CLASS& value, ForwardArgs<Args...> &args)
    {
        if constexpr (N < sizeof...(OBJ))
        {
            std::string type = value.names[N];
            this->remove("type");
            this->add("type", type);
            value.selected = N;

            PrintTraversalBranch("type", type);
            {
                TraversalBranchScope branch(*this, "type", type);
                pushPrefix(type);
                std::apply(
                    [&](auto &...parse_args)
                    {
                        std::get<N>(value.obj).Parse(
                            std::get<N>(value.obj), *this, parse_args...);
                    },
                    args.values);
                popPrefix();
            }

            static_polymorphism_traverser<CLASS, N+1, OBJ...>(
                value, args);
        }
    }

    template <typename CLASS, int N, typename... OBJ, typename... Args>
    void static_polymorphism_parser(CLASS& value, ForwardArgs<Args...> &args)
    {
        if constexpr (N == 0)
        {
            std::string type;
            query_default("type",type,value.names[0]);

            for (unsigned int i = 0; i < sizeof...(OBJ); i++)
                if (type == value.names[i])
                    value.selected = i;

            if (value.selected < 0)
                if (!InTraversalMode())
                    Util::ParmParseException(   INFO,getPrefix(),
                                                "Error reading " + getPrefix() +
                                                ", invalid type " + type);

            
            pushPrefix(type);
        }

        if constexpr (N < sizeof...(OBJ))
        {
            if (value.selected == N)
            {
                std::apply(
                    [&](auto &...parse_args)
                    {
                        std::get<N>(value.obj).Parse(
                            std::get<N>(value.obj), *this, parse_args...);
                    },
                    args.values);
                popPrefix();
                return;
            }
            else 
                return static_polymorphism_parser<CLASS, N+1, OBJ...>(
                    value, args);
        }
        else Util::Abort(INFO);
    }

public:

    template <int N>
    void query_exactly( std::vector<std::string> names, std::pair<std::string, Set::Scalar> values[N],
                        std::vector<Unit> units = std::vector<Unit>(),
                        const std::source_location &location = std::source_location::current())
    {
        for (const auto &name : names)
            RecordInputAndContinue(name, location, names);
        std::vector<std::string> unit_names;
        for (const auto &unit : units)
            unit_names.push_back(unit.normalized_unitstring());
        RecordConstraint("exactly", N, names, unit_names, location);

        if (InTraversalMode())
        {
            int cnt = 0;
            for (unsigned int n = 0; n < names.size() && cnt < N; n++)
            {
                if (!contains(names[n].c_str())) continue;

                values[cnt].first = names[n];
                if (units.size() == 0)
                    amrex::ParmParse::query(names[n].c_str(), values[cnt].second);
                else
                {
                    std::string strvalue;
                    amrex::ParmParse::query(names[n].c_str(), strvalue);
                    Unit read = Unit::Parse(strvalue);
                    values[cnt].second = read.normalized_value();
                }
                cnt++;
            }
            for (unsigned int n = 0; n < names.size() && cnt < N; n++)
            {
                if (contains(names[n].c_str())) continue;

                values[cnt].first = names[n];
                values[cnt].second = 1.0;
                cnt++;
            }
            return;
        }

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
            if (!InTraversalMode())
                Util::ParmParseException(INFO,names_str);
        }
    }
};
}
#endif
