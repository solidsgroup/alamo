#ifndef UTIL_UTIL_H
#define UTIL_UTIL_H

#include <iostream>
#include <iomanip>
#include <sys/ioctl.h>
#include <unistd.h>
#include <iterator>
#include <algorithm>
#include <complex>
#include <stdlib.h>

#include <string>

#include "AMReX.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include "Color.H"
#define INFO __FILE__, __func__, __LINE__
#define TEST(x) #x, x


#ifndef ALAMO_SINGLE_DEFINITION
#define ALAMO_SINGLE_DEFINITION AMREX_ATTRIBUTE_WEAK
#endif


/// \brief A collection of utility routines
namespace Util
{
/// \brief A collection of string operations
namespace String
{
int ReplaceAll(std::string &str, const std::string before, const std::string after);
int ReplaceAll(std::string &str, const char before, const std::string after);
std::string Wrap(std::string text, unsigned per_line);
std::vector<std::string> Split(std::string &str, const char token = ' ');
bool Contains(std::string &str, const std::string find);
template<class T>
T Parse(std::string);
template<>
std::complex<int> Parse(std::string);
}

extern std::string globalprefix;

std::string GetFileName();
void CopyFileToOutputDir(std::string a_path, bool fullpath = true, std::string prefix = "");

std::pair<std::string,std::string> GetOverwrittenFile();

void SignalHandler(int s);

void Initialize (int argc, char* argv[]);
void Initialize ();

void Finalize ();

void Terminate(const char * msg, int signal, bool backtrace);

void Abort (const char * msg);

template<typename... Args>
AMREX_FORCE_INLINE
void Assert (std::string file, std::string func, int line, std::string smt, bool pass, Args const &... args)
{
    if (pass) return;

    std::ostringstream infostream;
    infostream << globalprefix;
    infostream << Color::Bold << Color::FG::Red << "ABORT("<< amrex::ParallelDescriptor::MyProc()<<")" << Color::Reset << ":   ";
    infostream << Color::FG::Red << file << Color::Reset << ":" << line << " ";
    infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
    infostream << "Assertion failed: " << Color::BG::DarkGray << Color::Underlined << Color::Bold <<  Color::FG::LightGreen << smt << Color::Reset;
    std::string info = infostream.str();


    std::ostringstream messagestream;
    using List= int[];
    (void)List{0, ( (void)(messagestream << args), 0 ) ... };
    std::string message = messagestream.str();
    String::ReplaceAll(message,'\n',"\n"+info);

    std::cout << info << message << std::endl;

    Abort("Fatal Error");
}

template<typename... Args>
void Abort (std::string file, std::string func, int line, Args const &... args)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ostringstream infostream;
        infostream << globalprefix;
        infostream << Color::Bold << Color::FG::Red << "ABORT" << Color::Reset << ":   ";
        infostream << Color::FG::Red << file << Color::Reset << ":" << line << " ";
        infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
        std::string info = infostream.str();

        std::ostringstream messagestream;
        using List= int[];
        (void)List{0, ( (void)(messagestream << args), 0 ) ... };
        std::string message = messagestream.str();//String::Wrap(messagestream.str(),150);
        //std::string messageorg = message;
        String::ReplaceAll(message,'\n',"\n"+info);
        std::cout << info << message << std::endl;
    }
    Abort("Fatal Error");
}

template<typename... Args>
void ParallelAbort (std::string file, std::string func, int line, Args const &... args)
{

    std::ostringstream infostream;
    infostream << globalprefix;
    infostream << Color::Bold << Color::FG::Red << "ABORT("<< amrex::ParallelDescriptor::MyProc()<<")" << Color::Reset << ":   ";
    infostream << Color::FG::Red << file << Color::Reset << ":" << line << " ";
    infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
    std::string info = infostream.str();

    std::ostringstream messagestream;
    using List= int[];
    (void)List{0, ( (void)(messagestream << args), 0 ) ... };
    std::string message = messagestream.str();//String::Wrap(messagestream.str(),150);
    //std::string messageorg = message;
    String::ReplaceAll(message,'\n',"\n"+info);
    std::cout << info << message << std::endl;

    Abort("Fatal Error");
}

template<typename... Args>
void Message (std::string file, std::string func, int line, Args const &... args)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ostringstream infostream;
        infostream << globalprefix;
        infostream << Color::Bold << Color::FG::Blue << "MESSAGE" << Color::Reset << ": ";
        infostream << Color::FG::LightBlue << file << Color::Reset << ":" << line << " ";
        infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
        std::string info = infostream.str();

        std::ostringstream messagestream;
        using List= int[];
        (void)List{0, ( (void)(messagestream << args), 0 ) ... };
        std::string message = messagestream.str();//String::Wrap(messagestream.str(),150);

        String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;
    }
}

template<typename... Args>
void ParallelMessage (std::string file, std::string func, int line, Args const &... args)
{
    std::ostringstream infostream;
    infostream << globalprefix;
    infostream << Color::Bold << Color::FG::Blue << "MESSAGE("<< amrex::ParallelDescriptor::MyProc()<<")" << Color::Reset << ": ";
    infostream << Color::FG::LightBlue << file << Color::Reset << ":" << line << " ";
    infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
    std::string info = infostream.str();
    std::ostringstream messagestream;
    using List= int[];
    (void)List{0, ( (void)(messagestream << args), 0 ) ... };
    std::string message = messagestream.str();//String::Wrap(messagestream.str(),150);
    String::ReplaceAll(message,'\n',"\n"+info);
    std::cout << info << message << std::endl;
}

template<typename... Args>
void Warning (std::string file, std::string func, int line, Args const &... args)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ostringstream infostream;
        infostream << globalprefix;
        infostream << Color::Bold << Color::FG::Yellow << "WARNING" << Color::Reset << ": ";
        infostream << Color::FG::Yellow << file << Color::Reset << ":" << line << " ";
        infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
        std::string info = infostream.str();

        std::ostringstream messagestream;
        using List= int[];
        (void)List{0, ( (void)(messagestream << args), 0 ) ... };
        std::string message = messagestream.str();///String::Wrap(messagestream.str(),150);

        String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;
    }
}


template<typename... Args>
void Exception (std::string file, std::string func, int line, Args const &... args)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ostringstream infostream;
        infostream << globalprefix;
        infostream << Color::Bold << Color::FG::Red << "EXCEPTION" << Color::Reset << ": ";
        infostream << Color::FG::LightGray << file << ":" << line << " (" << func << Color::Reset << ") ";
        std::string info = infostream.str();

        std::ostringstream messagestream;
        messagestream << Color::FG::Red;
        using List= int[];
        (void)List{0, ( (void)(messagestream << args), 0 ) ... };
        messagestream << Color::Reset;

        std::string message = messagestream.str();

        Util::String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;

        std::throw_with_nested(std::runtime_error("IO::Exception"));
    }
}

template<typename... Args>
AMREX_FORCE_INLINE
void AssertException (std::string file, std::string func, int line, std::string smt, bool pass, Args const &... args)
{
    if (pass) return;

    std::ostringstream infostream;
    infostream << globalprefix;
    infostream << Color::Bold << Color::FG::Red << "ABORT("<< amrex::ParallelDescriptor::MyProc()<<")" << Color::Reset << ":   ";
    infostream << Color::FG::Red << file << Color::Reset << ":" << line << " ";
    infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
    infostream << "Assertion failed: " << Color::BG::DarkGray << Color::Underlined << Color::Bold <<  Color::FG::LightGreen << smt << Color::Reset;
    std::string info = infostream.str();


    std::ostringstream messagestream;
    using List= int[];
    (void)List{0, ( (void)(messagestream << args), 0 ) ... };
    std::string message = messagestream.str();
    String::ReplaceAll(message,'\n',"\n"+info);

    std::cout << info << message << std::endl;

    Exception(file,func,line,"Fatal Error");
}

template<typename... Args>
void ParmParseException (   std::string file, std::string func, int line,
                            std::string file2, std::string /*func2*/, int line2,
                            std::string fullname,
                            Args const &... args)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ostringstream infostream;
        infostream << globalprefix;
        infostream << Color::Bold << Color::FG::Red << "EXCEPTION" << Color::Reset << ": ";
        infostream << Color::FG::LightGray << file << ":" << line << " (" << func << Color::Reset << ") ";
        std::string info = infostream.str();

        std::ostringstream messagestream;
        messagestream << Color::FG::Red;
        using List= int[];
        (void)List{0, ( (void)(messagestream << args), 0 ) ... };
        messagestream << Color::Reset;

        const std::size_t nargs = sizeof...(Args);

        if (nargs) messagestream << "\n";
        if (fullname != "")
            messagestream << "while reading " << Color::FG::Yellow << fullname << Color::Reset;
        if (line2 >= 0)
        {
            messagestream << " at " << file2 << ":" << line2;
        }

        std::string message = messagestream.str();

        Util::String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;

        std::throw_with_nested(std::runtime_error("IO::ParmParse Exception"));
    }
}


std::pair<std::string,std::string>
CreateCleanDirectory (const std::string &path, bool callbarrier = true);


/// \brief A collection of I/O routines for Unit Tests
namespace Test
{
int Message(std::string testname);
int Message(std::string testname, int passed);
int SubMessage(std::string testname, int passed);
int SubFinalMessage(int failed);
}

//
// These two functions are not really necessary anymore.
// TODO: remove RealFillBoundary and all references to them.
//
template<class T>
AMREX_FORCE_INLINE
void RealFillBoundary(amrex::FabArray<amrex::BaseFab<T>> &a_mf,const amrex::Geometry &/*a_geom*/)
{
    BL_PROFILE("Util::RealFillBoundary");
    a_mf.setMultiGhost(true);
    a_mf.FillBoundary();
}
AMREX_FORCE_INLINE
void RealFillBoundary(amrex::MultiFab &a_mf, const amrex::Geometry &/*a_geom*/)
{
    BL_PROFILE("Util::RealFillBoundary");
    a_mf.setMultiGhost(true);
    a_mf.FillBoundary();
}


void AverageCellcenterToNode(amrex::MultiFab& node_mf, const int &dcomp, const amrex::MultiFab &cell_mf, const int &scomp, const int &ncomp/*, const int ngrow=0*/);




}
#endif
