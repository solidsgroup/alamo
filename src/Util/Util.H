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
#include <string_view>

#include <string>
#include "Util/String.H"

#include "AMReX.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_iMultiFab.H>

#include "Color.H"
#define INFO __FILE__, __func__, __LINE__
#define TEST(x) #x, x


#ifndef ALAMO_SINGLE_DEFINITION
#define ALAMO_SINGLE_DEFINITION AMREX_ATTRIBUTE_WEAK
#endif

namespace Util
{

extern std::string globalprefix;
extern bool initialized;
extern bool finalized;


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
    Util::String::ReplaceAll(message,'\n',"\n"+info);

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
        Util::String::ReplaceAll(message,'\n',"\n"+info);
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
    Util::String::ReplaceAll(message,'\n',"\n"+info);
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

        Util::String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;
    }
}

#ifdef AMREX_DEBUG
template<typename... Args>
void DebugMessage (std::string file, std::string func, int line, Args const &... args)
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

        Util::String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;
    }
}
#else
template<typename... Args>
void DebugMessage (std::string, std::string, int, Args const &... )
{
}
#endif




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
    Util::String::ReplaceAll(message,'\n',"\n"+info);
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

        Util::String::ReplaceAll(message,'\n',"\n"+info);

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

    if (amrex::ParallelDescriptor::IOProcessor())
    {
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
        Util::String::ReplaceAll(message,'\n',"\n"+info);

        std::cout << info << message << std::endl;

        Exception(file,func,line,"Fatal Error");
    }
}



template<typename... Args>
void ParmParseException (   std::string file, std::string func, int line,
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
void SubWarning(std::string testname);
int SubFinalMessage(int failed);
}

//
// These two functions are not really necessary anymore.
// TODO: remove RealFillBoundary and all references to them.
//
template<class T>
AMREX_FORCE_INLINE
void RealFillBoundary(amrex::FabArray<amrex::BaseFab<T>> &a_mf,const amrex::Geometry &/*a_geom*/, const int nghost=2)
{
    // We want the default to setMultiGhost(true) if nghost isn't explicitly passed, so nghost should default to a number greater than 1
    BL_PROFILE("Util::RealFillBoundary");
    if (nghost <= 1) {}
    else a_mf.setMultiGhost(true);
    a_mf.FillBoundary();
}
AMREX_FORCE_INLINE
void RealFillBoundary(amrex::MultiFab &a_mf, const amrex::Geometry &/*a_geom*/, const int nghost=2)
{
    // We want the default to setMultiGhost(true) if nghost isn't explicitly passed, so nghost should default to a number greater than 1
    BL_PROFILE("Util::RealFillBoundary");
    if (nghost <= 1) {}
    else a_mf.setMultiGhost(true);
    a_mf.FillBoundary();
}

template<class T>
AMREX_FORCE_INLINE
void Probe(
    std::string file, std::string func, int line,
    amrex::FabArray<amrex::BaseFab<T>> &a_mf, int i, int j, int k, int n)
{
    // Loop over boxes owned by this rank
    for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& validbox = mfi.validbox();
        const amrex::Box& box = mfi.grownnodaltilebox();

        if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
        {
            const auto& arr = a_mf.const_array(mfi);
            if (validbox.strictly_contains(i,j,k))
                Util::ParallelMessage(file,func,line, "inside: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else if (validbox.contains(i,j,k))
                Util::ParallelMessage(file,func,line, "border: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else
                Util::ParallelMessage(file,func,line, "ghost:  (",i,",",j,",",k,") ",arr(i, j, k, n));
        }
    }
}

AMREX_FORCE_INLINE
void Probe(
    std::string file, std::string func, int line,
    const amrex::MultiFab &a_mf, int _i, int _j, int /*_k*/, int n,
    int pad = 0)
{
    const int t_notexist = 0;
    const int t_inside = 1;
    const int t_border = 2;
    const int t_ghost = 3;
    

    for (int rank = 0; rank < amrex::ParallelDescriptor::NProcs() ; rank++)
    {
        amrex::ParallelDescriptor::Barrier();
        if (rank != amrex::ParallelDescriptor::MyRankInTeam()) continue;

        // Loop over boxes owned by this rank
        for (amrex::MFIter mfi(a_mf,false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& validbox = mfi.validbox();
            const amrex::Box& box = mfi.grownnodaltilebox();
            bool contains = false;

            std::vector<std::vector<Set::Scalar>> vals(2*pad + 1);
            std::vector<std::vector<int>> types(2*pad + 1);
            for (auto & val : vals) val.resize(2*pad + 1,NAN);
            for (auto & type : types) type.resize(2*pad + 1,t_notexist);

            for (int I = 0; I <= 2*pad; I++)
            {
                for (int J = 0; J <= 2*pad; J++)
                {
                    int i = (_i - pad) + I;
                    int j = (_j - pad) + J;
                
                    if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
                    {
                        const auto& arr = a_mf.const_array(mfi);
                        vals[I][J] = arr(i,j,0,n);
                        if (validbox.strictly_contains(i,j,0)) types[I][J] = t_inside;
                        else if (validbox.contains(i,j,0))     types[I][J] = t_border;
                        else types[I][J] = t_ghost;
                        contains = true;
                    }
                }
            }

            if (contains)
            {
                std::stringstream ss;
                ss << "  " << std::setw(10) << Color::FG::LightYellow;
                for (int I = 0; I <= 2*pad; I++)
                {
                    int i = (_i - pad) + I;
                    ss << "i=" << std::left << std::setw(13)  << i;
                }
                Util::ParallelMessage(file,func,line,ss.str());
                

                for (int J = 2*pad; J >= 0; J--)
                {
                    int j = (_j - pad) + J;
                
                    std::stringstream ss;
                    ss << Color::FG::LightYellow << "j=" << std::left << std::setw(5) << j << Color::Reset;
                    for (int I = 0; I <= 2*pad; I++)
                    {
                        if (types[I][J] == t_ghost) ss << Color::Dim;
                        if (types[I][J] == t_border) ss << Color::FG::Green;
                        if (types[I][J] == t_inside) ss << Color::FG::Blue;
                        if (types[I][J] == t_notexist) ss << Color::FG::Red;
                        ss << std::left << std::setw(15) << vals[I][J];
                        ss << Color::Reset;
                    }
                    Util::ParallelMessage(file,func,line,ss.str());
                }
            }
        }
    }
}


AMREX_FORCE_INLINE
void Probe(
    std::string file, std::string func, int line,
    amrex::iMultiFab &a_mf, int i, int j, int k, int n)
{
    // Loop over boxes owned by this rank
    for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& validbox = mfi.validbox();
        const amrex::Box& box = mfi.grownnodaltilebox();

        if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
        {
            const auto& arr = a_mf.const_array(mfi);
            if (validbox.strictly_contains(i,j,k))
                Util::ParallelMessage(file,func,line, "inside: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else if (validbox.contains(i,j,k))
                Util::ParallelMessage(file,func,line, "border: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else
                Util::ParallelMessage(file,func,line, "ghost:  (",i,",",j,",",k,") ",arr(i, j, k, n));
        }
    }
}

void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, std::string desc, const amrex::Box domain);

AMREX_FORCE_INLINE
void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, const amrex::Box domain)
{
    Compare(file,func,line,mf,"",domain);
}
AMREX_FORCE_INLINE
void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, std::string desc="")
{    
    Compare(file,func,line,mf,desc,amrex::Box::TheUnitBox());
}


void AverageCellcenterToNode(amrex::MultiFab& node_mf, const int &dcomp, const amrex::MultiFab &cell_mf, const int &scomp, const int &ncomp/*, const int ngrow=0*/);

// Prevent compiler from complaining about unused variables
template <class... Ts>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void IgnoreUnused (const Ts&...) {}



}

#endif
