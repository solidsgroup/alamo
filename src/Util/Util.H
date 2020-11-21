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

std::string GetFileName();

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
	a_mf.setMultiGhost(true);
	a_mf.FillBoundary();
}
AMREX_FORCE_INLINE
void RealFillBoundary(amrex::MultiFab &a_mf, const amrex::Geometry &/*a_geom*/)
{
	a_mf.setMultiGhost(true);
	a_mf.FillBoundary();
}



}
#endif
