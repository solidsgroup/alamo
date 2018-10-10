#ifndef UTIL_UTIL_H
#define UTIL_UTIL_H

#include <string>

#include "AMReX.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include "Color.H"
#define INFO __FILE__, __func__, __LINE__

namespace Util
{
namespace String
{
int ReplaceAll(std::string &str, const std::string before, const std::string after);
int ReplaceAll(std::string &str, const char before, const std::string after);
}


std::string GetFileName();

void SignalHandler(int s);

void Initialize (int argc, char* argv[]);

void Finalize ();

void Terminate(const char * msg, int signal, bool backtrace);

void Abort (const char * msg);

template<typename... Args>
void Abort (std::string file, std::string func, int line, Args const &... args)
{
	std::ostringstream infostream;
	infostream << Color::Bold << Color::FG::Red << "ABORT" << Color::Reset << ":   ";
	infostream << Color::FG::Red << file << Color::Reset << ":" << line << " ";
	infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
	std::string info = infostream.str();

	std::ostringstream messagestream;
	using List= int[];
	(void)List{0, ( (void)(messagestream << args), 0 ) ... };
	std::string message = messagestream.str();
	std::string messageorg = message;

	String::ReplaceAll(message,'\n',"\n"+info);

	if (amrex::ParallelDescriptor::IOProcessor())
		std::cout << info << message << std::endl;

	Abort(messageorg.c_str());
}

template<typename... Args>
void Message (std::string file, std::string func, int line, Args const &... args)
{
	std::ostringstream infostream;
	infostream << Color::Bold << Color::FG::Blue << "MESSAGE" << Color::Reset << ": ";
	infostream << Color::FG::LightBlue << file << Color::Reset << ":" << line << " ";
	infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
	std::string info = infostream.str();

	std::ostringstream messagestream;
	using List= int[];
	(void)List{0, ( (void)(messagestream << args), 0 ) ... };
	std::string message = messagestream.str();

	String::ReplaceAll(message,'\n',"\n"+info);

	if (amrex::ParallelDescriptor::IOProcessor())
		std::cout << info << message << std::endl;
}

template<typename... Args>
void Warning (std::string file, std::string func, int line, Args const &... args)
{
	std::ostringstream infostream;
	infostream << Color::Bold << Color::FG::Yellow << "WARNING" << Color::Reset << ": ";
	infostream << Color::FG::Yellow << file << Color::Reset << ":" << line << " ";
	infostream << "(" << Color::FG::LightGreen << func << Color::Reset << ") ";
	std::string info = infostream.str();

	std::ostringstream messagestream;
	using List= int[];
	(void)List{0, ( (void)(messagestream << args), 0 ) ... };
	std::string message = messagestream.str();

	String::ReplaceAll(message,'\n',"\n"+info);

	if (amrex::ParallelDescriptor::IOProcessor())
		std::cout << info << message << std::endl;
}

void CreateCleanDirectory (const std::string &path, bool callbarrier = true);

}
#endif
