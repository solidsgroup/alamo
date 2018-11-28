#include "Util.H"
#include "Color.H"

#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"

#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
#include "Color.H"
#include <chrono>

namespace Util
{

std::string filename = "";

std::string GetFileName()
{
	if (filename == "")
	{
		amrex::ParmParse pp;
		amrex::ParmParse pp_amr("amr");

		if (pp_amr.contains("plot_file") && pp.contains("plot_file"))
			Util::Abort("plot_file specified in too many locations");
		else if (pp_amr.contains("plot_file"))
		{
			if (amrex::ParallelDescriptor::IOProcessor())
				amrex::Warning("amr.plot_file will be depricated; use plot_file instead");
			pp_amr.query("plot_file", filename);

		}
		else if (pp.contains("plot_file"))
		{
			pp.query("plot_file", filename);
		}
		else
			if (amrex::ParallelDescriptor::IOProcessor())
				Util::Abort("No plot file specified! (Specify plot_file = \"plot_file_name\" in input file");
	}
	return filename;
}

void SignalHandler(int s)
{
	if (amrex::ParallelDescriptor::IOProcessor())
	{
		std::string filename = GetFileName();
		IO::Status status;
		if (s == SIGSEGV) status = IO::Status::Segfault;
		else if (s == SIGINT) status = IO::Status::Interrupt;
		if (s == SIGABRT) status = IO::Status::Abort;
		IO::WriteMetaData(filename,status);
	}

	amrex::BLBackTrace::handler(s);
}


void Initialize (int argc, char* argv[])
{
	srand (time(NULL));

	if (argc < 2)
	{
		std::cout << "No plot file specified!" << std::endl;
		exit(-1);
	}

	amrex::Initialize(argc, argv);

	amrex::ParmParse pp_amrex("amrex");
	pp_amrex.add("throw_exception",1);
	//amrex.throw_exception=1

	signal(SIGSEGV, Util::SignalHandler);
	signal(SIGINT,  Util::SignalHandler);
	signal(SIGABRT, Util::SignalHandler);

	std::string filename = GetFileName();

	if (amrex::ParallelDescriptor::IOProcessor())
	{
		Util::CreateCleanDirectory(filename, false);
		IO::WriteMetaData(filename);
	}
}

void Finalize()
{
	std::string filename = GetFileName();
	IO::WriteMetaData(filename,IO::Status::Complete);
	amrex::Finalize();
}



void
Abort (const char * msg) { Terminate(msg, SIGABRT, true); }

void
Terminate(const char * msg, int signal, bool /*backtrace*/)
{
	amrex::write_to_stderr_without_buffering(msg);
	SignalHandler(signal);
}

void
CreateCleanDirectory (const std::string &path, bool callbarrier)
{
	if(amrex::ParallelDescriptor::IOProcessor()) {
		if(amrex::FileExists(path)) {
			std::time_t t = std::time(0);
			std::tm * now = std::localtime(&t);
			int year = now->tm_year+1900;
			int month = now->tm_mon+1;
			int day = now->tm_mday;
			int hour = now->tm_hour;
			int minute = now->tm_min;
			int second = now->tm_sec;

			std::stringstream ss;
			ss << year
			   << std::setfill('0') << std::setw(2) << month
			   << std::setfill('0') << std::setw(2) << day
			   << std::setfill('0') << std::setw(2) << hour
			   << std::setfill('0') << std::setw(2) << minute
			   << std::setfill('0') << std::setw(2) << second;

			std::string newoldname(path + ".old." + ss.str());
			if (amrex::system::verbose) {
				amrex::Print() << "Util::CreateCleanDirectory():  " << path
					       << " exists.  Renaming to:  " << newoldname << std::endl;
			}
			std::rename(path.c_str(), newoldname.c_str());
		}
		if( ! amrex::UtilCreateDirectory(path, 0755)) {
			amrex::CreateDirectoryFailed(path);
		}
	}
	if(callbarrier) {
		// Force other processors to wait until directory is built.
		amrex::ParallelDescriptor::Barrier("amrex::UtilCreateCleanDirectory");
	}
}






namespace String
{
int ReplaceAll(std::string &str, const std::string before, const std::string after)
{
	size_t start_pos = 0;
	while((start_pos = str.find(before, start_pos)) != std::string::npos) {
		str.replace(start_pos, before.length(), after);
		start_pos += after.length();
	}
	return 0;
}
int ReplaceAll(std::string &str, const char before, const std::string after)
{
	size_t start_pos = 0;
	while((start_pos = str.find(before, start_pos)) != std::string::npos) {
		str.replace(start_pos, 1, after);
		start_pos += after.length();
	}
	return 0;
}


std::string Wrap(std::string text, unsigned per_line)
{
	unsigned line_begin = 0;

	while (line_begin < text.size())
	{
		const unsigned ideal_end = line_begin + per_line ;
		unsigned line_end = ideal_end <= text.size() ? ideal_end : text.size()-1;

		if (line_end == text.size() - 1)
			++line_end;
		else if (std::isspace(text[line_end]))
		{
			text[line_end] = '\n';
			++line_end;
		}
		else    // backtrack
		{
			unsigned end = line_end;
			while ( end > line_begin && !std::isspace(text[end]))
				--end;

			if (end != line_begin)                  
			{                                       
				line_end = end;                     
				text[line_end++] = '\n';            
			}                                       
			else                                    
				text.insert(line_end++, 1, '\n');
		}

		line_begin = line_end;
	}

	return text;
}



}


Set::Scalar Random()
{
	return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}

namespace Test
{
int Message(std::string testname, bool passed)
{
	winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	std::stringstream ss;
	if (passed)
		ss << "[" << Color::FG::Green << Color::Bold << "PASS" << Color::Reset << "]";
	else
		ss << "[" << Color::FG::Red << Color::Bold << "FAIL" << Color::Reset << "]";

	std::cout << std::left
		  << Color::FG::White << Color::Bold << testname << Color::Reset
		  << std::setw(w.ws_col - testname.size() + ss.str().size() - 6)  << std::right << std::setfill('.') << ss.str() << std::endl;
	if (passed) return 0;
	else return 1;
}
}



}
