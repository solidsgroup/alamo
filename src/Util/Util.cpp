#include "Util.H"
#include "Color.H"

#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"

#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
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
			Abort("plot_file specified in too many locations");
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
				Abort("No plot file specified! (Specify plot_file = \"plot_file_name\" in input file");
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
Abort (std::string msg) { Abort(msg.c_str()); }

void
Abort (const char * msg) { Terminate(msg, SIGABRT, true); }

void
Terminate(const char * msg, int signal, bool backtrace)
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



}
