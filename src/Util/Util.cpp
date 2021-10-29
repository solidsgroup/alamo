#include "Util.H"
#include "Color.H"

#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"

#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
#include "Color.H"
#include <chrono>

#include "Numeric/Stencil.H"

namespace Util
{

std::string filename = "";
std::pair<std::string,std::string> file_overwrite;

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
		IO::FileNameParse(filename);
		// else
		// 	if (amrex::ParallelDescriptor::IOProcessor())
		// 		Util::Abort("No plot file specified! (Specify plot_file = \"plot_file_name\" in input file");
	}
	return filename;
}

std::pair<std::string,std::string> GetOverwrittenFile()
{
	return file_overwrite;
}

void SignalHandler(int s)
{
	if (amrex::ParallelDescriptor::IOProcessor())
	{
		std::string filename = GetFileName();
		IO::Status status = IO::Status::Running;
		if (s == SIGSEGV) status = IO::Status::Segfault;
		else if (s == SIGINT) status = IO::Status::Interrupt;
		if (s == SIGABRT) status = IO::Status::Abort;
		if (filename != "")
			IO::WriteMetaData(filename,status);
	}

#ifdef MEME
	amrex::ParmParse pp;
	if (!pp.contains("nomeme"))
	{			
		time_t timer; time(&timer);
		std::stringstream cmd;
		cmd << "xdg-open " << BUILD_DIR << "/src/Util/Meme/cat0" << (1+((int)timer)%6) << ".gif &";
		std::system(cmd.str().c_str());
		std::cout << Color::Bold << Color::FG::Red << "PROGRAM FAILED!" << Color::Reset << " (Compile without -DMEME, or set nomeme = 1 in the input file to disable this!)";
	}
#endif 

	amrex::BLBackTrace::handler(s);
}


void Initialize ()
{
	int argc = 0;
	char **argv = nullptr;
	Initialize(argc,argv);
}
void Initialize (int argc, char* argv[])
{
	srand (time(NULL));

	// if (argc < 2)
	// {
	// 	std::cout << "No plot file specified!" << std::endl;
	// 	exit(-1);
	// }

	amrex::Initialize(argc, argv);

	amrex::ParmParse pp_amrex("amrex");
	pp_amrex.add("throw_exception",1);
	//amrex.throw_exception=1

	signal(SIGSEGV, Util::SignalHandler);
	signal(SIGINT,  Util::SignalHandler);
	signal(SIGABRT, Util::SignalHandler);

	std::string filename = GetFileName();

	if (amrex::ParallelDescriptor::IOProcessor() && filename != "")
	{
		file_overwrite = Util::CreateCleanDirectory(filename, false);
		IO::WriteMetaData(filename);
	}
}

void Finalize()
{
	std::string filename = GetFileName();
	if (filename != "")
	 	IO::WriteMetaData(filename,IO::Status::Complete);
	amrex::Finalize();
}



void
Abort (const char * msg) { Terminate(msg, SIGABRT, true); }

void
Terminate(const char * /* msg */, int signal, bool /*backtrace*/)
{
	SignalHandler(signal);
}

std::pair<std::string,std::string>
CreateCleanDirectory (const std::string &path, bool callbarrier)
{
	std::pair<std::string,std::string> ret("","");

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
			ret.first = path;
			ret.second = newoldname;
		}
		if( ! amrex::UtilCreateDirectory(path, 0755)) {
			amrex::CreateDirectoryFailed(path);
		}
	}
	if(callbarrier) {
		// Force other processors to wait until directory is built.
		amrex::ParallelDescriptor::Barrier("amrex::UtilCreateCleanDirectory");
	}
	return ret;
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
std::vector<std::string> Split(std::string &str, const char delim)
{
	std::vector<std::string> ret;
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, delim)) {
		ret.push_back(item);
	}
	return ret;
}
bool Contains(std::string &str, const std::string find)
{
  if (str.find(find) != std::string::npos) return true;
  else return false;
}

template<>
std::complex<int> Parse(std::string input)
{
	int re=0, im=0;

	ReplaceAll(input, "+", " +");
	ReplaceAll(input, "-", " -");
	std::vector<std::string> tokens = Split(input,' ');
	for (unsigned int i = 0; i < tokens.size(); i++)
	{
		if(tokens[i]=="") continue;
		if(Contains(tokens[i],"i"))
		{
			ReplaceAll(tokens[i],"i","");
			im += std::stoi(tokens[i]);
		}
		else 
		{
			re += std::stoi(tokens[i]);
		}
	}
	return std::complex<int>(re,im);
}


}


namespace Test
{
int Message(std::string testname)
{
	if (amrex::ParallelDescriptor::IOProcessor())
		std::cout << std::left
			  << Color::FG::White << Color::Bold << testname << Color::Reset << std::endl;
	return 0;
}
int Message(std::string testname, int failed)
{
	if (amrex::ParallelDescriptor::IOProcessor())
	{
		winsize w;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		std::stringstream ss;
		if (!failed)
			ss << "[" << Color::FG::Green << Color::Bold << "PASS" << Color::Reset << "]";
		else
			ss << "[" << Color::FG::Red << Color::Bold << "FAIL" << Color::Reset << "]";

		int terminalwidth = 80; //std::min(w.ws_col,(short unsigned int) 100);

		std::cout << std::left
			  << testname 
			  << std::setw(terminalwidth - testname.size() + ss.str().size() - 6)  << std::right << std::setfill('.') << ss.str() << std::endl;
	}
	return failed;
}
int SubMessage(std::string testname, int failed)
{
	if (amrex::ParallelDescriptor::IOProcessor())
	{
		winsize w;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		std::stringstream ss;
		if (!failed)
			ss << "[" << Color::FG::LightGreen << Color::Bold << "PASS" << Color::Reset << "]";
		else
			ss << "[" << Color::FG::Red << Color::Bold << "FAIL" << Color::Reset << "]";

		int terminalwidth = 80; 

		std::cout << std::left
			  << "  ├ "
			  << testname 
			  << std::setw(terminalwidth - testname.size() + ss.str().size() - 12)  << std::right << std::setfill('.') << ss.str() << std::endl;
	}
	return failed;
}
int SubFinalMessage(int failed)
{
	if (amrex::ParallelDescriptor::IOProcessor())
	{
		winsize w;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		std::stringstream ss;
		std::cout << std::left << "  └ ";

		if (!failed)
			std::cout << Color::FG::Green << Color::Bold << failed << " tests failed" << Color::Reset << std::endl;
		else
			std::cout << Color::FG::Red << Color::Bold << failed << " tests failed" << Color::Reset << std::endl;

		
	}
	return failed;
}

}

void AverageCellcenterToNode(amrex::MultiFab& node_mf, const int &dcomp, const amrex::MultiFab &cell_mf, const int &scomp, const int &ncomp/*, const int ngrow=0*/)
{
	Util::Assert(INFO,TEST(dcomp + ncomp <= node_mf.nComp()));
	Util::Assert(INFO,TEST(scomp + ncomp <= cell_mf.nComp()));
	//Util::Assert(INFO,TEST(cell_mf.boxArray() == node_mf.boxArray()));
	Util::Assert(INFO,TEST(cell_mf.DistributionMap() == cell_mf.DistributionMap()));
	Util::Assert(INFO,TEST(cell_mf.nGrow() > 0));
	for (amrex::MFIter mfi(node_mf,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
			amrex::Box bx = mfi.nodaltilebox();
			amrex::Array4<Set::Scalar>       const& node = node_mf.array(mfi);
			amrex::Array4<const Set::Scalar> const& cell = cell_mf.array(mfi);
			for (int n = 0; n < ncomp; n++)
				amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
										   node(i,j,k,dcomp+n) = Numeric::Interpolate::CellToNodeAverage(cell,i,j,k,scomp+n);
									   });
	}
}


}
