#include "Util.H"
#include "AMReX_Config.H"
#include "AMReX_DistributionMapping.H"
#include "AMReX_VisMF.H"
#include "Color.H"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"

#include "IO/ParmParse.H"
#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
#include "Color.H"
#include "Numeric/Stencil.H"
#include "Util/MPI.H"
#include <mpi.h>

namespace
{
void
ParseArgsError(const std::string &message)
{
    std::cerr << "ERROR: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

bool
IsInputDefinition(const std::string &arg)
{
    return arg == "input" || arg.rfind("input=", 0) == 0;
}

void
RejectParseArgsInputFiles(const std::vector<char*> &argv)
{
    for (std::size_t i = 1; i < argv.size(); ++i)
    {
        const std::string arg(argv[i]);

        if (arg == "--") break;

        if (IsInputDefinition(arg))
        {
            ParseArgsError("--parse-args does not accept input-file directives: " + arg);
        }

        if (arg == "=") continue;
        if (!arg.empty() && arg[0] == '-') continue;
        if (arg.find('=') != std::string::npos) continue;
        if (i + 1 < argv.size() && std::string(argv[i + 1]) == "=") continue;
        if (i > 1 && std::string(argv[i - 1]) == "=") continue;

        ParseArgsError("--parse-args does not accept input files or positional arguments: " + arg);
    }
}

void
InjectParseArgsDefaults()
{
    {
        amrex::ParmParse pp("amr");

        int max_level = 0;
        pp.queryAdd("max_level", max_level);

        std::vector<int> n_cell(AMREX_SPACEDIM, 1);
        pp.queryAdd("n_cell", n_cell);

        int max_grid_size = 1;
        pp.queryAdd("max_grid_size", max_grid_size);

        int blocking_factor = 1;
        pp.queryAdd("blocking_factor", blocking_factor);
    }

    {
        amrex::ParmParse pp("geometry");

        std::vector<double> prob_lo(AMREX_SPACEDIM, 0.0);
        pp.queryAdd("prob_lo", prob_lo, AMREX_SPACEDIM);

        if (!pp.contains("prob_hi") && !pp.contains("prob_extent"))
        {
            std::vector<double> prob_hi(AMREX_SPACEDIM, 1.0);
            pp.addarr("prob_hi", prob_hi);
        }

        std::vector<int> is_periodic(AMREX_SPACEDIM, 0);
        pp.queryAdd("is_periodic", is_periodic, AMREX_SPACEDIM);
    }

    {
        amrex::ParmParse pp;

        std::string stop_time = "1.0";
        pp.queryAdd("stop_time", stop_time);

        std::string timestep = "1.0";
        pp.queryAdd("timestep", timestep);
    }
}
}

namespace Util
{

std::string filename = "";
std::string globalprefix = "";
std::pair<std::string,std::string> file_overwrite;
bool initialized = false;
bool finalized = false;

std::string GetFileName()
{
    if (filename == "")
    {
        IO::ParmParse pp;

        pp.forbid("amr.plot_file","Depricated");

        // Output file path
        pp.query_default("plot_file", filename, "output"); // Name of directory containing all output data

        IO::FileNameParse(filename);
        // else
        // if (amrex::ParallelDescriptor::IOProcessor())
        // Util::Abort("No plot file specified! (Specify plot_file = \"plot_file_name\" in input file");
    }
    return filename;
}
void CopyFileToOutputDir(std::string a_path, bool fullpath, std::string prefix)
{
    if (IO::ParmParse::InTraversalMode()) return;

    try
    {
        if (filename == "")
            Util::Exception(INFO,"Cannot back up files yet because the output directory has not been specified");

        std::string basefilename = std::filesystem::path(a_path).filename();
        std::string absolutepath = std::filesystem::absolute(std::filesystem::path(a_path)).string();
        std::string abspathfilename = absolutepath;
        std::replace(abspathfilename.begin(),abspathfilename.end(),'/','_');
        if (prefix != "")
        {
            abspathfilename = prefix + "__" + abspathfilename;
            basefilename    = prefix + "__" + abspathfilename;
        }

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::string destinationpath;
            if (fullpath) destinationpath = filename+"/"+abspathfilename;
            else          destinationpath = filename+"/"+basefilename;

            // Copy the file where the file name is the absolute path, with / replaced with _
            if (std::filesystem::exists(destinationpath))
                Util::Exception(INFO,"Trying to copy ",destinationpath," but it already exists.");
            std::filesystem::copy_file(a_path,destinationpath);
        }
    }
    catch (std::filesystem::filesystem_error const& ex)
    {
        Util::Exception(INFO,
                        "file system error: \n",
                        "     what():  " , ex.what()  , '\n',
                        "     path1(): " , ex.path1() , '\n',
                        "     path2(): " , ex.path2() , '\n',
                        "     code().value():    " , ex.code().value() , '\n',
                        "     code().message():  " , ex.code().message() , '\n',
                        "     code().category(): " , ex.code().category().name());
    }
}

std::pair<std::string,std::string> GetOverwrittenFile()
{
    return file_overwrite;
}

void SignalHandler(int s)
{
    if (!IO::ParmParse::InTraversalMode() &&
        amrex::ParallelDescriptor::IOProcessor())
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
    IO::ParmParse pp;
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
    initialized = true;
}
void Initialize (int argc, char* argv[])
{
    srand (time(NULL));

    bool parse_args = false;
    std::string parse_args_output = "alamo-inputs.schema.json";
    std::vector<char*> amrex_argv;
    amrex_argv.reserve(argc > 0 ? argc : 0);
    for (int i = 0; i < argc; i++)
    {
        if (std::string(argv[i]) == "--parse-args")
        {
            parse_args = true;
            continue;
        }
        if (std::string(argv[i]) == "--parse-args-output")
        {
            if (i + 1 >= argc)
                ParseArgsError("--parse-args-output requires a file path");
            parse_args_output = argv[++i];
            continue;
        }
        amrex_argv.push_back(argv[i]);
    }

    if (parse_args) RejectParseArgsInputFiles(amrex_argv);

    int amrex_argc = static_cast<int>(amrex_argv.size());
    char **amrex_argv_ptr = amrex_argv.empty() ? nullptr : amrex_argv.data();

    IO::ParmParse::SetTraversalMode(parse_args);
    if (parse_args)
        IO::ParmParse::SetTraversalOutputFile(parse_args_output);

    amrex::Initialize(amrex_argc, amrex_argv_ptr);

    if (parse_args) InjectParseArgsDefaults();

    IO::ParmParse pp;
    pp.add("amrex.throw_exception",1);
    //amrex.throw_exception=1

    signal(SIGSEGV, Util::SignalHandler);
    signal(SIGINT,  Util::SignalHandler);
    signal(SIGABRT, Util::SignalHandler);

    std::string filename = GetFileName();

    if (!IO::ParmParse::InTraversalMode() &&
        amrex::ParallelDescriptor::IOProcessor() && filename != "")
    {
        file_overwrite = Util::CreateCleanDirectory(filename, false);
        IO::WriteMetaData(filename);
    }

    std::string length, time, mass, temperature, current, amount, luminousintensity;
    // Set the system length unit
    pp.query_default("system.length",length,"m");
    // Set the system time unit
    pp.query_default("system.time",time,"s");
    // Set the system mass unit
    pp.query_default("system.mass",mass,"kg");
    // Set the system temperature unit
    pp.query_default("system.temperature",temperature,"K");
    // Set the system current unit
    pp.query_default("system.current",current,"A");
    // Set the system amount unit
    pp.query_default("system.amount",amount,"mol");
    // Set the system luminous intensity unit
    pp.query_default("system.luminousintensity",luminousintensity,"cd");
    try
    {
        Unit::setLengthUnit(length);
        Unit::setTimeUnit(time);
        Unit::setMassUnit(mass);
        Unit::setTemperatureUnit(temperature);
        Unit::setCurrentUnit(current);
        Unit::setAmountUnit(amount);
        Unit::setLuminousIntensityUnit(luminousintensity);

        // Update Constants to desired system units
        Set::Constant::SetGlobalConstants();
    }
    catch (std::runtime_error &e)
    {
        Util::Exception(INFO, "Error in setting system units: ", e.what());
    }

    //
    // This is some logic to unit-ize the geometry.prob_lo, geometry.prob_hi input variables/
    // We also do some checking to make sure the geometry is valid.
    //
    // Note that here, unlike most places, we actually **replace and overwrite** the 
    // geom.prob_* variables, since they are read deep inside amrex infrastructure.
    //
    {
        IO::ParmParse pp("geometry");
        
        std::vector<Set::Scalar> prob_lo, prob_hi;
        // Location of the lower+left+bottom corner
        pp.queryarr_required("prob_lo", prob_lo, Unit::Length());
        // Location of the upper_right_top corner
        pp.queryarr_required("prob_hi", prob_hi, Unit::Length());
        pp.remove("prob_lo");
        pp.remove("prob_hi");

        Util::Assert(   INFO,TEST(prob_lo[0] < prob_hi[0]),
                        "Invalid domain specified: ", prob_lo[0], " < x < ", prob_hi[0], " is incorrect.");
        Util::Assert(   INFO,TEST(prob_lo[1] < prob_hi[1]),
                        "Invalid domain specified: ", prob_lo[0], " < y < ", prob_hi[0], " is incorrect.");
#if AMREX_SPACEDIM>2
        Util::Assert(   INFO,TEST(prob_lo[2] < prob_hi[2]),
                        "Invalid domain specified: ", prob_lo[0], " < z < ", prob_hi[0], " is incorrect.");
#endif

        Util::DebugMessage(INFO,"Domain lower left corner: ", Set::Vector(prob_lo.data()).transpose());
        Util::DebugMessage(INFO,"Domain upper right corenr: ", Set::Vector(prob_hi.data()).transpose());

        pp.addarr("prob_lo",prob_lo);
        pp.addarr("prob_hi",prob_hi);
    }


    // This allows the user to ignore certain arguments that
    // would otherwise cause problems.
    // Most generally this is used in the event of a "above inputs
    // specified but not used" error.
    // The primary purpose of this was to fix those errors that arise
    // in regression tests.

    {
        IO::ParmParse pp;
        std::vector<std::string> ignore;
        if (pp.contains("ignore")) Util::Message(INFO, "Ignore directive detected");
        pp.queryarr("ignore", ignore); // Space-separated list of entries to ignore
        for (unsigned int i = 0; i < ignore.size(); i++)
        {
            Util::Message(INFO, "ignoring ", ignore[i]);
            pp.remove(ignore[i].c_str());
        }
    }
}

void Finalize()
{
    if (IO::ParmParse::InTraversalMode())
    {
        IO::ParmParse::WriteInputTreeJsonFile(IO::ParmParse::TraversalOutputFile());
    }
    else
    {
        std::string filename = GetFileName();
        if (filename != "")
            IO::WriteMetaData(filename,IO::Status::Complete);
    }
    amrex::Finalize();
    finalized = true;
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
void SubWarning(std::string testname)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        winsize w;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
        std::stringstream ss;
        ss << "[" << Color::FG::LightYellow << Color::Bold << "WARN" << Color::Reset << "]";

        int terminalwidth = 80; 

        std::cout << std::left
            << "  ├ "
            << testname 
            << std::setw(terminalwidth - testname.size() + ss.str().size() - 12)  << std::right << std::setfill('.') << ss.str() << std::endl;
    }
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
