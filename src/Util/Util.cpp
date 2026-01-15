#include "Util.H"
#include "AMReX_Config.H"
#include "AMReX_DistributionMapping.H"
#include "AMReX_VisMF.H"
#include "Color.H"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"

#include "IO/ParmParse.H"
#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
#include "Color.H"
#include "Numeric/Stencil.H"
#include "Util/MPI.H"
#include "mpi.h"

namespace Util
{

std::string filename = "";
std::string globalprefix = "";
std::pair<std::string,std::string> file_overwrite;
bool initialized = false;
bool finalized = false;


std::map<std::string,int>  Compare_map;
int Compare_cnt = 0;

void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &a_mf, std::string desc,
    amrex::Box domain)
{
#ifndef AMREX_DEBUG
    return;
#endif

    Compare_cnt++;

    Set::Scalar tolerance = 1E-10;

    std::ostringstream messagestream;
    std::string filesan = file;
    Util::String::ReplaceAll(filesan,"/","_");
    messagestream << filesan << "_" << func << "_" << line << desc;
    std::string key = messagestream.str();

    int cntr = 0;
    if (Compare_map.find(key) == Compare_map.end())
        Compare_map[key] = 0;
    else
        Compare_map[key] ++;
    cntr = Compare_map[key];

    messagestream << "_" << cntr;

    std::string name = messagestream.str();
    std::vector<std::string> allnames;
    Util::MPI::AllGather(name,allnames);
    for (auto & othername : allnames)
    {
        if (name != othername)
            Util::ParallelAbort(INFO, "MPI Paths have diverged: I am trying to read ", 
                                name, " but some other process is trying to read ", 
                                othername);
    }

    amrex::BoxArray ba = a_mf.boxArray();
    amrex::DistributionMapping dm = a_mf.DistributionMap();
    int nghost = a_mf.nGrow();
    int ncomp = a_mf.nComp();
    if (domain != amrex::Box::TheUnitBox())
    {
        for (auto &b : ba.boxList()) b = b.grow(2) & domain;
        nghost = 0;
    }
    amrex::MultiFab mf(ba,dm,ncomp,nghost);
    amrex::MultiFab::Copy(mf,a_mf,0,0,ncomp,0);
    if (amrex::ParallelDescriptor::NProcs() == 1)
    {
        Util::ParallelMessage(INFO,"[",Compare_cnt,"] ","Storing ",name);
        Util::ParallelMessage(INFO,"[",Compare_cnt,"] ",name, " - ncomp ",ncomp);
        std::filesystem::create_directory("checks");
        amrex::VisMF::Write(mf, "checks/" + name);
    }
    else
    {
        Util::ParallelMessage(INFO,"[",Compare_cnt,"] ","Comparing ",name);
        {
            amrex::MultiFab mftmp;
            amrex::VisMF::Read(mftmp,"checks/" + name);
            if (ba != mftmp.boxArray())
            {
                Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","our boxarray ",ba);
                Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","saved boxarray ",mftmp.boxArray());
                Util::Abort(file,func,line,"[",Compare_cnt,"] ","different box arrays !! !! !!");
            }
        }

        amrex::MultiFab mforig(ba,mf.distributionMap,ncomp,nghost);
        amrex::VisMF::Read(mforig,"checks/" + name);

        amrex::MultiFab mfdiff(ba, mf.distributionMap, ncomp, nghost);
        amrex::MultiFab::Copy(mfdiff, mf,0,0,ncomp,nghost);
        amrex::MultiFab::Subtract(mfdiff,mforig,0,0,ncomp,nghost);
        mfdiff.abs(0,ncomp,nghost);

        for (int n = 0 ; n < ncomp; n++)
        {
            Set::Scalar max = mfdiff.max(n,nghost);
            auto argmax = mfdiff.maxIndex(n,nghost);
            if (max > tolerance)
            {
                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(INFO,"ours (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Probe(INFO,mf,argmax[0],argmax[1],0,n,3);

                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(INFO,"original (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Probe(INFO,mforig,argmax[0],argmax[1],0,n,3);

                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(INFO,"diff (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Probe(INFO,mfdiff,argmax[0],argmax[1],0,n,3);

                messagestream << "_diff";
                amrex::VisMF::Write(mfdiff,messagestream.str());
                Util::Warning(INFO,"[",Compare_cnt,"] ",name);
                Util::Abort(INFO,"[",Compare_cnt,"] ",max," at ", argmax);
            }
        }
    }
    //if (Compare_cnt == 169) Util::ParallelAbort(INFO,"exiting");
}


std::string GetFileName()
{
    if (filename == "")
    {
        IO::ParmParse pp;
        IO::ParmParse pp_amr("amr");

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
            pp_query("plot_file", filename); // Name of directory containing all output data
        }
        IO::FileNameParse(filename);
        // else
        // if (amrex::ParallelDescriptor::IOProcessor())
        // Util::Abort("No plot file specified! (Specify plot_file = \"plot_file_name\" in input file");
    }
    return filename;
}
void CopyFileToOutputDir(std::string a_path, bool fullpath, std::string prefix)
{
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

    amrex::Initialize(argc, argv);

    IO::ParmParse pp_amrex("amrex");
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

    IO::ParmParse pp;
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
        
        if (pp.contains("prob_lo"))
        {
            std::vector<Set::Scalar> prob_lo, prob_hi;
            // Location of the lower+left+bottom corner
            pp.queryarr("prob_lo", prob_lo, Unit::Length());
            // Location of the upper_right_top corner
            pp.queryarr("prob_hi", prob_hi, Unit::Length());
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
    }
}

void Finalize()
{
    std::string filename = GetFileName();
    if (filename != "")
        IO::WriteMetaData(filename,IO::Status::Complete);
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
