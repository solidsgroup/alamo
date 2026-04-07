#!/usr/bin/env python3
import sys

sys.path.append('./scripts')
from report import init_html, append_html, finalize_html, write_master_html

import argparse
import os, glob, shutil, subprocess
import configparser, io
from collections import OrderedDict
from datetime import datetime
import socket
import time
import re
import pathlib
import threading
import signal
import fnmatch

from sympy import capture

from concurrent.futures import ThreadPoolExecutor, TimeoutError


#
# Code snippit to suppress the "^C" when user presses control C
#
import sys
import termios
import atexit
def disable_ctrl_echo():
    if not sys.stdin.isatty():
        return None
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    new = termios.tcgetattr(fd)
    # Turn off echoing of control characters like ^C
    new[3] = new[3] & ~termios.ECHOCTL
    termios.tcsetattr(fd, termios.TCSADRAIN, new)
    return fd, old
def restore_ctrl_echo(state):
    if state is None:
        return
    fd, old = state
    termios.tcsetattr(fd, termios.TCSADRAIN, old)
# Disable at start
_term_state = disable_ctrl_echo()
# Ensure it gets restored no matter what
atexit.register(lambda: restore_ctrl_echo(_term_state))


class color:
    reset = "\033[0m"
    red   = "\033[31m"
    green   = "\033[32m"
    blue   = "\033[34m"
    magenta   = "\033[35m"
    boldblue   = "\033[1m\033[34m"
    boldgreen   = "\033[1m\033[32m"
    boldyellow   = "\033[1m\033[33m"
    bold = "\033[1m"
    lightgray = "\033[37m"
    darkgray = "\033[90m"

#
# RE tool to strip out color escapes
# 
ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
def clean(text,max_length=60):
    ret = ansi_escape.sub('',text)
    if len(ret) > max_length:
        return ret[:max_length-3] + "..."
    return ret

#
# Get a unique string ID to label all output files
#
now = datetime.now()
testid = now.strftime("output_%Y-%m-%d_%H.%M.%S_"+socket.gethostname())
os.makedirs("report",exist_ok=True)
init_html(html_path=f"report/{testid}.html")
print("Test ID = ",testid)

#
# Special order from SO - dictionary allows for keys to be specified multiple times and 
# config parser will read it all.
class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super().__setitem__(key, value)

#
# Convert metadata file to dictionary
#
def readMetadata(metadatafile):
    metadata = dict()
    for line in metadatafile.readlines():
        if line.startswith('#'): continue;
        if '::' in line:
            ### skip these for now ...
            continue
            #line = re.sub(r'\([^)]*\)', '',line)
            #line = line.replace(" :: ", " = ").replace('[','').replace(',','').replace(']','').replace(' ','')
        if len(line.split(' = ')) != 2: continue;
        col = line.split(' = ')[0]#.replace('.','_')
        val = line.split(' = ')[1].replace('\n','')#.replace('  ','').replace('\n','').replace(';','')
        metadata[col] = val
    return metadata

#
# Provide some simple command line arguments for filtering the types of 
# regression tests to run. For instance, if you need to run without
# mpirun, you can use the --serial flag.
#
parser = argparse.ArgumentParser(description='Configure ALAMO');
parser.add_argument('tests', default=None, nargs='*', help='Spatial dimension [3]')
parser.add_argument('--serial',action='store_true',default=False,help='Run in serial only (no mpi)')
parser.add_argument('--dim',default=None,type=int,help='Specify dimensions to run in')
parser.add_argument('--cmd',default=False,action='store_true',help="Print out the exact command used to run each test")
parser.add_argument('--sections',default=None, nargs='*', help='Specific sub-tests to run')
parser.add_argument('--exe',default=None, nargs='*', help='Run only certain executables')
parser.add_argument('--debug',default=False,action='store_true',help='Use the debug version of the code')
parser.add_argument('--profile',default=False,action='store_true',help='Use the profiling version of the code')
parser.add_argument('--perf',default=False,action='store_true',help='Use clang perf')
parser.add_argument('--coverage',default=False,action='store_true',help='Use the gcov version of the code for all tests')
parser.add_argument('--only-coverage',default=False,action='store_true',help='Gracefully skip non-coverage tests')
parser.add_argument('--only-non-coverage',default=False,action='store_true',help='Gracefully skip coverage tests')
parser.add_argument('--no-coverage',default=False,action='store_true',help='Prevent coverage version of the code from being used')
parser.add_argument('--memcheck',default=None,help='Run tests with memory checking. ')
parser.add_argument('--benchmark',default=socket.gethostname(),help='Current platform if testing performance')
parser.add_argument('--dryrun',default=False,action='store_true',help='Do not actually run tests, just list what will be run')
parser.add_argument('--comp', default="g++", help='Compiler. Options: [g++], clang++, icc')
parser.add_argument('--timeout', default=10000, help='Timeout value in seconds (default: 10000)')
parser.add_argument('--post', default=None, help='Use a post script to post results')
parser.add_argument('--clean', dest='clean', default=True, action='store_true', help='Clean up output files if test is successful (on by default)')
parser.add_argument('--no-clean', dest='clean', default=False, action='store_false', help='Keep all output files')
parser.add_argument('--permissive', dest='permissive', default=False, action='store_true', help='Option to run without erroring out (if at all possible)')
parser.add_argument('--permit-timeout', dest='permit_timeout', default=False, action='store_true', help='Permit timeouts without failing')
parser.add_argument('--permit-skips', dest='permit_skips', default=False, action='store_true', help='Permit skips without failing')
parser.add_argument('--no-backspace',default=False,dest="no_backspace",action='store_true',help="Avoid using backspace (For GH actions)")
parser.add_argument('--check-mpi',default=False,dest="check_mpi",action='store_true',help="Check if MPI is running correctly")
parser.add_argument('--mpirun-flags',dest="mpirun_flags",default="",help="Extra arguments to pass to mpirun (like --oversubscribe). All arguments must be in a string.")
parser.add_argument('--fft',dest="fft",default=False,action='store_true',help="Enable fft-based tests")
parser.add_argument('--fft-only',dest="fft_only",default=False,action='store_true',help="Run fft tests only")
parser.add_argument('--post-timeout', dest="post_timeout", default=10000, help='How long to wait before skipping results posting')
parser.add_argument('--python', default=False,action='store_true', help='Include python tests')
parser.add_argument('--only-python', default=False,action='store_true', help='Run python tests only')
args=parser.parse_args()

if args.coverage and args.no_coverage:
    raise Exception("Cannot specify both --coverage and --no-coverage")
if args.only_coverage and args.no_coverage:
    raise Exception("Cannot specify both --only-coverage and --no-coverage")
if args.memcheck and not args.debug:
    raise Exception("Debug must be enabled with memory check")
if args.memcheck and not args.serial:
    raise Exception("Memory check supported in serial only")

if args.only_python:
    args.python = True

if args.post:
    if not os.path.isfile(args.post):
        raise Exception(args.post,"is not a file")
    sys.path.append(str(pathlib.Path(args.post).parent))
    import post
    postdata = post.init()


class DryRunException(Exception):
    pass

def test(testdir):
    """
    Run tests using alamo on the input file stored in "testdir"
    We assume the following convention: for a test called "MyTest",

    test directory:      ./tests/MyTest          # Directory containig all test information
    alamo input file:    ./tests/MyTest/input    # Input file with config file comments 
    output file:         ./tests/MyTest/output   # All alamo output
    check script:        ./tests/test            # Executable script that returns 0 for successful test
    
    """

    # Counters to track failed tests, skipped tests, passed tests, passed checks
    fails = 0
    kills = 0
    skips = 0
    tests = 0
    checks = 0
    warnings = 0
    fasters = 0
    slowers = 0
    timeouts = 0
    post_fails = 0
    records = []


    
    # Formatting strings
    bs = "\b\b\b\b\b\b"


    if os.path.isfile(testdir + "/input.py"):
        print(f"RUN    {color.bold}{testdir}{color.reset}")
        cmd = f'python {testdir}/input.py'
        print("  │      Running test............................................",end="",flush=True)
        try:
            proc = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate(timeout=int(args.timeout))
            retcode = proc.returncode
            if retcode: raise subprocess.CalledProcessError(retcode, proc.args, output=stdout, stderr=stderr)
            print("[{}PASS{}]".format(color.boldgreen,color.reset))
            tests += 1
            print("  └ Python test complete ")
            return fails, kills, checks, warnings, tests, skips, fasters, slowers, timeouts, records

        except subprocess.CalledProcessError as e:
            print(bs+"[{}FAIL{}]".format(color.red,color.reset))
            print("  │      {}CMD   : {}{}".format(color.red,' '.join(e.cmd),color.reset))
            for line in e.stdout.decode('utf-8').split('\n'): print("  │      {}STDOUT: {}{}".format(color.red,clean(line,1000),color.reset))
            for line in e.stderr.decode('utf-8').split('\n'): print("  │      {}STDERR: {}{}".format(color.red,clean(line,1000),color.reset))
            fails += 1

        if os.path.isfile(f"{testdir}/input"):
            raise(Exception(f"Test {testdir} cannot have both input and input.py"))

        return fails, kills, checks, warnings, tests, skips, fasters, slowers, timeouts, records


    # Parse the input file ./tests/MyTest/input containing #@ comments.
    # Everything commeneted with #@ will be interpreted as a "config" file
    # The variable "config" is a dict of dicts where each item corresponds to
    # a test configuration.
    cfgfile = io.StringIO()
    input = open(testdir + "/input")
    for line in input.readlines():
        if line.startswith("#@"):
            cfgfile.write(line.replace("#@",""))
    cfgfile.seek(0)
    config = configparser.ConfigParser(dict_type=MultiOrderedDict,strict=False)
    config.read_file(cfgfile)

    #
    # Apply wildcard sections and then remove them from the list
    #
    sections = [s for s in config.sections() if "*" not in s]
    section_wildcards = [s for s in config.sections() if "*" in s]
    for desc in sections:
        desc_wildcards = [s for s in section_wildcards if fnmatch.fnmatch(desc,s)]
        for wc in desc_wildcards:
            new = dict(config[wc])
            newargs = None
            if "args" in config[desc].keys() and "args" in config[wc].keys():
                newargs = new["args"] + "\n" + config[desc]["args"]
            new.update(config[desc])
            if newargs: new["args"] = newargs
            config[desc] = new
    for wc in section_wildcards:
        config.remove_section(wc)

    if args.sections:
        if len(args.sections)>0:
            sections = list(set(config.sections()).intersection(set(args.sections)))

    # If there are no runs specified, then we will not test anything
    # in this directory, and will print an "ignore" message.
    # (Eventually, everything in ./tests should be tested!)
    if not len(sections):
        print("{}IGNORE {}{}".format(color.darkgray,testdir,color.reset))
        return 0,0,0,0,0,0,0,0,0,[]
        
    # Otherwise let the user know that we are in this directory
    print("RUN    {}{}{}".format(color.bold,testdir,color.reset))

    # Store a list of all the cleanup commands that will be run after
    # going through all sections
    cleanup_paths = []

    # Iterate through all test configurations
    for desc in sections:

        # Start off with the current environment as the default
        env = os.environ.copy()

        record = dict()
        record['testdir'] = testdir
        record['section'] = desc
        record['test-section'] = record['testdir']+'/'+record['section']
        record['section'] = desc
        record['testid']  = testid
        record['path'] = "{}/{}_{}".format(testdir,testid,desc)


        p = subprocess.check_output('git rev-parse --abbrev-ref HEAD'.split(),stderr=subprocess.PIPE)
        record['branch'] = p.decode('utf-8').replace('\n','')

        # In some cases we want to run the exe but can't check it.
        # Skipping the check can be done by specifying the "check" input.
        check = True
        if 'check' in config[desc].keys(): 
            if config[desc]['check'] in {"no","No","false","False","0"}:
                check = False
            elif config[desc]['check'] in {"yes","Yes","true","True","1"}:
                check = True
            else:
                raise(Exception("Invalid value for check: {}".format(config[desc]['check'])))
            config[desc].pop('check')

        # Determine if we want to use the coverage version of the code.
        coverage = args.coverage
        if 'coverage' in config[desc].keys():
            if config[desc]['coverage'] in {"no","No","false","False","0"}:
                coverage = False
            elif config[desc]['coverage'] in {"yes","Yes","true","True","1"}:
                coverage = True
            else:
                raise(Exception("Invalid value for coverage: {}".format(config[desc]['coverage'])))
            config[desc].pop('coverage')

        if args.only_coverage and not coverage:
            continue
        if args.only_non_coverage and coverage:
            continue
        if args.no_coverage:
            coverage = False

        timeout = int(args.timeout)
        if 'timeout' in config[desc].keys():
            timeout = int(config[desc]['timeout'])
            config[desc].pop('timeout')

        dobenchmark = False
        benchmark = None
        benchmark_name = "benchmark-{}".format(args.benchmark)
        if benchmark_name in config[desc].keys():
            dobenchmark = True
            benchmark = float(config[desc]["benchmark-{}".format(args.benchmark)])
            config[desc].pop(benchmark_name)

        # Build the command to run the script. This can be done in two ways:
        #
        # 1. with the 'cmd' input where 'cmd' is the precise run command.
        #    This is NOT PORTABLE and useful only for initial testing. You should
        #    generally use option 2.
        # 2. with 'dim', 'nprocs', 'args', etc keywords. See current tests for
        #    examples

        command = "" # The executable that gets called
        cmdargs = "" # Extra arguments to pass to alamo in addition to input file
        if 'cmd' in config[desc].keys():
            command = config[desc]['cmd']
            if len(config[desc].keys()) > 1:
                raise Exception("If 'cmd' is specified no other parameters can be set. Received " + ",".join(config[desc].keys))
            config[desc].pop('cmd')
        else:
            # If we are doing memory checking, cut off the simulation early
            if args.memcheck: cmdargs += " max_step=2 "

            exe = 'alamo'
            if 'exe' in config[desc].keys():
                exe = config[desc]['exe']
                config[desc].pop('exe')
            dim = 3 # Dimension of alamo to use
            if 'dim' in config[desc].keys():
                dim = int(config[desc]['dim'])
                config[desc].pop('dim')
            nprocs = 1 # Number of MPI processes, if 1 then will run without mpirun
            if 'nprocs' in config[desc].keys():
                nprocs = int(config[desc]['nprocs'])
                config[desc].pop('nprocs')
            if 'args' in config[desc].keys():
                cmdargs += config[desc]['args'].replace('\n',' ')
                config[desc].pop('args')

            cmdargs += " plot_file={}/{}_{}".format(testdir,testid,desc)

            if 'ignore' in config[desc].keys():
                cmdargs += " ignore={}".format(config[desc]['ignore'])
                config[desc].pop('ignore')

            # Quietly ignore this one if running in serial mode.
            if nprocs > 1 and args.serial: 
                continue

            # Skip all non-fft tests if --fft-only 
            if args.fft_only and 'fft' not in config[desc].keys():
                continue

            # Otherwise, skip fft tests unless passing in -fft flag
            if 'fft' in config[desc].keys():
                if config[desc]['fft'] in {"yes","Yes","true","True","1"}:
                    if not args.fft and not args.fft_only: continue
                config[desc].pop('fft')

            # Specify performance flag
            if args.perf:
                env["CPUPROFILE"] = "profile.prof"

            # If not running in serial, specify mpirun command
            if nprocs > 1: command += f"mpirun {args.mpirun_flags} -np {nprocs} "
            # Specify alamo command.
            
            exestr = "./bin/{}-{}d".format(exe,dim)
            if args.debug: exestr += "-debug"
            if args.memcheck: exestr += "-{}".format(args.memcheck)
            if args.profile: exestr += "-profile"
            if args.perf: exestr += "-perf"
            if coverage: exestr += "-coverage"
            exestr += "-"+args.comp
            
            #
            # Sometimes we don't have a coverage version built; in that case,
            # use the non-coverage version.
            #
            if not os.path.isfile(exestr):
                exestr=exestr.replace("-coverage","")

            # If we specified a CLI dimension that is different, quietly ignore.
            if args.dim and not args.dim == dim:
                continue

            if args.exe:
                if not exe in args.exe:
                    continue
                
            # Determine if we want to restart this simulation from a previous test.
            if 'restart' in config[desc].keys():
                restartfile = config[desc]['restart'].replace(r'{testid}',testid)
                if not os.path.isdir(restartfile):
                    print("  ├ {}{} (skipped - no restart file){}".format(color.boldyellow,desc,color.reset))
                    skips += 1
                    continue
                cmdargs += " restart=" + restartfile
                config[desc].pop('restart')
        

            # If the exestr doesn't exist, exit noisily.
            # The script will continue but will return a nonzero
            # exit code.
            if not os.path.isfile(exestr):
                print("  ├ {}{} (skipped - no {} executable){}".format(color.boldyellow,desc,exestr,color.reset))
                skips += 1
                continue

            # If we have enabled "skipping," exit noisily. The script will continue but will return a nonzero
            # exit code.
            if 'skip' in config[desc].keys():
                if config[desc]['skip'].lower() in ['true','yes','1']:
                    print("  ├ {}{} (skip indicated in input){}".format(color.boldyellow,desc,color.reset))
                    skips += 1
                    continue
                config[desc].pop('skip')

            command += exestr + " "
            command += "{}/input ".format(testdir)
            command += cmdargs
            record['cmd_test'] = command
        


        # 
        # [ R U N   T H E   T E S T ]
        #

        print("  ├ " + desc)
        if args.cmd: print("  ├      " + command)
        if args.no_backspace:
            print("  │      Running test............................................",end="")
            bs = ""
        else:
            print("  │      Running test............................................[----]",end="",flush=True)
            
        # Spawn the process and wait for it to finish before continuing.
        try:
            if args.dryrun: raise DryRunException()
            timeStarted = time.time()

            #
            # This is a thread that periodically checks the metadata file to determine
            # progress of the alamo run. It scans the metadata file for the alamo-computed
            # progress and prints out the current progress to the terminal.
            #
            def check_progress(proc,none):
                n = 0.0
                while True:
                    if proc.poll() is not None:
                        break
                    try: 
                        metadatafile = open("{}/{}_{}/metadata".format(testdir,testid,desc),"r")
                        metadata = readMetadata(metadatafile)
                        metadatafile.close()
                        status = int(metadata["Status"].split('(')[1].split('%')[0])
                        if not args.no_backspace:
                            print(bs + f"[{status:3d}%]",end="",flush=True)
                    except Exception as e:
                        True #do nothing
                    time.sleep(1)

            def write_log(stdout,stderr):
                try:
                    if stdout is None: stdout = ""
                    if isinstance(stdout,bytes): stdout = stdout.decode('utf-8')
                    fstdout = open("{}/{}_{}/stdout".format(testdir,testid,desc),"w")
                    fstdout.write(ansi_escape.sub('',stdout))
                    fstdout.close()
                except Exception:
                    pass
                try:
                    if stderr is None: stderr = ""
                    if isinstance(stderr,bytes): stderr = stderr.decode('utf-8')
                    fstderr = open("{}/{}_{}/stderr".format(testdir,testid,desc),"w")
                    fstderr.write(ansi_escape.sub('',stderr))
                    fstderr.close()
                except Exception:
                    pass


            # Start the run
            proc = subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                    env=env)

            # Start the check_progress thread
            monitor_thread = threading.Thread(target=check_progress, args=(proc,""))
            monitor_thread.start()

            # Now, wait for test to complete or error out
            stdout, stderr = proc.communicate(timeout=timeout)
            retcode = proc.returncode
            if retcode: raise subprocess.CalledProcessError(retcode, proc.args, output=stdout, stderr=stderr)

            if args.check_mpi:
                pattern = r"MPI initialized with (\d+) MPI processes"
                match = re.search(pattern, stdout.decode('utf-8'))
                if match:
                    actual_nprocs = int(match.group(1))
                    if nprocs != actual_nprocs:
                        raise subprocess.CalledProcessError(retcode, proc.args, output=stdout,
                                                            stderr=f"MPI is not working: should run with {nprocs} procs but is running with {actual_nprocs}!".encode('utf-8'))
                else:
                    raise subprocess.CalledProcessError(retcode, proc.args, output=stdout,
                                                        stderr=f"Could not tell if MPI is working!".encode('utf-8'))

            executionTime = time.time() - timeStarted
            record['executionTime'] = str(executionTime)
            print(bs+"[{}PASS{}]".format(color.boldgreen,color.reset), "({:.2f}s".format(executionTime),end="")
            record['runStatus'] = 'PASS'
            if dobenchmark:
                if abs(executionTime - benchmark) / (executionTime + benchmark) < 0.01: print(", no change)")
                elif abs(executionTime < benchmark):
                    print(",{} {:.2f}% faster{})".format(color.blue,100*(benchmark-executionTime)/executionTime,color.reset))
                    fasters += 1
                else:
                    print(",{} {:.2f}% slower{})".format(color.magenta,100*(executionTime-benchmark)/executionTime,color.reset))
                    slowers += 1
            else: print(")")

            write_log(stdout,stderr)
            tests += 1
        except KeyboardInterrupt:
            print(bs+"[{}KILL{}]".format(color.red,color.reset))
            record['runStatus'] = 'KILL'
            kills += 1
            continue

        # If an error is thrown, we'll go here. We will print stdout and stderr to the screen, but 
        # we will continue with running other tests. (Script will return an error unless permit-timout
        # has been enabled - usually for profiling)
        except subprocess.CalledProcessError as e:
            print(bs+"[{}FAIL{}]".format(color.red,color.reset))
            record['runStatus'] = 'FAIL'
            print("  │      {}CMD   : {}{}".format(color.red,' '.join(e.cmd),color.reset))
            stdout = e.stdout
            if isinstance(stdout,bytes): stdout = stdout.decode('utf-8')
            stderr = e.stderr
            if isinstance(stderr,bytes): stderr = stderr.decode('utf-8')
            for line in stdout.split('\n'): print("  │      {}STDOUT: {}{}".format(color.red,clean(line,1000),color.reset))
            for line in stderr.split('\n'): print("  │      {}STDERR: {}{}".format(color.red,clean(line,1000),color.reset))
            fails += 1
            write_log(stdout,stderr)
            append_html(record)
            continue
        # If we run out of time we go here. We will print stdout and stderr to the screen, but 
        # we will continue with running other tests. (Script will return an error unless --permit-timout is specified)
        except subprocess.TimeoutExpired as e:

            #
            # Attempt to kill the program gracefully ith SIGINT so that it has a chance to
            # write .gcda files (for coverage report)
            #
            if not args.no_backspace:
                print(bs+"[kill]",end="")
            proc.send_signal(signal.SIGINT)
            for i in range(10): # Give it about 20 seconds...
                if proc.poll() is None:
                    time.sleep(2)
                else: break
                    
            #
            # Time's up!
            #
            proc.kill()

            print(bs+"[{}TIME{}]".format(color.lightgray,color.reset))
            record['runStatus'] = 'TIMEOUT'
            try:
                stdout = e.stdout
                if isinstance(stdout,bytes): stdout = stdout.decode('utf-8')
                stderr = e.stderr
                if isinstance(stderr,bytes): stderr = stderr.decode('utf-8')
                write_log(stdout,stderr)
                append_html(record)
                stdoutlines =stdout.split('\n')
                if len(stdoutlines) < 10:
                    for line in stdoutlines:      print("  │      {}STDOUT: {}{}".format(color.lightgray,clean(line),color.reset))
                else:
                    for line in stdoutlines[:5]:  print("  │      {}STDOUT: {}{}".format(color.lightgray,clean(line),color.reset))
                    for i in range(3):            print("  │      {}        {}{}".format(color.lightgray,"............",color.reset))
                    for line in stdoutlines[-5:]: print("  │      {}STDOUT: {}{}".format(color.lightgray,clean(line),color.reset))
            except Exception as e1:
                for line in str(e1).split('\n'):
                    print("  │      {}EXCEPT: {}{}".format(color.red,line,color.reset))
            timeouts += 1
            append_html(record)
            continue
        except DryRunException as e:
            print("")
            record['runStatus'] = '----'
            
        # Catch-all handling so that if something else odd happens we'll still continue running.
        except Exception as e:
            print(bs+"[{}FAIL{}]".format(color.red,color.reset))
            record['runStatus'] = 'FAIL'
            for line in str(e).split('\n'): print("  │      {}{}{}".format(color.red,clean(line,1000),color.reset))
            fails += 1
            append_html(record)
            continue
        
        #
        # 
        # [ P R O F I L I N G ]
        #
        #

        if args.perf:

            profile_file = None

            if nprocs>1:
                profile_file = "profile.prof.rank-0"
            else:
                profile_file = "profile.prof"

            output_svg = f"{testdir}/{testid}_{desc}/flame.svg"

            # Step 1: Collapse the profile
            try:
                print("  │      Processing profiling data...............................",
                      end="",flush=True)
                pprof_cmd = "google-pprof" if shutil.which("google-pprof") else "pprof"
                result = subprocess.run(
                    [pprof_cmd, "--collapsed", exestr, profile_file],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=True, text=True)

                print("[{}DONE{}]".format(color.boldgreen,color.reset))
                print("  │      Generating flame plot...................................",
                      end="",flush=True)

                with open("out.folded", "w") as f:
                    f.write(result.stdout)
            
                result = subprocess.run(
                    ["ext/brendangregg/FlameGraph/flamegraph.pl", "out.folded"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=True, text=True )
            
                with open(output_svg, "w") as f:
                    f.write(result.stdout)
            
                print("[{}DONE{}]".format(color.boldgreen,color.reset))
                print(f"  │      {color.lightgray}{output_svg}{color.reset}")

            except subprocess.CalledProcessError as e:
                print("[{}FAIL{}]".format(color.red,color.reset))
                for line in e.stderr.split('\n'):
                    print("  │      {}ERROR: {}{}".format(color.red,clean(line,1000),color.reset))
                print("Error during flamegraph generation:")
                print(e.stderr)
                raise


            # remove profiling files
            os.remove("out.folded")
            for f in glob.glob("profile.prof*"):
                os.remove(f)


        #
        # 
        # [ V E R I F I C A T I O N   T E S T ]
        #
        #


        # If we have specified that we are doing a check, use the 
        # ./tests/MyTest/test 
        # script to determine if the run was successful.
        # The exception handling is basically the same as for the above test.
        if check and not args.memcheck:
            try:
                cmd = ["./test","{}_{}".format(testid,desc)]
                if args.cmd: 
                    print("  ├      " + ' '.join(cmd))
                print("  │      Checking result.........................................",end="",flush=True)

                if args.dryrun:
                    if ('check-file') in config[desc].keys(): config[desc].pop('check-file')
                    raise DryRunException()
                if "check-file" in config[desc ].keys():
                    cmd.append(config[desc]['check-file'])
                    config[desc].pop('check-file')

                proc = subprocess.Popen(cmd,cwd=testdir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                if retcode == 0:
                    if (stderr):
                        print("[{}WARN{}]".format(color.boldyellow,color.reset)) 
                        for line in stderr.decode('utf-8').split('\n')[:-1]: print("  │      {}STDERR: {}{}".format(
                                color.boldyellow,line,color.reset))
                        record['checkStatus'] = 'WARN'
                        warnings += 1
                    else:
                        print("[{}PASS{}]".format(color.boldgreen,color.reset)) 
                        record['checkStatus'] = 'PASS'
                    checks += 1
                else:
                    raise subprocess.CalledProcessError(retcode, proc.args, output=stdout, stderr=stderr)

            except subprocess.CalledProcessError as e:
                print("[{}FAIL{}]".format(color.red,color.reset))
                record['checkStatus'] = 'FAIL'
                print("  │      {}CMD   : {}{}".format(color.red,' '.join(e.cmd),color.reset))
                for line in e.stdout.decode('utf-8').split('\n'): print("  │      {}STDOUT: {}{}".format(color.red,clean(line,1000),color.reset))
                for line in e.stderr.decode('utf-8').split('\n'): print("  │      {}STDERR: {}{}".format(color.red,clean(line,1000),color.reset))
                fails += 1
                append_html(record)
                continue
            except DryRunException as e:
                print("[----]")
                record['checkStatus'] = "----"
            except Exception as e:
                print("[{}FAIL{}]".format(color.red,color.reset))
                record['checkStatus'] = 'FAIL'
                for line in str(e).split('\n'): print("  │      {}{}{}".format(color.red,line,color.reset))
                fails += 1
                append_html(record)
                continue
        else:
            record['checkStatus'] = 'NONE'


        #
        # Scan metadata file for insteresting things to include in the record
        #
        if args.post:
            print("  │      Posting result..........................................",end="")
            try:
                metadatafile = open("{}/{}_{}/metadata".format(testdir,testid,desc),"r")
                metadata = readMetadata(metadatafile)
                metadatafile.close()
                record['git_commit_hash'] = metadata['Git_commit_hash']
                record['platform'] = metadata['Platform']
                record['test-section'] = record['testdir'] + '/' + record['section']
                p = subprocess.run('git show --no-patch --format=%ci {}'.format(record['git_commit_hash'].split('-')[0]).split(),capture_output=True)
                record['git_commit_date'] = p.stdout.decode('utf-8').replace('\n','')

                with ThreadPoolExecutor(max_workers=1) as executor:
                    future = executor.submit(post.updateDatabase, postdata,[record])
                    future.result(timeout=int(args.post_timeout));
                print("[{}DONE{}]".format(color.boldgreen,color.reset))
            except TimeoutError as e:
                post_fails += 1
                print("[{}TIME{}]".format(color.lightgray,color.reset))
            except Exception as e:
                post_fails += 1
                print("[{}FAIL{}]".format(color.red,color.reset))
                for line in str(e).split('\n'):
                    print("  │      {}{}{}".format(color.red,line,color.reset))
                if not args.permissive:
                    raise Exception()
                True # permissive

        if not args.dryrun:
            records.append(record)
            append_html(record)

        #
        # Clean up all of the node and cell file 
        #
        ok_to_clean = False
        if args.clean and record['runStatus'] == 'PASS':
            if 'checkStatus' not in record.keys():
                ok_to_clean = True # successful run and no testing done
            elif record['checkStatus'] == 'NONE':
                ok_to_clean = True # successful run and no testing done
            elif record['checkStatus'] == 'PASS':
                ok_to_clean = True # successful run and testing completed
            else:
                ok_to_clean = False # successful run but testing failed
        else:
            ok_to_clean = False

        if ok_to_clean:
            path = "{}/{}_{}".format(testdir,testid,desc)
            cleanup_paths.append(path)
            
        
        #
        # CHECK FOR UNUSED INPUTS
        #
        for key in config[desc].keys():
            if "benchmark-" in key: config[desc].pop(key)
            if key == "check-file": config[desc].pop(key)
        if len(config[desc].keys()):
            raise Exception("The following were specified but not used: "+str(list(config[desc].keys())))
    
    if len(cleanup_paths):
        print(f"  ├ {color.lightgray}Cleaning up {len(cleanup_paths)} directories{color.reset}")
        for path in cleanup_paths:
            p = subprocess.run(f'rm -rf *cell *node',capture_output=True,cwd=path,shell=True)
        

            
    # Print a quick summary for this test family.
    summary = "  └ "
    sums = []
    if tests: sums.append("{}{} tests run{}".format(color.blue,tests,color.reset))
    if checks: sums.append("{}{} checks passed{}".format(color.green,checks,color.reset))
    if warnings: sums.append("{}{} warnings{}".format(color.boldyellow,checks,color.reset))
    if fails: sums.append("{}{} tests failed{}".format(color.red,fails,color.reset))
    if kills: sums.append("{}{} tests killed{}".format(color.red,kills,color.reset))
    if skips: sums.append("{}{} tests skipped{}".format(color.boldyellow,skips,color.reset))
    if timeouts: sums.append("{}{} tests timed out{}".format(color.lightgray,timeouts,color.reset))
    if post_fails: sums.append("{}{} tests failed to post{}".format(color.lightgray,post_fails,color.reset))
    print(summary + ", ".join(sums))
    return fails, kills, checks, warnings, tests, skips, fasters, slowers, timeouts, records

# We may wish to pass in specific test directories. If we do, then test those only.
# Otherwise look at everything in ./tests/
if args.tests: tests = sorted(args.tests)
else: tests = sorted(glob.glob("./tests/*"))

tests = [str(pathlib.Path(f)) for f in tests]

class stats:
    fails = 0   # Number of failed runs - script errors if this is nonzero
    kills = 0   # Number of jobs killed with Ctrl+C
    skips = 0   # Number of tests that were unexpectedly skipped - script errors if this is nonzero
    checks = 0  # Number of successfully passed checks
    warnings = 0
    tests = 0   # Number of successful checks
    fasters = 0
    slowers = 0
    timeouts = 0
    records = []

# Iterate through all test directories, running the above "test" function
# for each.
for testdir in tests:
    if (not os.path.isdir(testdir)):
        print(f"{color.darkgray}IGNORE {testdir} (no input){color.reset}")
        continue
    if os.path.isfile(testdir + "/input") and os.path.isfile(testdir + "/input.py"):
        print(f"{color.darkgray}IGNORE {testdir} (cannot specify both input and input.py){color.reset}")
        continue
    if os.path.isfile(testdir + "/input") and args.only_python:
        print(f"{color.darkgray}IGNORE {testdir} (running python tests only){color.reset}")
        continue
    if os.path.isfile(testdir + "/input.py") and not args.python:
        print(f"{color.darkgray}IGNORE {testdir} (did not specify --python or --only-python){color.reset}")
        continue

    f, k, c, w, t, s, fa, sl, to, re = test(testdir)
    stats.fails += f
    stats.kills += k
    stats.tests += t
    stats.checks += c
    stats.warnings += w
    stats.skips += s
    stats.fasters += fa
    stats.slowers += sl
    stats.timeouts += to
    stats.records += re

# Print a quick summary of all tests
print("\nTest Summary")
print("{}{} tests run{}".format(color.blue,stats.tests,color.reset))
print("{}{} tests run and verified{}".format(color.boldgreen,stats.checks,color.reset))
if not stats.fails: print("{}0 tests failed{}".format(color.boldgreen,color.reset))
else:         print("{}{} tests failed{}".format(color.red,stats.fails,color.reset))
if stats.kills:         print("{}{} tests killed{}".format(color.red,stats.kills,color.reset))
if stats.warnings: print("{}{} warnings{}".format(color.boldyellow,stats.warnings,color.reset))
if stats.skips: print("{}{} tests skipped{}".format(color.boldyellow,stats.skips,color.reset))
if stats.fasters: print("{}{} tests ran faster".format(color.blue,stats.fasters,color.reset))
if stats.slowers: print("{}{} tests ran slower".format(color.magenta,stats.slowers,color.reset))
if stats.timeouts: print("{}{} tests timed out".format(color.lightgray,stats.timeouts,color.reset))
print("")

return_code = stats.fails
if not args.permit_skips:
    return_code += stats.skips
if not args.permit_timeout:
    return_code += stats.timeouts

finalize_html()

# Return nonzero only if no tests failed or were unexpectedly skipped
exit(return_code)
