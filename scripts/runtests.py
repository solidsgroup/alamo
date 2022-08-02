#!/usr/bin/env python3
import argparse
import os, glob, subprocess
import configparser, io
from collections import OrderedDict
from datetime import datetime
import socket
import time
import re

from sympy import capture

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

#
# Get a unique string ID to label all output files
#
now = datetime.now()
testid = now.strftime("output_%Y-%m-%d_%H.%M.%S_"+socket.gethostname())
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
parser.add_argument('--debug',default=False,action='store_true',help='Use the debug version of the code')
parser.add_argument('--profile',default=False,action='store_true',help='Use the profiling version of the code')
parser.add_argument('--benchmark',default=socket.gethostname(),help='Current platform if testing performance')
args=parser.parse_args()

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
    skips = 0
    tests = 0
    checks = 0
    fasters = 0
    slowers = 0

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

    sections = config.sections()
    if args.sections:
        if len(args.sections)>0:
            sections = list(set(config.sections()).intersection(set(args.sections)))

    # If there are no runs specified, then we will not test anything
    # in this directory, and will print an "ignore" message.
    # (Eventually, everything in ./tests should be tested!)
    if not len(sections):
        print("{}IGNORE {}{}".format(color.darkgray,testdir,color.reset))
        return 0,0,0,0,0,0
        
    # Otherwise let the user know that we are in this directory
    print("RUN    {}{}{}".format(color.bold,testdir,color.reset))

    # Iterate through all test configurations
    for desc in sections:

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

        dobenchmark = False
        benchmark = None
        if "benchmark-{}".format(args.benchmark) in config[desc].keys():
            dobenchmark = True
            benchmark = float(config[desc]["benchmark-{}".format(args.benchmark)])

        # Build the command to run the script. This can be done in two ways:
        #
        # 1. with the 'cmd' input where 'cmd' is the precise run command.
        #    This is NOT PORTABLE and useful only for initial testing. You should
        #    generally use option 2.
        # 2. with 'dim', 'nprocs', 'args', etc keywords. See current tests for
        #    examples
        command = ""
        if 'cmd' in config[desc].keys():
            command = config[desc]['cmd']
            if len(config[desc].keys()) > 1:
                raise Exception("If 'cmd' is specified no other parameters can be set. Received " + ",".join(config[desc].keys))
        else:
            dim = 3 # Dimension of alamo to use
            if 'dim' in config[desc].keys(): dim = int(config[desc]['dim'])
            nprocs = 1 # Number of MPI processes, if 1 then will run without mpirun
            if 'nprocs' in config[desc].keys(): nprocs = int(config[desc]['nprocs'])
            cmdargs = "" # Extra arguments to pass to alamo in addition to input file
            if 'args' in config[desc].keys(): cmdargs = config[desc]['args'].replace('\n',' ')

            cmdargs += " plot_file={}/{}_{}".format(testdir,testid,desc)

            if 'ignore' in config[desc].keys():
                cmdargs += " ignore={}".format(config[desc]['ignore'])

            # Quietly ignore this one if running in serial mode.
            if nprocs > 1 and args.serial: 
                continue
            # If not running in serial, specify mpirun command
            if nprocs > 1: command += "mpirun -np {} ".format(nprocs)
            # Specify alamo command.
            
            if args.debug and args.profile: exe = "./bin/alamo-{}d-profile-debug-g++".format(dim)
            elif args.debug: exe = "./bin/alamo-{}d-debug-g++".format(dim)
            elif args.profile: exe = "./bin/alamo-{}d-profile-g++".format(dim)
            else: exe = "./bin/alamo-{}d-g++".format(dim)
            # If we specified a CLI dimension that is different, quietly ignore.
            if args.dim and not args.dim == dim:
                continue
            # If the exe doesn't exist, exit noisily. The script will continue but will return a nonzero
            # exit code.
            if not os.path.isfile(exe):
                print("  ├ {}{} (skipped - no executable){}".format(color.boldyellow,desc,color.reset))
                skips += 1
                continue
            command += exe + " "
            command += "{}/input ".format(testdir)
            command += cmdargs

        
        # Run the actual test.
        print("  ├ " + desc)
        if args.cmd: print("  ├      " + command)
        print("  │      Running test............................................",end="",flush=True)
        # Spawn the process and wait for it to finish before continuing.
        try:
            timeStarted = time.time()
            p = subprocess.run(command.split(),capture_output=True,check=True)
            executionTime = time.time() - timeStarted
            fstdout = open("{}/{}_{}/stdout".format(testdir,testid,desc),"w")
            fstdout.write(ansi_escape.sub('',p.stdout.decode('ascii')))
            fstdout.close()
            fstderr = open("{}/{}_{}/stderr".format(testdir,testid,desc),"w")
            fstderr.write(ansi_escape.sub('',p.stderr.decode('ascii')))
            fstderr.close()
            print("[{}PASS{}]".format(color.boldgreen,color.reset), "({:.2f}s".format(executionTime),end="")
            if dobenchmark:
                if abs(executionTime - benchmark) / (executionTime + benchmark) < 0.01: print(", no change)")
                elif abs(executionTime < benchmark):
                    print(",{} {:.2f}% faster{})".format(color.blue,100*(benchmark-executionTime)/executionTime,color.reset))
                    fasters += 1
                else:
                    print(",{} {:.2f}% slower{})".format(color.magenta,100*(executionTime-benchmark)/executionTime,color.reset))
                    slowers += 1
            else: print(")")
            tests += 1
        # If an error is thrown, we'll go here. We will print stdout and stderr to the screen, but 
        # we will continue with running other tests. (Script will return an error)
        except subprocess.CalledProcessError as e:
            print("[{}FAIL{}]".format(color.red,color.reset))
            print("  │      {}CMD   : {}{}".format(color.red,' '.join(e.cmd),color.reset))
            for line in e.stdout.decode('ascii').split('\n'): print("  │      {}STDOUT: {}{}".format(color.red,line,color.reset))
            for line in e.stderr.decode('ascii').split('\n'): print("  │      {}STDERR: {}{}".format(color.red,line,color.reset))
            fails += 1
            continue
        # Catch-all handling so that if something else odd happens we'll still continue running.
        except Exception as e:
            print("[{}FAIL{}]".format(color.red,color.reset))
            for line in str(e).split('\n'): print("  │      {}{}{}".format(color.red,line,color.reset))
            fails += 1
            continue
        
        # If we have specified that we are doing a check, use the 
        # ./tests/MyTest/test 
        # script to determine if the run was successful.
        # The exception handling is basically the same as for the above test.
        if check:
            print("  │      Checking result.........................................",end="",flush=True)
            try:
                cmd = ["./test","{}_{}".format(testid,desc)]
                if "check-file" in config[desc].keys():
                    cmd.append(config[desc]['check-file'])
                if args.cmd: 
                    print("  ├      " + ' '.join(cmd))
                p = subprocess.check_output(cmd,cwd=testdir,stderr=subprocess.PIPE)
                checks += 1
            except subprocess.CalledProcessError as e:
                print("[{}FAIL{}]".format(color.red,color.reset))
                print("  │      {}CMD   : {}{}".format(color.red,' '.join(e.cmd),color.reset))
                for line in e.stdout.decode('ascii').split('\n'): print("  │      {}STDOUT: {}{}".format(color.red,line,color.reset))
                for line in e.stderr.decode('ascii').split('\n'): print("  │      {}STDERR: {}{}".format(color.red,line,color.reset))
                fails += 1
                continue
            except Exception as e:
                print("[{}FAIL{}]".format(color.red,color.reset))
                for line in str(e).split('\n'): print("  │      {}{}{}".format(color.red,line,color.reset))
                fails += 1
                continue

            print("[{}PASS{}]".format(color.boldgreen,color.reset))
    
    # Print a quick summary for this test family.
    if fails: print("  └ {}{} tests failed{}".format(color.red,fails,color.reset),end="")
    else: print("  └ {}{} tests failed{}".format(color.boldgreen,0,color.reset),end="")
    if skips: print(", {}{} tests skipped{}".format(color.boldyellow,skips,color.reset))
    else: print("")
    return fails, checks, tests, skips, fasters, slowers

# We may wish to pass in specific test directories. If we do, then test those only.
# Otherwise look at everything in ./tests/
if args.tests: tests = sorted(args.tests)
else: tests = sorted(glob.glob("./tests/*"))

class stats:
    fails = 0   # Number of failed runs - script errors if this is nonzero
    skips = 0   # Number of tests that were unexpectedly skipped - script errors if this is nonzero
    checks = 0  # Number of successfully passed checks
    tests = 0   # Number of successful checks
    fasters = 0
    slowers = 0

# Iterate through all test directories, running the above "test" function
# for each.
for testdir in tests:
    if (not os.path.isdir(testdir)) or (not os.path.isfile(testdir + "/input")):
        print("{}IGNORE {} (no input){}".format(color.darkgray,testdir,color.reset))
        continue
    f, c, t, s, fa, sl = test(testdir)
    stats.fails += f
    stats.tests += t
    stats.checks += c
    stats.skips += s
    stats.fasters += fa
    stats.slowers += sl
    

# Print a quick summary of all tests
print("\nTest Summary")
print("{}{} tests run{}".format(color.green,stats.tests,color.reset))
print("{}{} tests run and verified{}".format(color.boldgreen,stats.checks,color.reset))
if not stats.fails: print("{}0 tests failed{}".format(color.boldgreen,color.reset))
else:         print("{}{} tests failed{}".format(color.red,stats.fails,color.reset))
if stats.skips: print("{}{} tests skipped{}".format(color.boldyellow,stats.skips,color.reset))
if stats.fasters: print("{}{} tests ran faster".format(color.blue,stats.fasters,color.reset))
if stats.slowers: print("{}{} tests ran slower".format(color.magenta,stats.slowers,color.reset))
print("")

# Return nonzero only if no tests failed or were unexpectedly skipped
exit(stats.fails + stats.skips)
