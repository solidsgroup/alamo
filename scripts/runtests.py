#!/usr/bin/env python3
import argparse
import os, glob, subprocess
import configparser, io

class color:
    reset = "\033[0m"
    red   = "\033[31m"
    green   = "\033[32m"
    boldgreen   = "\033[1m\033[32m"
    boldyellow   = "\033[1m\033[33m"
    bold = "\033[1m"
    lightgray = "\033[37m"
    darkgray = "\033[90m"

parser = argparse.ArgumentParser(description='Configure ALAMO');
parser.add_argument('tests', default=None, nargs='*', help='Spatial dimension [3]')
parser.add_argument('--serial',action='store_true',default=False,help='Run in serial only (no mpi)')
parser.add_argument('--dim',default=None,type=int,help='Specify dimensions to run in')
args=parser.parse_args()

def test(testdir):
    fails = 0
    skips = 0
    tests = 0
    checks = 0

    cfgfile = io.StringIO()
    input = open(testdir + "/input")
    for line in input.readlines():
        if line.startswith("#@"):
            cfgfile.write(line.replace("#@",""))
    cfgfile.seek(0)
    config = configparser.ConfigParser()
    config.read_file(cfgfile)

    if not len(config.sections()):
        print("{}IGNORE {}{}".format(color.darkgray,testdir,color.reset))
        return 0,0,0,0
    print("RUN    {}{}{}".format(color.bold,testdir,color.reset))
    for desc in config.sections():

        check = True

        if 'check' in config[desc].keys(): 
            if config[desc]['check'] in {"no","No","false","False","0"}:
                check = False
            elif config[desc]['check'] in {"yes","Yes","true","True","1"}:
                check = True
            else:
                raise(Exception("Invalid value for check: {}".format(config[desc]['check'])))

        command = ""
        if 'cmd' in config[desc].keys():
            command = config[desc]['cmd']
            if len(config[desc].keys()) > 1:
                raise Exception("If 'cmd' is specified no other parameters can be set. Received " + ",".join(config[desc].keys))
        else:
            dim = 3
            if 'dim' in config[desc].keys(): dim = int(config[desc]['dim'])
            nprocs = 1
            if 'nprocs' in config[desc].keys(): nprocs = int(config[desc]['nprocs'])
            cmdargs = ""
            if 'args' in config[desc].keys(): cmdargs = config[desc]['args']

            if nprocs > 1 and args.serial: 
                continue
            if nprocs > 1: command += "mpirun -np {} ".format(nprocs)
            exe = "./bin/alamo-{}d-g++".format(dim)
            if args.dim and not args.dim == dim:
                continue
            if not os.path.isfile(exe):
                print("  ├ {}{} (skipped - no executable){}".format(color.boldyellow,desc,color.reset))
                skips += 1
                continue
            command += exe + " "
            command += "{}/input ".format(testdir)
            command += cmdargs

        
        print("  ├ " + desc)
        #print("  │      " + command)
        print("  │      Running test............................................",end="",flush=True)
        try:
            p = subprocess.check_output(command.split(),stderr=subprocess.PIPE)
            print("[{}PASS{}]".format(color.boldgreen,color.reset))
            tests += 1
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

        if check:
            print("  │      Checking result.........................................",end="",flush=True)
            try:
                p = subprocess.check_output(["./test"],cwd=testdir,stderr=subprocess.PIPE)
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
        #for line in p.stdout.readlines(): print(line.decode('ascii'),end="")
        #for line in p.stderr.readlines(): print(line.decode('ascii'),end="")
    
    if fails: print("  └ {}{} tests failed{}".format(color.red,fails,color.reset),end="")
    else: print("  └ {}{} tests failed{}".format(color.boldgreen,0,color.reset),end="")
    if skips: print(", {}{} tests skipped{}".format(color.boldyellow,skips,color.reset))
    else: print("")
    return fails, checks, tests, skips







if args.tests: tests = sorted(args.tests)
else: tests = sorted(glob.glob("./tests/*"))

class stats:
    fails = 0
    skips = 0
    checks = 0
    tests = 0

for testdir in tests:
    if (not os.path.isdir(testdir)) or (not os.path.isfile(testdir + "/input")):
        print("{}IGNORE {}{}".format(color.darkgray,testdir,color.reset))
        continue
    f, c, t, s = test(testdir)
    stats.fails += f
    stats.tests += t
    stats.checks += c
    stats.skips += s
    

print("\nTest Summary")
print("{}{} tests run{}".format(color.green,stats.tests,color.reset))
print("{}{} tests run and verified{}".format(color.boldgreen,stats.checks,color.reset))
if not stats.fails: print("{}0 tests failed{}".format(color.boldgreen,color.reset))
else:         print("{}{} tests failed{}".format(color.red,stats.fails,color.reset))
if stats.skips: print("{}{} tests skipped{}".format(color.boldyellow,stats.skips,color.reset))
print("")

exit(stats.fails + stats.skips)
