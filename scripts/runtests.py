#!/usr/bin/env python3
import argparse
import os, glob, subprocess

class color:
    reset = "\033[0m"
    red   = "\033[31m"
    green   = "\033[32m"
    boldgreen   = "\033[1m\033[32m"
    boldyellow   = "\033[1m\033[33m"
    bold = "\033[1m"
    lightgray = "\033[37m"

def test(testdir):
    fails = 0

    input = open(testdir + "/input")
    desc = []
    exes = []
    for line in input.readlines():
        if not line.startswith("#@"): continue
        line = line.replace("#@","")
        cmd = line.split("]")[-1]#.replace("\n","")
        while cmd.startswith(" "): cmd = cmd[1:]
        des = line.split("]")[0].split("[")[-1]

        desc.append(des)            
        exes.append(cmd)

    print("RUN  {}{}{}".format(color.bold,testdir,color.reset))
    for description, command in zip(desc,exes):
        print("  ├ " + description)
        print("  │      Running test..................................................[    ]",end="",flush=True)
        p = subprocess.run(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #for line in p.stdout.readlines(): print(line.decode('ascii'),end="")
        if p.returncode:
            print("\b\b\b\b\b\b[{}FAIL{}]".format(color.red,color.reset))
            fails += 1
            continue
        print("\b\b\b\b\b\b[{}PASS{}]".format(color.boldgreen,color.reset))
        print("  │      Checking result...............................................[    ]",end="",flush=True)
        p = subprocess.run(["./test"],cwd=testdir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #p.wait()
        if p.returncode:
            print("\b\b\b\b\b\b[{}FAIL{}]".format(color.red,color.reset))
            fails += 1
            continue
        print("\b\b\b\b\b\b[{}PASS{}]".format(color.boldgreen,color.reset))
        #for line in p.stdout.readlines(): print(line.decode('ascii'),end="")
        #for line in p.stderr.readlines(): print(line.decode('ascii'),end="")
    
    if fails: print("  └ {}{} tests failed{}".format(color.red,fails,color.reset))
    else: print("  └ {}{} tests failed{}".format(color.boldgreen,0,color.reset))
    return fails





parser = argparse.ArgumentParser(description='Configure ALAMO');
parser.add_argument('tests', default=None, nargs='*', help='Spatial dimension [3]')
args=parser.parse_args()


if args.tests: tests = sorted(args.tests)
else: tests = sorted(glob.glob("./tests/*"))

fails = 0
skips = 0
for testdir in tests:
    if (not os.path.isdir(testdir)) or (not os.path.isfile(testdir + "/input")) or (not os.path.isfile(testdir + "/test")):
        print("{}SKIP {}{}".format(color.boldyellow,testdir,color.reset))
        skips += 1
        continue
    fails += test(testdir)

print("\nTest Summary")
if not fails: print("{}0 tests failed{}".format(color.boldgreen,color.reset))
else:         print("{}{} tests failed{}".format(color.red,fails,color.reset))
if skips: print("{}{} tests skipped{}".format(color.boldyellow,skips,color.reset))
print("")

exit(fails)
