#!/usr/bin/env python3
import os, glob, subprocess

class color:
    reset = "\033[0m"
    red   = "\033[31m"
    green   = "\033[32m"
    boldgreen   = "\033[1m\033[32m"
    bold = "\033[1m"


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

    print("{}{}{}".format(color.bold,testdir,color.reset))
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






fails = 0
for testdir in glob.glob("./tests/*"):
    if not os.path.isdir(testdir): continue
    if not os.path.isfile(testdir + "/input"): continue
    if not os.path.isfile(testdir + "/test"): continue
    fails += test(testdir)




