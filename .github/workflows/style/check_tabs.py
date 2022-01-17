#!/usr/bin/env python3
import os, sys

offenses = 0
for dirname, subdirlist, filelist in os.walk("src"):
    for f in filelist:
        if not (f.endswith(".H") or f.endswith(".cpp") or f.endswith(".cc") or f.endswith(".py")): continue
        #print(dirname,f)
        lines = open(dirname + "/" + f).readlines()
        for i in range(len(lines)):
            if '\t' in lines[i]:
                print("Tab found in " + dirname + "/" + f + ":" + str(i+1).zfill(4)+"" + lines[i].replace("\t", '\033[41m    \033[0m'),end="")
                offenses += 1

if offenses:
    print("")
    print("Failed: there were " + str(offenses) + " offending tabs")
    print("")
    print("This error was thrown because your code is formatted using tabs instead of spaces.")
    print("To fix this error, replace each tab with 4 spaces.")
    print("You can then execute")
    print("")
    print("      .github/workflows/style/check_tabs.py")
    print("")
    print("to run this script locally to make sure all tabs are gone before your next commit.")
    print("")
    print("If you are using an approved editor, it should automatically format according")
    print("to the rules specified in .editorconfig.")
    print("")
    print("To automatically fix all files, you can install eclint (https://www.npmjs.com/package/eclint)")
    print("Once it is instaled, run")
    print("")
    print("      eclint fix src/*")
    print("")
    print("in the source directory to autofix all formatting errors.")
    print("")
else:
    print("Format check OK")

sys.exit(offenses)

