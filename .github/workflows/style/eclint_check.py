#!/usr/bin/env python3
import os, sys, subprocess

offenses = 0
for dirname, subdirlist, filelist in os.walk("src"):
    for f in filelist:
        if not (f.endswith(".H") or f.endswith(".cpp") or f.endswith(".cc") or f.endswith(".py")): continue
        
        p = subprocess.run(['eclint','check',dirname + "/" + f], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.stdout: print(p.stdout)
        if p.stderr: print(p.stderr)
        if p.returncode:
            offenses += 1
        else:
            print(dirname + "/" + f, " OK")

if offenses:
    print("")
    print("Failed: there were " + str(offenses) + " offending files")
    print("")
    print("This error was thrown because your code does not comply to the editorconfig file.")
    print("To fix this error, use an editorconfig-compliant editor to autoformat your file.")
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

