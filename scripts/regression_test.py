import os
import re
from glob import glob
import subprocess
from datetime import datetime
import filecmp

timestamp = datetime.today().strftime('%Y%m%d_%H%M%S')

tests = [ { "name" : "TrigTest" } ,
          { "name" : "Eshelby"  } ]

#benchmark = False

for test in tests:
    print("")
    print("Test: ", test['name'])
    test_dir = "tests/" + test['name'] + "/"
    rt_plot_dir = "rt-" + timestamp + "/" + test['name'] + "/"
    bm_plot_dir = "bm-" + timestamp + "/" + test['name'] + "/"

    subprocess.run(["mkdir", "-p", rt_plot_dir])

    print("Running ...... ",end="")
    if os.path.isfile(test_dir + "/input"):
        fstdout = open(rt_plot_dir+"/"+"stdout","w")
        fstderr = open(rt_plot_dir+"/"+"stderr","w")
        ret = subprocess.run(["./bin/alamo-2d-g++", "{}/input".format(test_dir), "plot_file={}/output".format(rt_plot_dir)],stdout=fstdout,stderr=fstderr)
        if (ret.returncode == 0): print("[PASS]")
        else:                     print("[FAIL]")
    
    
    
    if False:
        subprocess.run(["mkdir", "-p", bm_plot_dir])
        subprocess.run(["cp", "{}/output".format(rt_plot_dir), "{}/output".format(bm_plot_dir), "-rf" ])

    bm_plot_dir = None
    for bm_dir in sorted(glob("bm-*"),reverse=True):
        if (os.path.isdir(bm_dir + "/" + test['name'])): 
            print(bm_dir + "/" + test['name'])
            bm_plot_dir = bm_dir + "/" + test['name']
            break

    print("Comparing: [", rt_plot_dir, "] <==> [", bm_plot_dir, "]")

    #
    # Do a direct file-by-file comparison
    #
    match = True
    for rt in sorted(glob(rt_plot_dir+"/output/**",recursive=True)):
        if not os.path.isfile(rt): continue
        if os.path.basename(rt) in ["output", "metadata", "diff.html"]: continue
        bm = rt.replace(rt_plot_dir,bm_plot_dir)
        if not filecmp.cmp(bm,rt):
            print(bm, " does not match ", rt)
            match = False
    if (match) : print("OK - files match")
    else: print("Error - files do not match")
            

    #
    # Check timing
    #
    rt_metadata_f = open(rt_plot_dir+"/output/metadata","r")
    rt_run_time = float(re.findall("Simulation_run_time = (\S*)","".join(rt_metadata_f.readlines()))[0])

    bm_metadata_f = open(bm_plot_dir+"/output/metadata","r")
    bm_run_time = float(re.findall("Simulation_run_time = (\S*)","".join(bm_metadata_f.readlines()))[0])

    print("Run time: ",rt_run_time, " (this test) compared to ", bm_run_time," (benchmark)")
    if (rt_run_time > bm_run_time): print(" --> Performance increase!")
    if (rt_run_time < bm_run_time): print(" --> Performance decrease :-(")
    

    



    
    



    