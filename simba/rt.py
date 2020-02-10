#!/usr/bin/env python3
import os
import re
from glob import glob
import subprocess
from datetime import datetime
import filecmp
import argparse
import configparser
import sqlite3
import simba

timestamp = datetime.today().strftime('%Y%m%d_%H%M%S')

parser = argparse.ArgumentParser(description='Sift through outputs')
parser.add_argument('inifile', help='Configuration file')
args = parser.parse_args()

db = sqlite3.connect('regtest.db')
db.text_factory = str
cur= db.cursor()
types = dict()

config = configparser.ConfigParser()
config.read(args.inifile)
print(config)
print(config.sections())
for s in config.sections():
    simba.updateTable(cur,s,types)
    print(config[s])
    for r in config[s]:
        print(r,config[s][r])

tests = [ { "name" : "TrigTest" } ,
          { "name" : "Eshelby"  } ]

#benchmark = False

for test in config.sections():

    print("[",test,"]")

    
    if 'input' in config[test]:
        input_file = config[test]['input']
    else:
        input_file = "tests/"+test+"/input"
    print("------------- Input file:  ", input_file)
    
    if 'dim' in config[test]:
        dim = int(config[test]['dim'])
    else:
        dim = 3
    print("------------- Dimnension:  ", dim)
    

    #test_dir = "tests/" + test['name'] + "/"
    rt_plot_dir = "rt-" + timestamp + "/" + test + "/"
    bm_plot_dir = "bm-" + timestamp + "/" + test + "/"
    print("------------- rt_plot_dir: ", rt_plot_dir)
    print("------------- bm_plot_dir: ", bm_plot_dir)

    subprocess.run(["mkdir", "-p", rt_plot_dir])

    print("Running ...... ",end="")
    #if os.path.isfile(test_dir + "/input"):
    fstdout = open(rt_plot_dir+"/"+"stdout","w")
    fstderr = open(rt_plot_dir+"/"+"stderr","w")
    ret = subprocess.run(["./bin/alamo-2d-g++", input_file, "plot_file={}/output".format(rt_plot_dir)],stdout=fstdout,stderr=fstderr)
    if (ret.returncode == 0): print("[PASS]")
    else:                     print("[FAIL]")
    
    
    
    if True:
        subprocess.run(["mkdir", "-p", bm_plot_dir])
        subprocess.run(["cp", "{}/output".format(rt_plot_dir), "{}/output".format(bm_plot_dir), "-rf" ])

    bm_plot_dir = None
    for bm_dir in sorted(glob("bm-*"),reverse=True):
        if (os.path.isdir(bm_dir + "/" + test)): 
            print(bm_dir + "/" + test)
            bm_plot_dir = bm_dir + "/" + test
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
    
    data = simba.parse(rt_plot_dir+'/output')
    types = simba.getTypes(data)
    simba.updateTable(cur,test,types)
    simba.updateRecord(cur,test,data)



db.commit()
db.close()
#exit()
    



    