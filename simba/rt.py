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
import ansi2html



timestamp = datetime.today().strftime('%Y%m%d_%H%M%S')

parser = argparse.ArgumentParser(description='Sift through outputs')
parser.add_argument('inifile', help='Configuration file')
parser.add_argument('--benchmark',action='store_true',default=False,help='Set this run as benchmark for all tests')
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
    simba.updateTable(cur,s,types,verbose=False)
    print(config[s])
    for r in config[s]:
        print(r,config[s][r])

simba.updateRegTestTable(cur,"regtest",verbose=False)

for test in config.sections():

    status = simba.Status()

    run_id = timestamp

    #
    # Settings from Config File
    # 

    if 'input' in config[test]: input_file = config[test]['input']
    else:                       input_file = "tests/"+test+"/input"
    
    if 'dim' in config[test]:   dim = int(config[test]['dim'])
    else:                       dim = 3
    
    if 'nprocs' in config[test]: nprocs = int(config[test]['nprocs'])
    else:                        nprocs = 1



    #test_dir = "tests/" + test['name'] + "/"
    rt_plot_dir = "rt-" + run_id + "/" + test + "/"
    bm_plot_dir = "bm-" + run_id + "/" + test + "/"
    print("------------- rt_plot_dir: ", rt_plot_dir)
    print("------------- bm_plot_dir: ", bm_plot_dir)

    subprocess.run(["mkdir", "-p", rt_plot_dir])

    print("Running ...... ",end="")
    #if os.path.isfile(test_dir + "/input"): ## TODO enable this check
    fstdout = open(rt_plot_dir+"/"+"stdout","w")
    fstderr = open(rt_plot_dir+"/"+"stderr","w")
    ret = subprocess.run(["mpirun", "-np", str(nprocs),
                            "./bin/alamo-{}d-g++".format(dim), 
                            input_file, 
                            "plot_file={}/output".format(rt_plot_dir)],
                            stdout=fstdout,stderr=fstderr)
    fstdout.close()
    fstderr.close()
    conv = ansi2html.Ansi2HTMLConverter()
    for outfiletype in ["stdout","stderr"]:
        f = open(rt_plot_dir+"/"+outfiletype,"r")
        ansi = "".join(f.readlines())
        f.close()
        f = open(rt_plot_dir+"/output/"+outfiletype,"w")
        f.writelines(conv.convert(ansi))
        f.close()

    if (ret.returncode == 0): print("[PASS]")
    else:                     print("[FAIL]")
    status.runcode = ret.returncode
    
    
    if False:
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
        if os.path.basename(rt) in ["output", "metadata", "diff.html", "stdout", "stderr"]: continue
        bm = rt.replace(rt_plot_dir,bm_plot_dir)
        if not os.path.isfile(bm):
            print("Error - mismatched files")
            match = False
            break
        if not filecmp.cmp(bm,rt):
            print(bm, " does not match ", rt)
            match = False
            break
    if (match) : print("OK - files match")
    else: print("Error - files do not match")
    if (match) : status.compare = "YES"
    else : status.compare = "NO"
            
    #
    # Check timing
    #
    rt_metadata_f = open(rt_plot_dir+"/output/metadata","r")
    rt_run_time = float(re.findall("Simulation_run_time = (\S*)","".join(rt_metadata_f.readlines()))[0])

    bm_metadata_f = open(bm_plot_dir+"/output/metadata","r").readlines()
    bm_hash       = re.findall("HASH = (\S*)","".join(bm_metadata_f))[0]
    print("Benchmark hash is ", bm_hash)
    bm_run_time   = float(re.findall("Simulation_run_time = (\S*)","".join(bm_metadata_f))[0])

    status.performance = (rt_run_time - bm_run_time)/bm_run_time if bm_run_time>0 else 0

    
    data = simba.parse(rt_plot_dir+'/output')
    types = simba.getTypes(data)
    simba.updateTable(cur,test,types,verbose=False)
    simba.updateRecord(cur,test,data,verbose=False)

    simba.updateRegTestRecord(cur,"regtest",data['HASH'],run_id,test,status,bm_hash)
    db.commit()


db.close()
#exit()
    



    