#!/usr/bin/env python3

import argparse
import os
from glob import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()

print("This script is still in progress!")

abspath = os.path.abspath(args.file)
print(abspath)

rootdir = os.path.dirname(abspath)
print(rootdir)

for leveldir in sorted(glob(rootdir+"/Level*")):
    print(leveldir.replace(rootdir+"/",""))
    f = open(leveldir + "/Cell_H")
    nodes = 0
    patches = 0
    for line in f.readlines():
        exp = re.compile('\(\(([0-9]*),([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*),([0-9]*)\)\)')
        res = exp.findall(line)
        if len(res) == 0 : continue
        #print(res[0])
        dimx = int(res[0][3])-int(res[0][0])
        dimy = int(res[0][4])-int(res[0][1])
        dimz = int(res[0][5])-int(res[0][2])
        #print(dimx,dimy,dimz,dimx*dimy*dimz)
        patches +=1 
        nodes += dimx*dimy*dimz
    print(nodes," nodes")
    print(patches, " patches")

    #exit()




