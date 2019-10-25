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

f = open(args.file)
f.readline() #Hyperclaw stuff
nfields = int(f.readline());  print("# fields:     " + str(nfields))
for n in range(0,nfields): f.readline()
dimension = int(f.readline()) # Spatial dimension maybe?
print("Dimension: ", dimension)
f.readline() # Time maybe?
namrlevs = int(f.readline()); print("# AMR Levels: " + str(namrlevs))
domainlo = [float(d) for d in f.readline().replace(' \n','').split(' ')]; print("Domain Lo Vect: ", domainlo)
domainhi = [float(d) for d in f.readline().replace(' \n','').split(' ')]; print("Domain Hi Vect: ", domainhi)
rootvolume = (domainhi[0]-domainlo[0])
if dimension > 1: rootvolume *= (domainhi[1]-domainlo[1])
if dimension > 2: rootvolume *= (domainhi[2]-domainlo[2])
f.readline() # List of ref ratios for all AMR levels
f.readline() # List of int domains for all AMR levels
f.readline() # List of refinements for all AMR levels
dx = []
for n in range(0,namrlevs+1): 
    dx.append([float(d) for d in f.readline().replace(' \n','').split(' ')])

rootdir = os.path.dirname(abspath)
print(rootdir)

#for leveldir in sorted(glob(rootdir+"/Level*")):
print('{:>10}'.format("Level"),end="")
print('{:>15}'.format("# Nodes"),end="")
print('{:>15}'.format("# Patches"),end="")
print('{:>20}'.format("Volume"),end="")
print('{:>20}'.format("% Volume"),end="")
print()

for amrlev in range(0,namrlevs+1):
    leveldir = rootdir+"/Level_"+str(amrlev)+"/"
    #print(leveldir)
    f = open(leveldir + "/Cell_H")
    nodes = 0
    patches = 0 
    volume = 0
    for line in f.readlines():
        
        if dimension == 2: 
            exp = re.compile('\(\(([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*)\)\)')
        if dimension == 3:
            exp = re.compile('\(\(([0-9]*),([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*),([0-9]*)\) \(([0-9]*),([0-9]*),([0-9]*)\)\)')
        res = exp.findall(line)
        if len(res) == 0 : continue
        if dimension == 2:
            dimx = int(res[0][2])-int(res[0][0])
            dimy = int(res[0][3])-int(res[0][1])
            volume += dimx*dx[amrlev][0] * dimy*dx[amrlev][1]
            nodes += dimx*dimy
        if dimension == 3:
            dimx = int(res[0][3])-int(res[0][0])
            dimy = int(res[0][4])-int(res[0][1])
            dimz = int(res[0][5])-int(res[0][2])
            volume += dimx*dx[amrlev][0] * dimy*dx[amrlev][1] * dimz*dx[amrlev][2]
            nodes += dimx*dimy*dimz
        patches +=1 
        
    print('{:>10}'.format(amrlev),end="")
    print('{:>15}'.format(nodes),end="")
    print('{:>15}'.format(patches),end="")
    print('{:>20.5}'.format(volume),end="")
    print('{:>20.5}'.format(100*volume/rootvolume),end="")
    print()



