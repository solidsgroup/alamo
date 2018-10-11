#!/usr/bin/python
import os 
from os import listdir
from os.path import isfile, join
# for root, dirs, files in os.walk("../src/"):
#     path = root.split(os.sep)
#     for file in files:
#         if (file.endswith(".H") and file.replace(".H","") != path[-1]):
#             print("source/API/"+"/".join(path[2:])+"/"+file.replace(".H",".rst"))

base_dir = "test";

def fn(lev,path):
    contents = listdir(path)
    subdirs = []
    files = []
    for item in contents:
        if isfile(join(path,item)):
            if (not item.endswith(".H")): continue
            if (path.endswith(item.replace(".H",""))): continue
            print("   "*lev + item);
            files.append(item.replace(".H","")+".rst")
        else:
            print("   "*lev + item + "/" + item + ".rst");
            fn(lev+1,join(path,item))
            files.append(item+"/"+item+".rst")
        #print (subdirs, files)
    f = open(base_dir+"/"+path+".rst","w")

root = "../src";
fn(0,root);
