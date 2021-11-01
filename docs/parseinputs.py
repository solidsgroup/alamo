#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join


def getParseDef(line):
    return "static void Parse" in line


def getParmParseDef(line):
    info = re.findall("ParmParse pp\(\"(.*)\"\)", line)
    if not len(info): return None
    #if len(info) != 1: return None
    return info[0]
    
def getParmParseInfo(line):
    # Try "query" or "queryarr"
    info = re.findall(r'pp.(queryarr\b|query\b)\(\"(.*)\",(.*)\);', line)
    if info: return info[0]

    # Try "queryclass"
    info = re.findall(r'pp.(queryclass\b)\(\"(.*)\",', line)
    if info:
        return [info[0][0],info[0][1],None]

    return None


def getDocs(lines,i):
    ret = None
    for j in range(0,len(lines)-i):
        if j>0 and getParseDef(lines[i-j]): break
        if j>0 and getParmParseDef(lines[i-j]): break
        if j>0 and getParmParseInfo(lines[i-j]): break
        docs = re.findall(r'.?\/\/(.*)', lines[i-j])
        if len(docs): 
            if not ret: ret = ""
            ret = docs[0] + ret
        elif j==0: continue
        else: break
    return ret


def extract(filename):
    sourcefile = open(filename)
    rets = dict()
    #rets[None] = dict()
    lines = sourcefile.readlines()
    prefix = None
    for i in range(len(lines)):
        if ("ParmParse pp" in lines[i]):
            prefix = getParmParseDef(lines[i])
            
            if prefix in rets.keys(): continue
                
            # Initialize the record for this ParmParser
            rets[prefix] = dict()
            rets[prefix]["items"] = []

            docs = getDocs(lines,i)
            if docs: rets[prefix]["docs"] = docs
            
        if (getParseDef(lines[i])):
            prefix = "[prefix]"
            rets[prefix] = dict()
            rets[prefix]["items"] = []

            docs = getDocs(lines,i)
            if docs: rets[prefix]["docs"] = docs
        if ('pp.query' in lines[i]):
            ret = dict()
            docs = re.findall(r'.?\/\/(.*)', lines[i-1])
            docs = getDocs(lines,i);
            if docs: ret["docs"] = docs
            info = getParmParseInfo(lines[i])
            if not info: continue
            ret["query"], ret["string"], ret["variable"] = info
            if prefix: ret["string"] = prefix + "." + ret["string"]
            rets[prefix]['items'].append(ret)
    return rets

#extract("../src/Integrator/Flame.cpp")
#exit()

docfile    = open("./source/Inputs.rst","w")
docfile.write(r"""
Inputs
------

""")

for dirname, subdirlist, filelist in os.walk("../src/"):
    for f in filelist:
        if f.endswith(".H") or f.endswith(".cpp"): 
            inputs = extract(dirname + "/" + f)
            if not len(inputs): continue

            classname = dirname.replace("../src/","").replace("/","::") + "::" + f.replace(".H","").replace(".cpp","")
            docfile.write(classname + "\n")
            docfile.write("".ljust(len(classname),"*")+"\n\n")
            docfile.write(".. flat-table:: \n")
            docfile.write("    :widths: 20 10 70\n")
            docfile.write("    :header-rows: 1\n\n")
            #docfile.write("    :widths: 1 1 1\n\n")
            docfile.write("    * - Parameter name\n")
            docfile.write("      - Type\n")
            docfile.write("      - Description\n")

            for prefix in inputs:
                if len(inputs[prefix]['items']) == 0: continue

                if ('docs' in inputs[prefix].keys()):
                    docfile.write("    * - {}  \n".format(inputs[prefix]['docs']))

                for item in inputs[prefix]['items']: 
                    docfile.write("    * - :code:`{}`\n".format(item['string']))
                    docfile.write("      - {}\n".format(item['query']))
                    docfile.write("      - {}\n".format(item['docs'] if 'docs' in item.keys() else ''))
                docfile.write("\n")
            docfile.write("\n")





