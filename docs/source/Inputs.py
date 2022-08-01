#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join


def getdocumentation(filename):
    sourcefile = open(filename+".H")
    ret = ""
    for line in sourcefile.readlines():
        if line.startswith(r"///"): # special provision for legacy doxygen comments
            ret += line.split(r"///")[1]
        elif line.startswith(r"// "):
            ret += line.split(r"// ")[1]
        elif line.startswith(r"//"):
            ret += line.split(r"//")[1]
        else:
            return ret
    return ret


def getParseDef(line):
    return "static void Parse" in line


def getParmParseDef(line):
    info = re.findall("ParmParse pp\(\"(.*)\"\)", line.split(r'//')[0])
    if not len(info): return None
    #if len(info) != 1: return None
    return info[0]
    
def getParmParseInfo(line):
    # Try "query" or "queryarr"
    info = re.findall(r'pp.(queryarr\b|query\b)\(\"(.*)\",(.*)\);', line.split(r'//')[0])
    if info: return info[0]

    # Try "queryclass"
    info = re.findall(r'pp.(queryclass\b)\(\"(.*)\",', line.split(r'//')[0])
    if info:
        return [info[0][0],info[0][1],None]

    return None


def getDocs(lines,i):
    ret = None
    for j in range(0,len(lines)-i):
        if j>0 and r'//' not in lines[i-j]: break
        if j>0 and getParseDef(lines[i-j]): break
        if j>0 and getParmParseDef(lines[i-j]): break
        if j>0 and getParmParseInfo(lines[i-j]): break
        #if j>0 and not lines[i-j].split(r'/')[0].isspace(): break
        docs = re.findall(r'.?\/\/(.*)', lines[i-j])
        if len(docs):
            if not ret: ret = ""
            if docs[0].endswith('/'): ret = docs[0] + ret
            else:                     ret = docs[0] + "\n" + ret
        elif j==0: continue
        else: break
    return ret


def extract(basefilename):
    rets = dict()
    for filename in [basefilename+".H",basefilename+".cpp"]:
        if not os.path.isfile(filename): continue
        sourcefile = open(filename)
        lines = sourcefile.readlines()
        prefix = None
        for i in range(len(lines)):
            if ("ParmParse pp" in lines[i].split(r'//')[0]):
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
                docs = getDocs(lines,i);
                if docs: ret["docs"] = docs
                info = getParmParseInfo(lines[i])
                if not info: continue
                ret["query"], ret["string"], ret["variable"] = info
                if prefix: ret["string"] = prefix + "." + ret["string"]
                rets[prefix]['items'].append(ret)
    return rets

docfile    = open("Inputs.rst","w")
docfile.write(r"""
======
Inputs
======

""")

headerchar = ["=","*","-","~","."]
written_headers = []

num_tot = 0
num_doc = 0

for dirname, subdirlist, filelist in os.walk("../../src/"):
    hdrname = dirname.replace("../../src/","").replace("/","::")
    depth = len(hdrname.split("::")) 

    write_header = True
    
    srcfilelist = set()
    for f in filelist:
        if f.endswith(".cpp"): srcfilelist.add(f.replace(".cpp",""))
        if f.endswith(".H"): srcfilelist.add(f.replace(".H",""))
    
    for f in sorted(srcfilelist):
        
        if True: 
            inputs = extract(dirname+"/"+f)
            documentation = getdocumentation(dirname+"/"+f)
            if not len(inputs) and not documentation:
                continue

            classname = dirname.replace("../../src/","").replace("/","::") + "::" + f.replace(".H","").replace(".cpp","")

            for i in range(len(classname.split('::'))):
                subhdr = '::'.join(classname.split('::')[:i])
                if subhdr not in written_headers:
                    docfile.write(subhdr+"\n")
                    docfile.write("".ljust(len(subhdr),headerchar[i-1]))
                    docfile.write("\n\n\n")
                    written_headers.append(subhdr)


            docfile.write(classname + "\n")
            lev = len(classname.split('::'))-1
            docfile.write("".ljust(len(classname),headerchar[lev])+"\n\n")
            
            if documentation:
                docfile.write(documentation)
            
            if not len(inputs): continue

            docfile.write(".. flat-table:: \n")
            docfile.write("    :widths: 20 10 70\n")
            docfile.write("    :header-rows: 1\n\n")
            docfile.write("    * - Parameter name\n")
            docfile.write("      - Type\n")
            docfile.write("      - Description\n")

            for prefix in inputs:
                #if len(inputs[prefix]['items']) == 0: continue

                if ('docs' in inputs[prefix].keys()):
                    docfile.write("    * - {}  \n".format(inputs[prefix]['docs'].replace("\n", "\n        ")))

                for item in inputs[prefix]['items']:
                    num_tot += 1
                    docfile.write("    * - :code:`{}`\n".format(item['string']))
                    docfile.write("      - {}\n".format(item['query']))
                    if 'docs' in item.keys():
                        num_doc += 1
                        docfile.write("      - {}\n".format(item['docs'].replace('\n','\n        ')))
                    else:
                        docfile.write("      - \n")
                docfile.write("\n")
            docfile.write("\n")


print("\n{} of {} inputs documented\n".format(num_doc,num_tot))


docfile.close()


