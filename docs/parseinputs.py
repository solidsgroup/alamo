#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join




def extract(sourcefile):
    rets = []
    lines = sourcefile.readlines()
    prefix = None
    for i in range(len(lines)):
        if ("ParmParse pp" in lines[i]):
            info = re.findall("ParmParse pp\(\"(.*)\"\)", lines[i])
            if len(info): prefix = info[0]
            else: prefix = None
        if ("static void Parse" in lines[i]):
            prefix = "[prefix]"
        if ('pp.query' in lines[i]):
            ret = dict()
            docs = re.findall(r'.?\/\/(.*)', lines[i-1])
            info = re.findall(r'pp.(queryarr\b|query\b)\(\"(.*)\",(.*)\);', lines[i])
            #if not len(info): print(lines[i],i)
            if len(info):
                if len(info[0]) == 3:
                    ret["query"], ret["string"], ret["variable"] = info[0]
                    if prefix: ret["string"] = prefix + "." + ret["string"]
                    ret["docs"] = docs[0] if len(docs) else "None"
                    rets.append(ret)
    return rets



docfile    = open("./source/Inputs.rst","w")
docfile.write(r"""
Inputs
------

.. list-table::
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
""")

for dirname, subdirlist, filelist in os.walk("../src/"):
    #print(filelist)
    for f in filelist:
        if f.endswith(".H"): 
            inputs = extract(open(dirname + "/" + f))
            if not inputs: continue
            classname = dirname.replace("../src/","").replace("/","::") + "::" + f.replace(".H","")
            docfile.write(classname + "\n")
            docfile.write("".ljust(len(classname),"*")+"\n\n")

            docfile.write(".. list-table:: \n")
            docfile.write("    :header-rows: 1\n\n")
            docfile.write("    * - Parameter name\n")
            docfile.write("      - Type\n")
            docfile.write("      - Description\n")
            
            for input in inputs: 
                docfile.write("    * - :code:`{}`\n".format(input['string']))
                docfile.write("      - {}\n".format("Single value" if input['query'] == "query" else "Array"))
                docfile.write("      - {}\n".format(input['docs']))

                #print(input)



#sourcefile = open("../src/IC/Laminate.H","r")

#lines = sourcefile.readlines()
#for i in range(len(lines)):
#    if ('pp.query' in lines[i]):
#        docs = re.findall(r'.?\/\/(.*)', lines[i-1])
#        print(docs)
#        info = re.findall(r'pp.(queryarr|query)\(\"(.*)\",(.*)\);', lines[i])[0]
#        if len(info) == 3:
#            query, string, variable = info
#            print(query, string, variable)
#            docfile.write("    * - :code:`{}`\n".format(string))
#            docfile.write("      - {}\n".format("Single value" if query == "query" else "Array"))
#            docfile.write("      - {}\n".format(docs[0] if docs else "None"))
#
##            docfile.write(query + string + variable)        
#        
#    if ("pp.query") in lines[i]:
#        print(lines[i])



