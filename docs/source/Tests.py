#!/usr/bin/python
import os 
import glob
import re
from os import listdir
from os.path import isfile, join
import io
import configparser
from collections import OrderedDict


#
# Special order from SO - dictionary allows for keys to be specified multiple times and 
# config parser will read it all. (Copy/pasted from runtests.sh)
class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super().__setitem__(key, value)


#def icon(str)

docfile    = open("Tests.rst","w")
docfile.write(r"""
=====
Tests
=====

""")

headerchar = ["=","*","-","~","."]
written_headers = []

num_tot = 0
num_doc = 0


docfile.write(r"""


""")

#docfile.write(r"- A regular icon: :material-outlined:`data_exploration;2em`, some more text")
#docfile.write(r""":raw:'<span class="material-symbols-outlined">search</span>'""")
#docfile.write(r""":html:'<span class="material-symbols-outlined">add</span>'""")

docfile.write("\n\n")

docfile.write(".. flat-table:: \n")
docfile.write("    :widths: 3 15 10 10 10\n")
docfile.write("    :header-rows: 1\n\n")
docfile.write("    * - Status\n")
docfile.write("      - Name\n")
docfile.write("      - Sections\n")
docfile.write("      - Dimension\n")
docfile.write("      - Validation\n")

if not os.path.isdir("Tests"):
    os.mkdir("Tests")

toctreestr  = ".. toctree::\n"
toctreestr += "   :hidden:\n\n"

for testdirname in sorted(glob.glob("../../tests/*")):
    if not os.path.isdir(testdirname): continue
    
    testname = os.path.basename(testdirname)

    if not os.path.isfile(testdirname+"/input"):
        docfile.write("    * - :icon-red:`error`\n\n")
        docfile.write("      - {}\n\n".format(testname))
        continue




    # Parse the input file ./tests/MyTest/input containing #@ comments.
    # Everything commeneted with #@ will be interpreted as a "config" file
    # The variable "config" is a dict of dicts where each item corresponds to
    # a test configuration.
    cfgfile = io.StringIO()
    input = open(testdirname + "/input")
    for line in input.readlines():
        if line.startswith("#@"):
            cfgfile.write(line.replace("#@",""))
    cfgfile.seek(0)
    config = configparser.ConfigParser(dict_type=MultiOrderedDict,strict=False)
    config.read_file(cfgfile)

    if len(config) <= 1:
        docfile.write("    * - :icon-gray:`warning`\n\n")
        docfile.write("      - {}\n\n".format(testname))
        continue

    docfile.write("    * - :icon-green:`check_circle`\n\n")
    docfile.write("      - :ref:`{}`\n".format(testname))
    docfile.write("      - {}\n".format(str(len(config)-1)))

    has2D = False
    has3D = False
    for c in config:
        if "dim" in config[c]:
            if config[c]["dim"] == "2": has2D = True
            if config[c]["dim"] == "3": has3D = True
    
    dimstr = ""
    if has2D: dimstr += ":icon:`2d` "
    if has3D: dimstr += ":icon:`3d_rotation` "
    docfile.write("      - {}\n".format(dimstr))
    
    if os.path.isfile(testdirname+"/test"):
        docfile.write("      - :icon-green:`verified`\n")
    docfile.write("\n")

    with open("Tests/{}.rst".format(testname),"w") as testdocfile:
        toctreestr += "   Tests/{}\n".format(testname)

        testdocfile.write(testname + "\n")
        testdocfile.write("="*len(testname) + "\n")

        for c in config:
            if c == "DEFAULT": continue
            testdocfile.write(c+"\n")
            testdocfile.write("-"*len(c)+"\n")

            testdocfile.write(".. flat-table:: \n")
            testdocfile.write("    :widths: 10 90\n")
            testdocfile.write("    :header-rows: 0\n\n")

            #
            # DIMENSION
            #
            if config[c]["dim"] == "2":
                testdocfile.write("    * - :icon:`2d`\n")
                testdocfile.write("      - Two-dimensional\n")
            else:
                config[c]["dim"] = "3"
                testdocfile.write("    * - :icon:`3d_rotation`\n")
                testdocfile.write("      - Three-dimensional\n")

            #
            # PARALLELISM
            #
            if "nprocs" in config[c] and int(config[c]["nprocs"]) > 1:
                testdocfile.write("    * - :icon:`grid_view`\n")
                testdocfile.write("      - Parallel ({} procs)\n".format(config[c]["nprocs"]))
            else:
                testdocfile.write("    * - :icon:`square`\n")
                testdocfile.write("      - Serial\n")
            
            #
            # TESTING OR NOT TESTING
            #
            if not os.path.isfile("{}/test".format(testdirname)) or ("checK" in config[c] and config[c]["check"] in {"no","No","false","False","0"}):
                testdocfile.write("    * - :icon:`report_off`\n")
                testdocfile.write("      - No testing\n")
            else:
                testdocfile.write("    * - :icon:`verified`\n")
                testdocfile.write("      - Testing\n")

            #
            # BENCHMARK TIME
            #
            if any(["benchmark-" in key for key in config[c]]):
                testdocfile.write("    * - :icon:`timer`\n")
                testdocfile.write("      - ")
                for key in config[c]:
                    if "benchmark-" in key:
                        testdocfile.write(config[c][key] + "s ({}) ".format(key.replace("benchmark-","")))
                testdocfile.write("\n")

            testdocfile.write("    * - :icon:`play_circle`\n")
            cmd = "./bin/alamo-{}d-g++".format(config[c]["dim"])
            cmd += " {}/input".format(testdirname.replace("../../",""))
            if "args" in config[c]:
                cmd += " "
                cmd += " ".join([arg for arg in config[c]["args"].split()])
            if "ignore" in config[c]:
                cmd += " ignore={}".format(config[c]["ignore"])

            testdocfile.write("      - :code:`{}`\n".format(cmd))



            testdocfile.write("\n\n")
        
        #print(os.path.isfile("../../../{}/input".format(testdirname)))
        testdocfile.write(".. literalinclude:: ../{}/input\n".format(testdirname))
        testdocfile.write("   :caption: Input file ({}/input)\n".format(testdirname))
        testdocfile.write("   :language: makefile\n")

        
        
docfile.write("\n\n")
docfile.write(toctreestr)

    

exit(0)

for dirname, subdirlist, filelist in os.walk("../../tests/"):
    hdrname = dirname.replace("../../src/","").replace("/","::")
    depth = len(hdrname.split("::")) 

    continue

    write_header = True
    
    srcfilelist = set() 
    for f in filelist:
        if f.endswith(".cpp"): srcfilelist.add(f.replace(".cpp",""))
        if f.endswith(".H"): srcfilelist.add(f.replace(".H",""))
    
    for f in sorted(srcfilelist):
        
        if True: 
            try:
                inputs = extract(dirname+"/"+f)
            except Exception as e:
                print("ERROR: problem reading",dirname)
                raise
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


