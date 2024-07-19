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

.. _tests:

==============================
:fas:`flask;fa-fw` Tests
==============================

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
        #docfile.write("    * - :icon-red:`error`\n\n")
        docfile.write("    * - :fas:`circle-xmark;sd-text-danger fa-fw fa-lg`\n\n")
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
        #docfile.write("    * - :icon-gray:`warning`\n\n")
        docfile.write("    * - :fas:`triangle-exclamation;sd-text-secondary fa-fw fa-lg`\n\n")
        if os.path.isfile(testdirname+"/Readme.rst"):
            docfile.write("      - :ref:`{}`\n".format(testname))
        else: 
            docfile.write("      - {}\n\n".format(testname))
    else:
        #docfile.write("    * - :icon-green:`check_circle`\n\n")
        docfile.write("    * - :fas:`circle-check;sd-text-success fa-fw fa-lg`\n\n")
        docfile.write("      - :ref:`{}`\n".format(testname))
        docfile.write("      - {}\n".format(str(len(config)-1)))
    
        has2D = False
        has3D = False
        for c in config:
            if "dim" in config[c]:
                if config[c]["dim"] == "2": has2D = True
                if config[c]["dim"] == "3": has3D = True
        
        dimstr = ""
        #if has2D: dimstr += ":icon:`2d` "
        if has2D: dimstr += ":fas:`maximize;fa-fw fa-lg sd-text-secondary` "
        #if has3D: dimstr += ":icon:`3d_rotation` "
        if has3D: dimstr += ":fab:`unity;fa-fw fa-lg sd-text-secondary` "
        docfile.write("      - {}\n".format(dimstr))
        
        if os.path.isfile(testdirname+"/test"):
            #docfile.write("      - :icon-green:`verified`\n")
            docfile.write("      - :fas:`medal;fa-fw fa-lg sd-text-secondary`\n")
        docfile.write("\n")

    if len(config) <= 1 and not os.path.isfile(testdirname+"/Readme.rst"):
        continue
    with open("Tests/{}.rst".format(testname),"w") as testdocfile:
        toctreestr += "   Tests/{}\n".format(testname)

        testdocfile.write(testname + "\n")
        testdocfile.write("="*len(testname) + "\n")

        if os.path.isfile(testdirname+"/Readme.rst"):
            testdocfile.write(".. include:: ../{}/Readme.rst\n\n\n".format(testdirname))

        for c in config:
            if c == "DEFAULT": continue
            #testsectionname = "[{}] {}".format(testname,c)
            testsectionname = c
            testdocfile.write(testsectionname+"\n")
            testdocfile.write("-"*len(testsectionname)+"\n")

            testdocfile.write(".. flat-table:: \n")
            testdocfile.write("    :widths: 10 90\n")
            testdocfile.write("    :header-rows: 0\n\n")

            #
            # DIMENSION
            #
            if config[c]["dim"] == "2":
                #testdocfile.write("    * - :icon:`2d`\n")
                testdocfile.write("    * - :fas:`maximize;fa-fw fa-lg`\n")
                testdocfile.write("      - Two-dimensional\n")
            else:
                config[c]["dim"] = "3"
                #testdocfile.write("    * - :icon:`3d_rotation`\n")
                testdocfile.write("    * - :fab:`unity;fa-fw fa-lg`\n")
                testdocfile.write("      - Three-dimensional\n")

            #
            # PARALLELISM
            #
            if "nprocs" in config[c] and int(config[c]["nprocs"]) > 1:
                #testdocfile.write("    * - :icon:`grid_view`\n")
                testdocfile.write("    * - :fas:`cubes;fa-fw fa-lg`\n")
                testdocfile.write("      - Parallel ({} procs)\n".format(config[c]["nprocs"]))
            else:
                #testdocfile.write("    * - :icon:`square`\n")
                testdocfile.write("    * - :fas:`cube;fa-fw fa-lg`\n")
                testdocfile.write("      - Serial\n")
            
            #
            # TESTING OR NOT TESTING
            #
            if not os.path.isfile("{}/test".format(testdirname)) or ("checK" in config[c] and config[c]["check"] in {"no","No","false","False","0"}):
                #testdocfile.write("    * - :icon:`report_off`\n")
                testdocfile.write("    * - :fas:`question;fa-fw fa-lg`\n")
                testdocfile.write("      - Not validated\n")
            else:
                #testdocfile.write("    * - :icon:`verified`\n")
                testdocfile.write("    * - :fas:`medal;fa-fw fa-lg`\n")
                testdocfile.write("      - Validated using check script\n")

            #
            # BENCHMARK TIME
            #
            if any(["benchmark-" in key for key in config[c]]):
                #testdocfile.write("    * - :icon:`timer`\n")
                testdocfile.write("    * - :fas:`stopwatch;fa-fw fa-lg`\n")
                testdocfile.write("      - ")
                for key in config[c]:
                    if "benchmark-" in key:
                        if "\n" in config[c][key]:
                            raise Exception("Error reading benchmark time for test {} section {}".format(testname,c))
                        testdocfile.write(config[c][key] + "s ({}) ".format(key.replace("benchmark-","")))
                testdocfile.write("\n")

            #testdocfile.write("    * - :icon:`play_circle`\n")
            testdocfile.write("    * - :fas:`circle-play;fa-fw fa-lg`\n")
            cmd = ""
            if "nprocs" in config[c] and int(config[c]["nprocs"]) > 1:
                cmd += "mpiexec -np {} ".format(config[c]["nprocs"])
            exe = "alamo"
            if "exe" in config[c]: exe = config[c]["exe"]
            cmd += "./bin/{}-{}d-g++".format(exe,config[c]["dim"])
            cmd += " {}/input".format(testdirname.replace("../../",""))
            if "args" in config[c]:
                cmd += " "
                cmdargs = [s.replace("= ","=").replace(" =","=") for s in config[c]["args"].split("\n")]
                for s in cmdargs:
                    if len(s.split('=')) == 2:
                        cmd += ' {}="{}"'.format(s.split('=')[0], s.split('=')[1])
            if "ignore" in config[c]:
                cmd += ' ignore="{}"'.format(config[c]["ignore"])

            testdocfile.write("      - .. code-block:: bash \n\n             {}\n".format(cmd))


            testdocfile.write("\n\n")
        
        #print(os.path.isfile("../../../{}/input".format(testdirname)))
        testdocfile.write(".. literalinclude:: ../{}/input\n".format(testdirname))
        testdocfile.write("   :caption: Input file ({}/input)\n".format(testdirname))
        testdocfile.write("   :language: makefile\n")

        
        
docfile.write("\n\n")
docfile.write(toctreestr)    
docfile.close()


