#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join
from glob import glob

from scraper import getdocumentation, geticon, extract, scrapeInputs

src2url = {}
try:
	doxsourcefiles = sorted(glob("../build/html/doxygen/*source.html"))
	for doxsourcefile in doxsourcefiles:
	    with open(doxsourcefile) as f:
	        for line in f.readlines():
	            if r"<title>" in line:
	                line = line.replace(r"<title>Alamo: ","")
	                line = line.replace(r" Source File</title>","")
	                line = line.replace("\n","")
	                src2url[line] = doxsourcefile.replace("../build/html/","")
	                continue
except Exception as e:
    print(e)

inputsheader = r"""
.. _inputs: 
	
=============================
:fas:`cube;fa-fw` Inputs
=============================
	
"""


inputsearchheader = r"""
.. _inputssearch: 

===========================================
:fas:`magnifying-glass;fa-fw` Inputs Search
===========================================

Use the following box to search over all alamo inputs and descriptions.
This is the same list as in the :ref:`Inputs` section as generated by the autodoc system.

:bdg-primary-line:`Note` If searching for documentation on an alamo command that you found
in an input file, remember that the prefix may not be included in this table. For instance,
if you are looking for documentation on the following inputs

.. code-block:: bash

    hc.heat.alpha = 1.0
    el.nmodels = 2

the prefixes :code:`hc` and :code:`el` are not necessarily included, and you will not find them
if you do an exact search.
Instead, do a reverse search, starting with :code:`alpha` and :code:`nmodels`, then work out 
the prefixes out.


.. raw:: html

    <div class="input-group mb-3">
      <div class="input-group-prepend">
        <span class="input-group-text" id="inputGroup-sizing-default">
            <span class="fas fa-magnifying-glass fa-fw"></span>
        </span>
      </div>
      <input id='myInput' onkeyup='searchTable()' type="text" class="form-control" aria-label="Default" aria-describedby="inputGroup-sizing-default">
    </div>

    <!--input id='myInput' onkeyup='searchTable()' type='text'-->

    <script>
    function searchTable() {
        var input, filter, found, table, tr, td, i, j;
        input = document.getElementById("myInput");
        filter = input.value.toUpperCase();
        table = document.getElementsByClassName("api-inputs-table")[0];
        body = table.getElementsByTagName("tbody")[0];
        tr = body.getElementsByTagName("tr");
        for (i = 0; i < tr.length; i++) {
            td = tr[i].getElementsByTagName("td");
            for (j = 0; j < td.length; j++) {
                if (td[j].innerHTML.toUpperCase().indexOf(filter) > -1) {
                    found = true;
                }
            }
            if (found) {
                tr[i].style.display = "";
                found = false;
            } else {
                tr[i].style.display = "none";
            }
        }
    }
    </script>

    <style>
    table.api-inputs-table thead tr th:nth-child(1)
    {
        width:30%!important;
    }
    </style>



.. rst-class:: api-inputs-table

.. flat-table:: 
    :header-rows: 1

    * - Parameter
      - Namespace / Class
      - Description
"""

def scrapeInputs(root="../../src/", writeFiles=True):

    srcfiles = set()
    for dirname, subdirlist, filelist in sorted(os.walk(root)):
        for f in filelist:
            if f.endswith(".cpp"): srcfiles.add(dirname+"/"+f.replace(".cpp",""))
            if f.endswith(".H"): srcfiles.add(dirname+"/"+f.replace(".H",""))
    srcfiles = list(srcfiles)
    
    if writeFiles: docfilesearch    = open("InputsSearch.rst","w")
    if writeFiles: docfilesearch.write(inputsearchheader)
    
    
    if writeFiles: docfile    = open("Inputs.rst","w")
    if writeFiles: docfile.write(inputsheader)
    
    headerchar = ["=","*","-","~","."]
    written_headers = []
    
    global num_tot, num_doc
    num_tot = 0
    num_doc = 0

    for dirname, subdirlist, filelist in sorted(os.walk(root)):
        hdrname = dirname.replace(root,"").replace("/","::")
        depth = len(hdrname.split("::")) 
    
        srcfileset = set()
        for f in filelist:
            if f.endswith(".cpp"): srcfileset.add(f.replace(".cpp",""))
            if f.endswith(".H"): srcfileset.add(f.replace(".H",""))
        srcfilelist = list(srcfileset)
        
        #
        # This function makes sure pure abstract classes get
        # listed first.
        #
        def alphabetize_with_abstract_first(key):
            if key == hdrname.split("::")[-1]:
                return "0"
            return(key[0])
        for f in sorted(srcfilelist,key=alphabetize_with_abstract_first):

            try:
                inputs = extract(dirname+"/"+f)
            except Exception as e:
                print("ERROR: problem reading",dirname)
                raise
            documentation = getdocumentation(dirname+"/"+f)
            if not len(inputs) and not documentation:
                continue
    
            classname = dirname.replace(root,"").replace("/","::") + "::" + f.replace(".H","").replace(".cpp","")
    
            subhdr = ""
            for i in range(len(classname.split('::'))):
                subhdr = '::'.join(classname.split('::')[:i])
                if subhdr not in written_headers:
                    if writeFiles:
                        if '::' not in subhdr and subhdr != "":
                            docfile.write("--------------------\n\n\n")
                            docfile.write(".. _{}:\n\n".format(subhdr))
                            docfile.write(geticon(subhdr) + subhdr+"\n")
                            docfile.write("".ljust(len(geticon(subhdr)+subhdr),headerchar[i-1]))
                        else:
    	                    docfile.write("\n" + subhdr+"\n")
    	                    docfile.write("".ljust(len(subhdr),headerchar[i-1]))
                        docfile.write("\n\n\n")
                    written_headers.append(subhdr)
    
            if classname.split("::")[-1] != classname.split("::")[-2]:
                if writeFiles: docfile.write(classname + "\n")
                lev = len(classname.split('::'))-1
                if writeFiles: docfile.write("".ljust(len(classname),headerchar[lev])+"\n\n")
                subhdr = classname
                
            if writeFiles:
                srcfile = dirname.replace("../../","")+"/"+f+".cpp"
                hdrfile = dirname.replace("../../","")+"/"+f+".H"
                if srcfile in src2url.keys() or hdrfile in src2url.keys():
                    docfile.write("\n\n")
                    
                    if srcfile in src2url.keys():
                        docfile.write(r":bdg-link-primary-line:`{} <{}>`".format(srcfile,src2url[srcfile])+"\n")
                    if hdrfile in src2url.keys():
                        docfile.write(r":bdg-link-secondary-line:`{} <{}>`".format(hdrfile,src2url[hdrfile])+"\n")
                    docfile.write("\n\n")


            if documentation and writeFiles:
                docfile.write(documentation)
                
            if not len(inputs): continue
            if len(inputs) == 1 and not list(inputs)[0]: continue
    
    
            if writeFiles:
                docfile.write("\n\n")
                docfile.write(".. rst-class:: api-inputs-table\n\n")
                docfile.write(".. flat-table:: \n")
                docfile.write("    :header-rows: 1\n\n")
                docfile.write("    * - Parameter\n")
                docfile.write("      - Type\n")
                docfile.write("      - Values\n")
    
            def writeInput(input,lev,prefix):
                prefix = list(filter(lambda x: x != "", prefix))
                if (input["type"]=="group"):
                    if writeFiles: docfile.write("    * - {}  \n".format(input['doc'].replace("\n", "\n        ")))
                    for subinput in input["inputs"]:
                        writeInput(subinput,lev+1,prefix + [input["prefix"]])
                if (input["type"] in ["query","queryarr","query_validate",
                                      "query_default","queryarr_default",
                                      "query_required","queryarr_required",
                                      "query_file", "select", "select_default"]):
                    global num_tot, num_doc
                    if input["parsefn"]: prefix = ["[prefix]"] + prefix
                    num_tot += 1

                    # update the selects to derive the "type" that they should be
                    # and set the name to "type"
                    if input["type"] in ["select","select_default"]:
                        input['possibles'] = [str(cl).split('::')[-1].lower() for cl in input["classes"]]
                        input["string"] += ".type"
                        
                        


                    if writeFiles:
                        codetarget = None
                        try:
                            filename = src2url[input["file"].replace("../../","")]
                            linenumber = "l"+str(input["line"]).zfill(5)
                            codetarget = filename+"#"+linenumber
                        except Exception as e:
                            print(e)

                        if codetarget:
                            docfile.write(      "    * - :bdg-link-secondary:`{}<{}>`".format('.'.join(prefix+[input['string']]),codetarget) + "\n")
                            docfilesearch.write("    * - :bdg-link-secondary:`{}<{}>`".format('.'.join(prefix+[input['string']]),codetarget) + "\n")
                        else:
                            docfile.write(      "    * - :bdg-danger:`{}`".format('.'.join(prefix+[input['string']])) + "\n")
                            docfilesearch.write("    * - :bdg-danger:`{}`".format('.'.join(prefix+[input['string']])) + "\n")
                    if input["doc"] != "":
                        num_doc += 1
                        if writeFiles:
                            docfile.write(      "      - {}\n".format(input['doc'].replace('\n','\n        ')))
                            docfilesearch.write("      - {}\n".format(input['doc'].replace('\n','\n        ')))
                            if input["type"] in ["query_default","queryarr_default"]:
                                docfile.write(      "      - :bdg-success-line:`{}`".format(input['default'].strip()))
                                docfilesearch.write("      - :bdg-success-line:`{}`".format(input['default'].strip()))
                            if "_required" in input["type"]:
                                docfile.write(      "      - :bdg-danger-line:`required`")
                                docfilesearch.write("      - :bdg-danger-line:`required`")
                            if "_file" in input["type"]:
                                docfile.write(      "      - :bdg-secondary-line:`file path`")
                                docfilesearch.write("      - :bdg-secondary-line:`file path`")
                            if "_validate" in input["type"]:
                                things = [d.replace('"',"").replace("'","").strip() for d in input['possibles'].split(',')]
                                string = ":bdg-success-line:`{}` ".format(things[0])
                                string += " ".join([":bdg-primary-line:`{}`".format(t) for t in things[1:]]) 
                                docfile.write      (      "      - {}".format(string))
                                docfilesearch.write(      "      - {}".format(string))
                            if input["type"] in ["select","select_default"]:
                                things = input['possibles']
                                if input["type"] == "select_default":
                                    string = ":bdg-success-line:`{}` ".format(things[0])
                                    string += " ".join([":bdg-primary-line:`{}`".format(t) for t in things[1:]])
                                elif input["type"] == "select":
                                    string = " ".join([":bdg-primary-line:`{}`".format(t) for t in things])
                                docfile.write      (      "      - {}".format(string))
                                docfilesearch.write(      "      - {}".format(string))
                                
                    else:
                        print(input['file'],':',input['line'],' ',input['string'],' missing documentation')
                        if writeFiles:
                            docfile.write("      - \n")
                            docfilesearch.write("      - \n")
    
                if writeFiles:
                    docfile.write("\n")
                    docfilesearch.write("\n")
    
    
            for input in inputs:
                writeInput(input,0,[])
            
            if writeFiles:
                docfile.write("\n")
                docfilesearch.write("\n")
    
    
    print("\n{} of {} inputs documented\n".format(num_doc,num_tot))
    
    
    if writeFiles:
        docfile.close()
        docfilesearch.close()

    return num_doc, num_tot

