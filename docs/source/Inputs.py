#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join


def geticon(classname):
    if classname.startswith("BC"): return ":fas:`border-top-left;fa-fw` "
    if classname.startswith("IC"): return ":fas:`circle-right;fa-fw` "
    if classname.startswith("IO"): return ":fas:`print;fa-fw` "
    if classname.startswith("Integrator"): return ":fas:`gear;fa-fw` "
    if classname.startswith("Model"): return ":fas:`panorama;fa-fw` "
    if classname.startswith("Numeric"): return ":fas:`calculator;fa-fw` "
    if classname.startswith("Operator"): return ":far:`map;fa-fw` "
    if classname.startswith("Set"): return ":fas:`braille;fa-fw` "
    if classname.startswith("Solver"): return ":fas:`diamond-turn-right;fa-fw` "
    if classname.startswith("Util"): return ":fas:`sliders;fa-fw` "
    else: return ""

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

def extract(basefilename):
    rets = list()
    class inputdoc: pass
    for filename in [basefilename+".H",basefilename+".cpp"]:
        if not os.path.isfile(filename): continue
        sourcefile = open(filename)
        lines = sourcefile.readlines()
        group = None
        parsefn = False
        switch = None
        def reset():
            nonlocal rets, group, parsefn, switch
            if group:
                rets.append(group)
                group = None
            if switch:
                rets.append(switch)
                switch = None
        def insert(val):
            nonlocal rets, group, parsefn, switch
            val["parsefn"] = parsefn
            if group: group["inputs"].append(val)
            else: rets.append(val)

        for i, line in enumerate(lines):
            
            # Catch standard pp.query and pp.queryarr inputs
            match = re.findall('^\s*pp.(query[arr]*[_required]*)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
            if match:
                #print(match)
                query = dict()
                query["type"] = match[0][0]
                query["string"] = match[0][1]
                query["doc"] = match[0][2]
                query["file"] = filename
                query["line"] = i+1
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall('^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch standard pp.query_default and pp.queryarr_default inputs
            match = re.findall('^\s*pp.(query[arr]*_default*)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,\s*"*([^"]+)"*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
            if match:
                #print(match)
                query = dict()
                query["type"] = match[0][0]
                query["string"] = match[0][1]
                query["default"] = match[0][2]
                query["doc"] = match[0][3]
                query["file"] = filename
                query["line"] = i+1
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall('^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch standard pp.query_default and pp.queryarr_default inputs
            match = re.findall('^\s*pp.(query_validate)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,\s*\{(.*)\}\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
            if match:
                #print(match)
                query = dict()
                query["type"] = match[0][0]
                query["string"] = match[0][1]
                query["possibles"] = match[0][2]
                query["doc"] = match[0][3]
                query["file"] = filename
                query["line"] = i+1
                query["default"] = True
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall('^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch pp.queryclass inputs
            match = re.findall('^\s*pp.queryclass(?:<(.*)>)?\s*\(\s*"([^"]*)"(?:.*static_cast\s*<\s*(.*)\s*>.*)?[^)]*\);\s*(?:\/\/\s*(.*)$)?',line)
            if match:
                queryclass = dict()
                queryclass["type"] = "queryclass"
                queryclass["class"] = match[0][0]+match[0][2]
                queryclass["string"] = match[0][1]
                queryclass["doc"] = match[0][3]
                queryclass["file"] = filename
                queryclass["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall('^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: queryclass["doc"] = match[0] + " " + queryclass["doc"]
                    else: break

                insert(queryclass)
                continue

            # Catch ParmParse groups. This is old-fashioned but still supported.
            # All subsequent inputs will be nested under a group
            match = re.findall('^\s*(?:IO|amrex)::ParmParse\s*pp\s*(?:\(\s*"([^"]*)"\s*\))?\s*;\s*(?:\/\/\s*(.*))?',line)
            if match:
                reset()
                group = dict()
                group["type"] = "group"
                group["prefix"] = match[0][0]
                group["doc"] = match[0][1]
                if group["doc"] == "":
                    for j in reversed(range(0,i)):
                        match = re.findall('^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                        if match: group["doc"] = match[0] + " " + group["doc"]
                        else: break
                group["inputs"] = list()

            # Catch definition of a Parser function.
            if re.match('^\s*(?!\/\/).*Parse\(.*(?:IO|amrex)::ParmParse\s*&\s*(pp)\)(?!\s*;)',line):
                if parsefn: raise Exception(filename,i,"Multiple Parse functions cannot be declared in a single file")
                parsefn = True

            # Catch definition of a switch group
            match =re.findall('^\s*\/\/\s*\[\s*switch\s*(\S*)\s*]\s*(.*)',line)
            if match:
                reset()
                switch = dict()
                switch["type"] = "switch"
                switch["string"] = match[0][0]
                switch["doc"] = match[0][1]
                switch["inputs"] = list()
                continue

            if switch:
                # If we are inside a switch block, look for matches like
                #    else if (ic_type == "KEY") value.ic = new CLASS(args,"PREFIX")
                match = re.findall('^\s*(?:else)*\s*if\s*\(\S*\s*==\s*"([^"]*)"\s*\).*new\s*([A-Z,:,a-z,0-9,_]*)[^"]*"(.*)"',line)
                if match:
                    input=dict()
                    input["type"]="switchitem"
                    input["string"]=match[0][0]
                    input["class"]=match[0][1]
                    input["prefix"]=match[0][2]
                    switch["inputs"].append(input)

            # Close a switch group
            if re.match('^\s*\/\/\s*\[\s*end\s*(?:switch)?\s*]',line):
                reset()
            
        reset()

    return rets


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
                
            if documentation and writeFiles:
                docfile.write(documentation)
                
            if not len(inputs): continue
            if len(inputs) == 1 and not list(inputs)[0]: continue
    
    
            if writeFiles:
                docfile.write("\n\n")
                docfile.write(".. rst-class:: api-inputs-table\n\n")
                docfile.write(".. flat-table:: \n")
                docfile.write("    :widths: 40 10 40 10\n")
                docfile.write("    :header-rows: 1\n\n")
                docfile.write("    * - Parameter\n")
                docfile.write("      - Type\n")
                docfile.write("      - Description\n")
                docfile.write("      - Values\n")
    
            def writeInput(input,lev,prefix):
                prefix = list(filter(lambda x: x != "", prefix))
                if (input["type"]=="group"):
                    if writeFiles: docfile.write("    * - {}  \n".format(input['doc'].replace("\n", "\n        ")))
                    for subinput in input["inputs"]:
                        writeInput(subinput,lev+1,prefix + [input["prefix"]])
                if (input["type"] in ["query","queryarr","query_validate",
                                      "query_default","queryarr_default",
                                      "query_required","queryarr_required"]):
                    global num_tot, num_doc
                    if input["parsefn"]: prefix = ["[prefix]"] + prefix
                    num_tot += 1
                    if writeFiles:
                        docfile.write("    * - :code:`{}`\n".format('.'.join(prefix+[input['string']])))
                        docfilesearch.write("    * - :code:`{}`\n".format('.'.join(prefix+[input['string']])))
                        docfile.write("      - {}\n".format(input["type"]))
                        docfilesearch.write("      - :ref:`{}`\n".format(subhdr))
                    if input["doc"] != "":
                        num_doc += 1
                        if writeFiles:
                            docfile.write(      "      - {}\n".format(input['doc'].replace('\n','\n        ')))
                            docfilesearch.write("      - {}\n".format(input['doc'].replace('\n','\n        ')))
                            if "_default" in input["type"]:
                                docfile.write(      "      - :bdg-success:`{}`".format(input['default'].strip()))
                                docfilesearch.write("      - :bdg-success:`{}`".format(input['default'].strip()))
                            if "_required" in input["type"]:
                                docfile.write(      "      - :bdg-danger-line:`required`")
                                docfilesearch.write("      - :bdg-danger-line:`required`")
                            if "_validate" in input["type"]:
                                things = [d.replace('"',"").replace("'","").strip() for d in input['possibles'].split(',')]
                                string = ":bdg-success:`{}` ".format(things[0])
                                string += " ".join([":bdg-primary:`{}`".format(t) for t in things[1:]]) 
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


