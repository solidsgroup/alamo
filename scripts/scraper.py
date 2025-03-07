import os 
import re
from os import listdir
from os.path import isfile, join
from glob import glob

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
            match = re.findall(r'^\s*pp.(query[arr]*[_required]*[_file]*)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,*\s*[INFO]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
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
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch standard pp.query_default and pp.queryarr_default inputs
            match = re.findall(r'^\s*pp.(query[arr]*_default*)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,\s*"*([^"^,]+)"*\s*,*\s*[INFO]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
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
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch standard pp.query_default and pp.queryarr_default inputs
            match = re.findall(r'^\s*pp.(query_validate)\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,\s*\{(.*)\}\s*,*\s*[INFO]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',lines[i])
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
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                insert(query)
                continue

            # Catch pp.queryclass inputs
            match = re.findall(r'^\s*pp.queryclass(?:<(.*)>)?\s*\(\s*"([^"]*)"(?:.*static_cast\s*<\s*(.*)\s*>.*)?[^)]*,*\s*[INFO]*\s*\);\s*(?:\/\/\s*(.*)$)?',line)
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
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: queryclass["doc"] = match[0] + " " + queryclass["doc"]
                    else: break

                insert(queryclass)
                continue

            # Catch ParmParse groups. This is old-fashioned but still supported.
            # All subsequent inputs will be nested under a group
            match = re.findall(r'^\s*(?:IO|amrex)::ParmParse\s*pp\s*(?:\(\s*"([^"]*)"\s*\))?\s*;\s*(?:\/\/\s*(.*))?',line)
            if match:
                reset()
                group = dict()
                group["type"] = "group"
                group["prefix"] = match[0][0]
                group["doc"] = match[0][1]
                if group["doc"] == "":
                    for j in reversed(range(0,i)):
                        match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                        if match: group["doc"] = match[0] + " " + group["doc"]
                        else: break
                group["inputs"] = list()

            # Catch definition of a Parser function.
            if re.match(r'^\s*(?!\/\/).*Parse\(.*(?:IO|amrex)::ParmParse\s*&\s*(pp)\)(?!\s*;)',line):
                if parsefn: raise Exception(filename,i,"Multiple Parse functions cannot be declared in a single file")
                parsefn = True


            # Catch definition of a select function:
            match = re.findall(r'pp\.(select[_default]*)<([^>]+)>\s*\("([^"]+)"\s*,[^)]*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = match[0][0]
                input["classes"] = match[0][1].replace(' ','').split(',')
                input["string"] = match[0][2].replace(' ','')
                input["doc"] = match[0][3]
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                insert(input)

            # Catch a queryclass 
            match = re.findall(r'pp\.queryclass<([^>]+)>\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,*\s*[INFO]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "queryclass"
                input["class"] = match[0][0].replace(' ','') 
                input["string"] = match[0][1].replace(' ','')
                input["doc"] = match[0][2]

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                insert(input)
                continue

            # Catch a queryclass with no ID
            match = re.findall(r'pp\.queryclass<([^>]+)>\s*\([^\)]*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            match = re.findall(r'pp\.queryclass<([^\)_]+)>\s*\([^\)]*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "querysubclass"
                input["class"] = match[0][0].replace(' ','') 
                input["doc"] = match[0][1]

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                insert(input)
            


            # Catch definition of a switch group
            match =re.findall(r'^\s*\/\/\s*\[\s*switch\s*(\S*)\s*]\s*(.*)',line)
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
                match = re.findall(r'^\s*(?:else)*\s*if\s*\(\S*\s*==\s*"([^"]*)"\s*\).*new\s*([A-Z,:,a-z,0-9,_]*)[^"]*"(.*)"',line)
                if match:
                    input=dict()
                    input["type"]="switchitem"
                    input["string"]=match[0][0]
                    input["class"]=match[0][1]
                    input["prefix"]=match[0][2]
                    switch["inputs"].append(input)

            # Close a switch group
            if re.match(r'^\s*\/\/\s*\[\s*end\s*(?:switch)?\s*]',line):
                reset()
            
        reset()

    return rets

def scrape(root="../src/"):
    headerchar = ["=","*","-","~","."]
    written_headers = []
    
    global num_tot, num_doc
    num_tot = 0
    num_doc = 0

    data=dict()

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
            classname = dirname.replace(root,"").replace("/","::") + "::" + f.replace(".H","").replace(".cpp","")
            data[classname] = dict()

            try:
                data[classname]['inputs'] = extract(dirname+"/"+f)
            except Exception as e:
                print("ERROR: problem reading",dirname)
                raise
            data[classname]['documentation'] = getdocumentation(dirname+"/"+f)

            data[classname]['srcfile'] = None
            if os.path.isfile(f'{dirname}/{f}.cpp'):
                data[classname]['srcfile'] = f'src/{dirname.replace(root,"")}/{f}.cpp'

            data[classname]['hdrfile'] = None
            if os.path.isfile(f'{dirname}/{f}.H'):
                data[classname]['hdrfile'] = f'src/{dirname.replace(root,"")}/{f}.H'

    return data

for s in scrape():
    print(s)
