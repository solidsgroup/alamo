import os 
import re
from os import listdir
from os.path import isfile, join
from glob import glob


_PP_CALL = re.compile(
    r'\bpp(?P<separator>[._])(?P<method>[A-Za-z_]\w*)'
    r'(?:\s*<\s*(?P<template>[^<>]+?)\s*>)?\s*\('
)


def _matching_parenthesis(text, opening):
    """Return the closing parenthesis for a C++ call, or None if unbalanced."""
    depth = 0
    quote = None
    escaped = False

    for i in range(opening, len(text)):
        char = text[i]

        if quote:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
            continue

        if char in ['"', "'"]:
            quote = char
        elif char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            if depth == 0:
                return i

    return None


def _split_cpp_arguments(arguments):
    """Split C++ call arguments on commas outside (), [], {}, and strings."""
    ret = []
    start = 0
    depths = {'(': 0, '[': 0, '{': 0}
    closing = {')': '(', ']': '[', '}': '{'}
    quote = None
    escaped = False

    for i, char in enumerate(arguments):
        if quote:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
            continue

        if char in ['"', "'"]:
            quote = char
        elif char in depths:
            depths[char] += 1
        elif char in closing:
            opener = closing[char]
            depths[opener] = max(0, depths[opener] - 1)
        elif char == ',' and not any(depths.values()):
            ret.append(arguments[start:i].strip())
            start = i + 1

    ret.append(arguments[start:].strip())
    return ret


def _parse_pp_call(line, methods):
    """Parse a complete pp call while allowing nested C++ argument expressions."""
    for match in _PP_CALL.finditer(line):
        if match.group('method') not in methods:
            continue

        opening = match.end() - 1
        closing = _matching_parenthesis(line, opening)
        if closing is None:
            continue

        trailing = re.fullmatch(
            r'\s*;\s*(?://\s*(.*?))?\s*', line[closing + 1:], re.DOTALL
        )
        if not trailing:
            continue

        return {
            'method': match.group('method'),
            'template': match.group('template'),
            'arguments': _split_cpp_arguments(line[opening + 1:closing]),
            'doc': trailing.group(1) or '',
        }

    return None


def _string_argument(argument):
    match = re.fullmatch(r'\s*"([^"]+)"\s*', argument)
    return match.group(1) if match else None


def allclassnames(root = "../src/"):
    allclassnames = set()
    for dirname, subdirlist, filelist in os.walk(root):
        for f in filelist:
            f = dirname + '/' + f
            f = f.replace(".cpp","")
            f = f.replace(".H","")
            f = f.replace(".cc","")
            f = f.replace(root,"")
            allclassnames.add('::'.join(f.split('/')))
    return allclassnames


def getdocumentation(filename):
    sourcefile = None
    if os.path.isfile(filename+".H"):
        sourcefile = open(filename+".H")
    elif os.path.isfile(filename+".cc"):
        sourcefile = open(filename+".cc")
    else:
        return None
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
    for filename in [basefilename+".H",basefilename+".cpp",basefilename+".cc"]:
        if not os.path.isfile(filename): continue
        sourcefile = open(filename)
        lines = sourcefile.readlines()

        #group = None
        parsefn = False


        i = -1
        line = ""
        multiline = True

        for _i, _line in enumerate(lines):
            #
            # Logic to deal with multi-line expressions
            #
            if multiline: 
                if ";" in _line:
                    line += _line
                    multiline = False
                else:
                    line += _line.split('//')[0].replace('\n','')
                    continue
            elif (("pp." in _line.split('//')[0] or "pp_" in _line.split('//')[0]) and
                  "pp.contains" not in _line.split('//')[0]): # special case - don't care about contains stmts
                i = _i
                if ";" in _line:
                    line = _line
                else:
                    line = _line.split('//')[0].replace('\n','')
                    multiline = True
                    continue
            else:
                i = _i
                line = _line
            while "  " in line:
                line = line.replace("  "," ")

            

            pp = r'pp[._]'
            stringmatch = r'\s*"([^"]+)"\s*'
            variablematch = r',\s*[\w.]+\s*'
            docmatch = r'\s*;\s*(?:\/\/\s*(.*))?$'
            templatematch = r'\s*<([^>]+)>\s*'
            stringarraymatch = r',\s*\{(.*)\}\s*'
            nargs = r',*\s*([^)]*)'


            # Skip of pp. is commented out
            if "//" in line and line.find("//") < line.find("pp."): continue
            if "//" in line and line.find("//") < line.find("pp_"): continue


            # Catch standard pp.query and pp.queryarr inputs
            parsed = _parse_pp_call(line, {
                "query", "queryarr", "query_required", "queryarr_required", "query_file"
            })
            if parsed and len(parsed["arguments"]) in [2, 3]:
                name = _string_argument(parsed["arguments"][0])
            else:
                name = None
            if name is not None:
                query = dict()
                query["type"] = parsed["method"]
                query["string"] = name
                query["unit"] = parsed["arguments"][2] if len(parsed["arguments"]) == 3 else ""
                query["doc"] = parsed["doc"]
                query["file"] = filename
                query["line"] = i+1
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                rets.append(query)
                continue

            # Catch standard pp.query_default and pp.queryarr_default inputs
            parsed = _parse_pp_call(line, {"query_default", "queryarr_default"})
            if parsed and len(parsed["arguments"]) in [3, 4]:
                name = _string_argument(parsed["arguments"][0])
            else:
                name = None
            if name is not None:
                query = dict()
                query["type"] = parsed["method"]
                query["string"] = name
                query["default"] = parsed["arguments"][2]
                query["unit"] = parsed["arguments"][3] if len(parsed["arguments"]) == 4 else ""
                query["doc"] = parsed["doc"]
                query["file"] = filename
                query["line"] = i+1
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                rets.append(query)
                continue

            # Catch standard pp.query_validate
            match = re.findall(rf'^\s*{pp}query_validate\s*\({stringmatch}{variablematch}{stringarraymatch}\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                query = dict()
                query["type"] = "query_validate"
                query["string"] = match[0][0]
                query["possibles"] = match[0][1]
                query["doc"] = match[0][2]
                query["file"] = filename
                query["line"] = i+1
                query["default"] = True
                
                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                rets.append(query)
                continue

            # Catch pp.query_exactly, including nested unit expressions.
            parsed = _parse_pp_call(line, {"query_exactly"})
            if parsed and parsed["template"] and len(parsed["arguments"]) in [2, 3]:
                number = re.fullmatch(r'\s*(\d+)\s*', parsed["template"])
                possibles = re.fullmatch(r'\s*\{(.*)\}\s*', parsed["arguments"][0], re.DOTALL)
            else:
                number = None
                possibles = None
            if number and possibles:
                query = dict()
                query["type"] = "query_exactly"
                query["number"] = number.group(1)
                query["possibles"] = possibles.group(1)
                query["doc"] = parsed["doc"]
                query["file"] = filename
                query["line"] = i+1
                query["default"] = True

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    match = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if match: query["doc"] = match[0] + " " + query["doc"]
                    else: break
                rets.append(query)
                continue


            # Catch pp.queryclass inputs
            match = re.findall(rf'^\s*{pp}queryclass(?:<(.*)>)?\s*\(\s*"([^"]*)"(?:.*static_cast\s*<\s*(.*)\s*>.*)?[^)]*\s*\);\s*(?:\/\/\s*(.*)$)?',line)
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
                rets.append(queryclass)
                continue

            # Catch definition of a select function:

            match = re.findall(rf'{pp}(select[_default]*){templatematch}\({stringmatch},.*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
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
                rets.append(input)
                continue


            # Catch a queryclass_enumerate
            match = re.findall(rf'{pp}select_enumerate{templatematch}\({stringmatch}{variablematch}{nargs}\){docmatch}',line)
            if match:
                input = dict()
                input["type"] = "select_enumerate"
                input["class"] = match[0][0].replace(' ','') 
                input["string"] = match[0][1].replace(' ','')
                input["doc"] = match[0][2]
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch definition of a select_main function:
            match = re.findall(rf'{pp}select_main\s*<\s*([^.]+)>\s*\(.*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "select_main"
                input["classes"] = match[0][0].replace(' ','').split(',')
                input["string"] = None 
                input["doc"] = None 
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch definition of a select_only function:
            match = re.findall(rf'{pp}select_only\s*<\s*([^.]+)>\s*\(.*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "select_only"
                input["class"] = match[0][0]
                input["string"] = None 
                input["doc"] = None 
                input["file"] = filename
                input["line"] = i+1

                # # Check if previous lines have simple comments. Ignores "///" comments and
                # # any comment beginning with [
                # for j in reversed(range(0,i)):
                #     docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                #     if docmatch:
                #         input["doc"] = docmatch[0] + " " + input["doc"]
                #     else: break
                rets.append(input)
                continue

            # Catch definition of a queryclass function:
            match = re.findall(rf'{pp}queryclass\s*<\s*([^.]+)>\s*\(.*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "queryclass"
                input["class"] = match[0][0]
                input["string"] = None 
                input["doc"] = None 
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch a queryclass 
            match = re.findall(rf'{pp}queryclass<([^>]+)>\s*\("([^"]+)"\s*,\s*[a-z,A-Z,0-9,_,.]*\s*,*\s*[INFO]*\s*\)\s*;\s*(?:\/\/\s*(.*))?$',line)
            if match:
                input = dict()
                input["type"] = "queryclass"
                input["class"] = match[0][0].replace(' ','') 
                input["string"] = match[0][1].replace(' ','')
                input["doc"] = match[0][2]
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch a query_enumerate or queryarr_enumerate
            match = re.findall(rf'{pp}(query_enumerate|queryarr_enumerate)\s*\({stringmatch}{variablematch}{nargs}\){docmatch}',line)
            if match:
                input = dict()
                input["type"] = match[0][0]
                input["string"] = match[0][1].replace(' ','')
                input["doc"] = match[0][2]
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch a queryclass_enumerate
            match = re.findall(rf'{pp}queryclass_enumerate{templatematch}\({stringmatch}{variablematch}{nargs}\){docmatch}',line)
            if match:
                input = dict()
                input["type"] = "queryclass_enumerate"
                input["class"] = match[0][0].replace(' ','') 
                input["string"] = match[0][1].replace(' ','')
                input["doc"] = match[0][2]
                input["file"] = filename
                input["line"] = i+1

                # Check if previous lines have simple comments. Ignores "///" comments and
                # any comment beginning with [
                for j in reversed(range(0,i)):
                    docmatch = re.findall(r'^\s*\/\/(?!\/)(?!\s*\[)\s*(.*)',lines[j])
                    if docmatch:
                        input["doc"] = docmatch[0] + " " + input["doc"]
                    else: break
                rets.append(input)
                continue

            # Catch a queryclass with no ID
            match = re.findall(rf'{pp}queryclass{templatematch}\({stringmatch}{nargs}\){docmatch}$',line)
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
                rets.append(input)
                continue
                


            if "pp.contains" in line: continue
            if "pp.remove" in line:   continue
            if "pp.forbid" in line:   continue
            if "pp_forbid" in line:   continue
            if "pp.ignore" in line:   continue
            if "pp.add" in line:      continue
            if "pp.getEntries" in line:      continue
            if "pp.getPrefix" in line:      continue
            if "pp.prefix" in line:      continue
            if "pp.dumpTable" in line:      continue
            if "pp.countval" in line:      continue
            if "pp.AllUnusedInputs" in line:      continue
            if "pp.AnyUnusedInputs" in line:      continue
            if "pp.m_table" in line:      continue
            if "queryclass" in line and "*this" in line: continue
            if "query" and "c_str" in line: continue
            if "query" and "data()" in line: continue

            line = line.split('//')[0]
            line = line.split('#')[0]
            if "pp." in line or "pp_" in line:
                print("WARNING: ", line,end="")
                print("         ", basefilename, "\n")

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
            if f.endswith(".cc"): srcfileset.add(f.replace(".cc",""))
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
            path = []
            if dirname.replace(root,"") != "":
                path += dirname.replace(root,"").split('/')
            basefilename = f.replace(".H","").replace(".cpp","").replace(".cc","")
            path += [basefilename]
            classname = '::'.join(path)
            data[classname] = dict()

            try:
                data[classname]['inputs'] = extract(dirname+"/"+basefilename)
            except Exception as e:
                print("ERROR: problem reading",dirname)
                raise
            data[classname]['documentation'] = getdocumentation(dirname+"/"+basefilename)

            data[classname]['srcfile'] = None
            if os.path.isfile(f'{dirname}/{f}.cpp'):
                data[classname]['srcfile'] = f'src/{dirname.replace(root,"")}/{basefilename}.cpp'

            data[classname]['hdrfile'] = None
            if os.path.isfile(f'{dirname}/{f}.H'):
                data[classname]['hdrfile'] = f'src/{dirname.replace(root,"")}/{basefilename}.H'

            data[classname]['mainfile'] = None
            if os.path.isfile(f'{dirname}/{f}.cc'):
                data[classname]['mainfile'] = f'src/{dirname.replace(root,"")}/{basefilename}.cc'

    return data
