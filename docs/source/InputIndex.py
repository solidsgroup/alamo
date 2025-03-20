import sys
import os
sys.path.append(os.path.abspath('../../scripts/'))
import glob
import recurse
import re
import scraper

src2url = {}
try:
	doxsourcefiles = sorted(glob.glob("../build/html/doxygen/*source.html"))
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
    print("Error processing doxygen file...",e)

def codetarget(file,line):
    target = ""
    try:
        filename = src2url[file.replace('//','/').replace("../../","")]
        linenumber = "l"+str(line).zfill(5)
        target = filename+"#"+linenumber
    except Exception as e:
        print("Could not find URL for this file: ",file)
    return target



class HTMLPrinter:
    tbody_cntr = 0

    f = None
    intable = False
    mypr = "     "
    def __init__(self,exe):
        self.f = open(f"InputIndex_{exe}.rst","w",encoding='utf-8')
        print(exe,file=self.f)
        print("--------------------------",file=self.f)
        print("",file=self.f)
        print("",file=self.f)
        print(scraper.getdocumentation(f"../../src/{exe}"), file=self.f)
        print("",file=self.f)
        print("",file=self.f)
        print(".. raw:: html",file=self.f)
        print("",file=self.f)

        print(self.mypr,"""<script>""",file=self.f)
        print(self.mypr,"""function showTab(allName,tabName, button) {""",file=self.f)
        print(self.mypr,"""    document.querySelectorAll('.' + allName).forEach(tbody => tbody.classList.remove('active')); """,file=self.f)
        print(self.mypr,"""    document.querySelector('.' + tabName).classList.add('active'); """,file=self.f)
        print(self.mypr,"""    console.log(tabName) """,file=self.f)
        print(self.mypr,"""    document.querySelectorAll('.btn-' + allName).forEach(btn => btn.classList.remove('active')); """,file=self.f)
        print(self.mypr,"""    button.classList.add('active');""",file=self.f)
        print(self.mypr,"""}""",file=self.f)
        print(self.mypr,"""</script>""",file=self.f)


    def __del__(self):
        if self.intable:
            print(self.mypr,"</table>", file=self.f)
        print(r"",file=self.f) 
        # apparently sphinx won't include mathjax unless there is an RST math directive.
        # So we have to put this in to make sure mathjax renders.
        print(r":math:`\quad`",file=self.f) 
        print(r"",file=self.f) 
        self.f.close()
    def starttable(self):
        if not self.intable:
            print(self.mypr,"<table class='api-inputs-table'>", file=self.f)
            print(self.mypr,f"<thead><tr>", file=self.f)
            print(self.mypr,f"<th>Input name</th>", file=self.f)
            print(self.mypr,f"<th>Description</th>", file=self.f)
            print(self.mypr,f"<th>Value</th>", file=self.f)
            print(self.mypr,f"</tr></thead>", file=self.f)
        self.intable=True
    def endtable(self):
        if self.intable:
            print(self.mypr,"</table><br/><br/>",file=self.f)
        self.intable=False
    def printinput(self,input,prefix,lev,classes=[]):
        if input['string']:
            name = f'.'.join(prefix + [input['string']])
            input_classes = "sd-sphinx-override sd-badge sd-bg-secondary sd-bg-text-secondary reference external"
            bdg_success   = "sd-sphinx-override sd-badge sd-outline-success sd-text-success"
            bdg_primary   = "sd-sphinx-override sd-badge sd-outline-primary sd-text-primary"
            bdg_secondary = "sd-sphinx-override sd-badge sd-outline-secondary sd-text-secondary"
            bdg_danger    = "sd-sphinx-override sd-badge sd-outline-danger sd-text-danger"

            srcfile = input['file']
            line = input['line']

            # Convert RST math directeves to plain mathjax
            processed_doc = re.sub(r":math:`(.*?)`", r"\\(\1\\)", input['doc'].replace('\n',''))
            # Convert RST code directives to <code>
            processed_doc = re.sub(r":code:`(.*?)`", r"<code>\1</code>", processed_doc)


            print(self.mypr,"  "*lev,f"<tr>", file=self.f)
            print(self.mypr,"  "*lev,f"  <td class='input-col' style='padding-left: {10*(lev+1)}px'>", file=self.f)
            print(self.mypr,"  "*lev,f'    <a href="{codetarget(srcfile,line)}" class="{input_classes}"><span>{name}</span></a>',file=self.f)
            print(self.mypr,"  "*lev,f"  </td>", file=self.f)
            print(self.mypr,"  "*lev,f"  <td>", file=self.f)
            print(self.mypr,"  "*lev,f"    {processed_doc}", file=self.f)
            print(self.mypr,"  "*lev,f"  </td>", file=self.f)


            if input['type'] in ["query","queryarr"]:
                print(self.mypr,"  "*lev,f"  <td>", file=self.f)
                msg = "Input does not have default or requirement indicator and may be undefind if not specified."
                print(self.mypr,"  "*lev,f"    <p><span title='{msg}' class='fas fa-exclamation-triangle fa-fw'></span> </p>", file=self.f)
                print(self.mypr,"  "*lev,f"  </td>", file=self.f)

            if input['type'] in ["query_required","queryarr_required"]:
                print(self.mypr,"  "*lev,f"  <td>", file=self.f)
                print(self.mypr,"  "*lev,f"    <p><span class='{bdg_danger}'>required</span></p>", file=self.f)
                print(self.mypr,"  "*lev,f"  </td>", file=self.f)
                
            if input['type'] in ["query_file"]:
                print(self.mypr,"  "*lev,f"  <td>", file=self.f)
                print(self.mypr,"  "*lev,f"    <p><span class='{bdg_secondary}'>file path</span></p>", file=self.f)
                print(self.mypr,"  "*lev,f"  </td>", file=self.f)


            if input['type'] in ["query_default","queryarr_default"]:
                print(self.mypr,"  "*lev,f"  <td>", file=self.f)
                print(self.mypr,"  "*lev,f"    <p><span class='{bdg_success}'>{input['default']}</span></p>", file=self.f)
                print(self.mypr,"  "*lev,f"  </td>", file=self.f)

            if input['type'] in ["query_validate"]:
                things = [d.replace('"',"").replace("'","").strip() for d in input['possibles'].split(',')]
                print(self.mypr,"  "*lev,f"  <td><p>", file=self.f)
                print(self.mypr,"  "*lev,f"    <span class='{bdg_success}'>{things[0]}</span>", file=self.f)
                for thing in things[1:]:
                    print(self.mypr,"  "*lev,f"    <span class='{bdg_primary}'>{thing}</span>", file=self.f)
                print(self.mypr,"  "*lev,f"  </p></td>", file=self.f)
            

            print(self.mypr,"  "*lev,f"</tr>", file=self.f)

    def printconditionalstart(self,input,prefix,lev,classes=[]):
        self.tbody_cntr += 1

        name = f'.'.join(prefix + [input['string']])
        input_classes = "sd-sphinx-override sd-badge sd-bg-secondary sd-bg-text-secondary reference external"
        bdg_success   = "sd-sphinx-override sd-badge sd-outline-success sd-text-success"
        bdg_primary   = "sd-sphinx-override sd-badge sd-outline-primary sd-text-primary"
        bdg_secondary = "sd-sphinx-override sd-badge sd-outline-secondary sd-text-secondary"
        bdg_danger    = "sd-sphinx-override sd-badge sd-outline-danger sd-text-danger"
        
        input['possibles'] = [str(cl).split('::')[-1].lower() for cl in input["classes"]]
        srcfile = input['file']
        line = input['line']

        print(self.mypr,"  "*lev,f"<tr class='conditional-start-first'>", file=self.f)
        print(self.mypr,"  "*lev,f"  <td rowspan=2 style='padding-left: {10*(lev+1)}px'>", file=self.f)
        print(self.mypr,"  "*lev,f'    <a href="{codetarget(srcfile,line)}" class="{input_classes}"><span>{name}.type</span></a>',file=self.f)
        print(self.mypr,"  "*lev,f"  </td>", file=self.f)
        print(self.mypr,"  "*lev,f"  <td colspan=2>", file=self.f)
        print(self.mypr,"  "*lev,f"    <p>{input['doc']}</p>", file=self.f)
        #print(self.mypr,"  "*lev,f"  </td>", file=self.f)
        #print(self.mypr,"  "*lev,f"<tr>", file=self.f)

        things = input['possibles']
        #print(self.mypr,"  "*lev,f"<tr class='conditional-start-second'>", file=self.f)

        #print(self.mypr,"  "*lev,f"  <td colspan=2><p>", file=self.f)

        for thing in things:
            conditional_class_all = "n" + str(self.tbody_cntr) + "-" + name.replace('.','-') + "-type"
            conditional_class_thing = "n" + str(self.tbody_cntr) + "-" + name.replace('.','-') + "-type-" + thing
            jscript_cmd = f"""showTab("{conditional_class_all}\",\"{conditional_class_thing}",this)"""
            print(self.mypr,"  "*lev,f"     <button class='btn-{conditional_class_all} {bdg_primary}' onclick='{jscript_cmd}'>{thing}</button>", file=self.f)
        print(self.mypr,"  "*lev,f"  </p></td>", file=self.f)
        print(self.mypr,"  "*lev,f"</tr>", file=self.f)



    def printconditional(self,inputname,inputvalue,lev,classes=[]):
        input_classes = "sd-sphinx-override sd-badge sd-bg-secondary sd-bg-text-secondary reference external"
        bdg_success   = "sd-sphinx-override sd-badge sd-outline-success sd-text-success"
        bdg_primary   = "sd-sphinx-override sd-badge sd-outline-primary sd-text-primary"
        bdg_secondary = "sd-sphinx-override sd-badge sd-outline-secondary sd-text-secondary"
        bdg_danger    = "sd-sphinx-override sd-badge sd-outline-danger sd-text-danger"

        conditional_class_all = "n" + str(self.tbody_cntr) + "-" + inputname.replace('.','-')
        conditional_class_thing = "n" + str(self.tbody_cntr) + "-" + inputname.replace('.','-') + "-" + inputvalue

        print(self.mypr,"  "*lev,f"<tbody class='conditional_class_inputs {conditional_class_all} {conditional_class_thing}'>",file=self.f)

    def printconditionalend(self,inputname,inputvalue,lev,classes=[]):
        input_classes = "sd-sphinx-override sd-badge sd-bg-secondary sd-bg-text-secondary reference external"
        bdg_success   = "sd-sphinx-override sd-badge sd-outline-success sd-text-success"
        bdg_primary   = "sd-sphinx-override sd-badge sd-outline-primary sd-text-primary"
        bdg_secondary = "sd-sphinx-override sd-badge sd-outline-secondary sd-text-secondary"
        bdg_danger    = "sd-sphinx-override sd-badge sd-outline-danger sd-text-danger"

        print(self.mypr,"  "*lev,f"</tbody>",file=self.f)

    def printtablename(self,inputname,lev):
        sanitizedname = inputname.replace('<','&lt;').replace('>','&gt;')
        print(self.mypr,"  "*lev,f"<h3> {sanitizedname}  </h3>", file=self.f)

for exe in ["alamo","mechanics","hydro","sfi","thermoelastic","topop"]:
    try:
        printer = HTMLPrinter(exe)
        recurse.recurse("../../src/",exe,printer=printer)
        del printer
    except Exception as e:
        del printer
        print("Problem writing file for ",exe)
        print(e)
        raise(e)


