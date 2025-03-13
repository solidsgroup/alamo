import scraper
import json
import myast
import os
import glob

schema = dict()

allclassnames = None
templateclasses = None

def recurse(root,srcfile,f=None):
    
    intable = False
    def starttable():
        nonlocal intable
        if not intable:
            print("<table class='api-inputs-table'>", file=f)
            print(f"<thead><tr>", file=f)
            print(f"<th>Input name</th>", file=f)
            print(f"<th>Description</th>", file=f)
            print(f"</tr></thead>", file=f)

        intable=True
    def endtable():
        nonlocal intable
        if intable:
            print("</table><br/><br/>",file=f)
        intable=False

    global allclassnames, templateclasses

    if not allclassnames:
        allclassnames = scraper.allclassnames(root)
    
    def resolve(subname,classname):
        if type(subname) != list:
            subname = [subname]
        classname = classname.split("::")
        while len(classname):
            classname = classname[:-1]
            tempname = '::'.join(classname + subname)
            if tempname in allclassnames:
                return tempname
        return ""
    
    if not templateclasses:
        templateclasses = myast.scan(f"{root}/{srcfile}.cc")

    def extractTemplates(cl,fullclassname):
        classname = cl
        templates = {}
        if "<" in cl:
            classname = resolve(cl.split("<")[0],fullclassname)
            values = cl.split("<")[1].split(">")[0].replace(" ","").split(',')
            keys = templateclasses[classname]['templates']
            for key, value in zip(keys, values):
                templates[key] = value
        else:
            classname = resolve(classname,fullclassname)
        return classname, templates #, template_keys
    
    def templateReplace(string,template_dict):
        for key in template_dict:
            string = string.replace(key,template_dict[key])
        return string
    
    def getInputs(root,src,prefix=[],templates={},lev=0):
        classname = src.replace(root,"")
        classname = classname.replace('/','::')
        for input in scraper.extract(f"{root}/{src}"):
            if input['type'] == 'select_only':
                for cl in [input['class']]:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
    
                    inputname = subclassname
                    if (len(subclasstemplates)):
                        inputname += f"<{','.join(subclasstemplates[c] for c in subclasstemplates)}>"
    
                    endtable()

                    print("  "*lev,"<h3>" + inputname.replace('<','&lt;').replace('>','&gt;') + "</h3>", file=f)

                    starttable()

                    getInputs(root,subclassname.replace("::","/"),prefix, subclasstemplates,lev+1)
    
            elif input['type'] == 'select' or input['type'] == 'select_default':
    
                inputname = '.'.join(prefix + [input['string'],"type"])
                for cl in input['classes']:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
                    inputvalue = subclassname.split('::')[-1].lower()
                    ####print("  "*lev,f"[ if type = {name} ]")
    
                    print("  "*lev,f"<tr>", file=f)
                    print("  "*lev,f"  <td style='padding-left: {10*lev}px' colspan=2>", file=f)
                    print("  "*lev,f"     <b> if </b>", file=f)
                    print("  "*lev,f"     {inputname}", file=f)
                    print("  "*lev,f"     = {inputvalue}", file=f)
                    print("  "*lev,f"  </td>", file=f)
                    print("  "*lev,f"</tr>", file=f)
    
                    getInputs(root,
                              subclassname.replace("::","/"),
                              prefix + [input['string'], inputvalue],
                              subclasstemplates,lev+1)
    
            elif input['type'] == 'queryclass':
                subclassname, subclasstemplates = extractTemplates(templateReplace(input['class'],templates),classname)
                if input['string']:
                    getInputs(root,subclassname.replace("::","/"), prefix + [input['string']], subclasstemplates, lev)
                else:
                    getInputs(root,subclassname.replace("::","/"), prefix, subclasstemplates, lev)
    
            else:
                if 'string' in input:
                    if input['string']:
                        name = f'.'.join(prefix + [input['string']])
                        print("  "*lev,f"<tr>", file=f)
                        print("  "*lev,f"  <td style='padding-left: {10*lev}px'>", file=f)
                        print("  "*lev,f"    {name}", file=f)
                        print("  "*lev,f"  </td>", file=f)
                        print("  "*lev,f"  <td>", file=f)
                        print("  "*lev,f"    {input['doc']}", file=f)
                        print("  "*lev,f"  </td>", file=f)
                        print("  "*lev,f"</tr>", file=f)
    

    starttable()

    getInputs(root,srcfile)
    
    endtable()
