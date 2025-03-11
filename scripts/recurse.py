import scraper
import json
import ast
import os

schema = dict()


def recurse(root,srcfile, printer):

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
    
    templateclasses = ast.scan(f"{root}/{srcfile}.cc")
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
        classname = src.replace("../src/","")
        classname = classname.replace('/','::')
        for input in scraper.extract(f"{root}/{src}"):
            if input['type'] == 'select_only':
                for cl in [input['class']]:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
    
                    inputname = subclassname
                    if (len(subclasstemplates)):
                        inputname += f"<{','.join(subclasstemplates[c] for c in subclasstemplates)}>"
    
                    print("  "*lev,f"<thead>")
                    print("  "*lev,f"<th colspan=2>")
                    print("  "*lev,inputname.replace('<','&lt;').replace('>','&gt;'))
                    print("  "*lev,f"</th>")
                    getInputs(root,subclassname.replace("::","/"),prefix, subclasstemplates,lev+1)
    
            elif input['type'] == 'select' or input['type'] == 'select_default':
    
                inputname = '.'.join(prefix + [input['string'],"type"])
                for cl in input['classes']:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
                    inputvalue = subclassname.split('::')[-1].lower()
                    ####print("  "*lev,f"[ if type = {name} ]")
    
                    print("  "*lev,f"<tr>")
                    print("  "*lev,f"  <td style='padding-left: {10*lev}px'>")
                    print("  "*lev,f"     {inputname}")
                    print("  "*lev,f"     = {inputvalue}")
                    print("  "*lev,f"  </td>")
                    print("  "*lev,f"</tr>")
    
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
                    printer.printinput(input)
                    if input['string']:
                        name = f'.'.join(prefix + [input['string']])
                        print("  "*lev,f"<tr>")
                        print("  "*lev,f"  <td style='padding-left: {10*lev}px'>")
                        print("  "*lev,f"    {name}")
                        print("  "*lev,f"  </td>")
                        print("  "*lev,f"  <td>")
                        print("  "*lev,f"    {input['doc']}")
                        print("  "*lev,f"  </td>")
                        print("  "*lev,f"</tr>")
    
    print("<html>")
    print("<table>")
    
    getInputs(root,srcfile)
    
    print("</table>")
    print("</html>")



class htmlprinter:
    def printinput(input):
        print(input)
    

recurse("../src/","alamo",htmlprinter)

