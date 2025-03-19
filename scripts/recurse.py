import scraper
import json
import myast
import os
import glob

schema = dict()

allclassnames = None
templateclasses = None

def recurse(root,srcfile,printer,f=None):

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
    
    def getInputs(root,src,prefix=[],templates={},lev=0, classes=[]):
        classname = src.replace(root,"")
        classname = classname.replace('/','::')
        for input in scraper.extract(f"{root}/{src}"):
            if input['type'] == 'select_only':
                for cl in [input['class']]:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
    
                    inputname = subclassname
                    if (len(subclasstemplates)):
                        inputname += f"<{','.join(subclasstemplates[c] for c in subclasstemplates)}>"
    
                    printer.endtable()

                    printer.printtablename(inputname,lev)

                    printer.starttable()

                    getInputs(root,subclassname.replace("::","/"),prefix, subclasstemplates,lev+1)
    
            elif input['type'] == 'select' or input['type'] == 'select_default':
    
                inputname = '.'.join(prefix + [input['string'],"type"])

                printer.printconditionalstart(input,prefix,lev)

                for cl in input['classes']:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
                    inputvalue = subclassname.split('::')[-1].lower()
    
                    printer.printconditional(inputname,inputvalue,lev)

                    #print(">>>",inputname,inputvalue)
    
                    getInputs(root,
                              subclassname.replace("::","/"),
                              prefix + [input['string'], inputvalue],
                              subclasstemplates,lev+1)

                    printer.printconditionalend(input,prefix,lev)

    
            elif input['type'] == 'queryclass':
                subclassname, subclasstemplates = extractTemplates(templateReplace(input['class'],templates),classname)
                if input['string']:
                    getInputs(root,subclassname.replace("::","/"), prefix + [input['string']], subclasstemplates, lev)
                else:
                    getInputs(root,subclassname.replace("::","/"), prefix, subclasstemplates, lev)
    
            else:
                if 'string' in input:
                    printer.printinput(input,prefix,lev)    

    printer.starttable()

    getInputs(root,srcfile)
    
    printer.endtable()
