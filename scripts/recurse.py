import scraper
import json
import myast
import os
import glob

schema = dict()

allclassnames = None
templateclasses = None

class DefaultPrinter:
    def __init__(self):
        True
    def __del__(self):
        True
    def starttable(self):
        True
    def endtable(self):
        True
    def printinput(self,input,prefix,lev,classes=[]):
        True
    def printconditionalstart(self,input,prefix,lev,classes=[]):
        True
    def printconditional(self,inputname,inputvalue,lev,classes=[]):
        True
    def printconditionalend(self,inputname,inputvalue,lev,classes=[]):
        True
    def printtablename(self,inputname,lev):
        True


def recurse(root,srcfile,printer=DefaultPrinter(),f=None):

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
    
    templateclasses = myast.scan(root,f"{srcfile}.cc")

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
    
            elif input['type'] in ['select', 'select_default']:
    
                inputname = '.'.join(prefix + [input['string'],"type"])

                printer.printconditionalstart(input,prefix,lev)

                for cl in input['classes']:
                    subclassname, subclasstemplates = extractTemplates(cl,classname)
                    inputvalue = subclassname.split('::')[-1].lower()
    
                    printer.printconditional(inputname,inputvalue,lev)
    
                    getInputs(root,
                              subclassname.replace("::","/"),
                              prefix + [input['string'], inputvalue],
                              subclasstemplates,lev+1)

                    printer.printconditionalend(input,prefix,lev)
    
            elif input['type'] in ['queryclass']:
                subclassname, subclasstemplates = extractTemplates(templateReplace(input['class'],templates),classname)
                if input['string']:
                    getInputs(root,subclassname.replace("::","/"), prefix + [input['string']], subclasstemplates, lev)
                else:
                    getInputs(root,subclassname.replace("::","/"), prefix, subclasstemplates, lev)
                    
            elif input['type'] in ['queryclass_enumerate']:
                subclassname, subclasstemplates = extractTemplates(templateReplace(input['class'],templates),classname)
                getInputs(root,subclassname.replace("::","/"), prefix + [input['string']+r'#'], subclasstemplates, lev)

            else:
                printer.printinput(input,prefix,lev)    

    printer.starttable()

    getInputs(root,srcfile)
    
    printer.endtable()
