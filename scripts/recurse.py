import scraper
import json
import my_ast
import os

#data = scraper.scrape("../src/")
#templateclasses = my_ast.ast_classes("../src/alamo.cc")
#print(templateclasses.keys())

schema = dict()
#def getInputs(cl,prefix=[]):
#
#    classname = cl
#    template_values = []
#    template_keys = []
#
#    if "<" in cl:
#        classname = cl.split("<")[0]
#        template_values = cl.split("<")[1].split(">")[0].replace(" ","").split(',')
#        template_keys = templateclasses[classname]['templates']
#        print("TEMPLATE:", classname, template_keys)
#        print(f"         - values = {template_values}")
#        print(f"         - keys   = {template_keys}")
#
#
#    for input in data[classname]["inputs"]:
#        if input['type'] in ["query","queryarr","queryarr_required",
#                             "query_default","query_required","query_validate","query_file"]:
#            name = '.'.join(prefix + [input['string']])
#            entry = {"description": input['doc']}
#            schema[name] = entry
#        elif input['type'] == "querysubclass":
#            namespace = '::'.join(classname.split('::')[:-1])
#            getInputs(namespace + '::' + input['class'].split('<')[0],prefix)
#        elif input['type'] == "queryclass" and input['class']:
#            if input['class'] in data.keys():
#                theclass = input['class']
#            else:
#                namespace = '::'.join(classname.split('::')[:-1])
#                theclass = namespace + '::' + input['class'].split('<')[0]
#            getInputs(theclass,prefix + [input['string']])
#        elif input['type'] in ["select","select_default"]:
#            for subcl in input['classes']:
#                getInputs(subcl, prefix + [input['string'], subcl.split('::')[-1].lower()] )
#        elif input['type'] in ["select_main"]:
#            for subcl in input['classes']:
#                getInputs(subcl, prefix)
#        elif input['type'] == "queryclass":
#            False # implicit substitution - not currently allowing


allnamespaces = set()
for dirname, subdirlist, filelist in os.walk("../src/"):
    for f in filelist:
        f = dirname + '/' + f
        f = f.replace(".cpp","")
        f = f.replace(".H","")
        f = f.replace(".cc","")
        f = f.replace("../src/","")
        allnamespaces.add('::'.join(f.split('/')))
def resolve(subname,classname):
    classname = classname.split("::")
    while len(classname):
        classname = classname[:-1]
        tempname = '::'.join(classname + [subname])
        if tempname in allnamespaces:
            return tempname
    return ""

templateclasses = my_ast.ast_classes("../src/alamo.cc")
def extractTemplates(cl,fullclassname):
    classname = cl
    templates = {}
    #template_values = []
    #template_keys = []
    
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
        if input['type'] == 'select_main':
            for cl in input['classes']:
                subclassname, subclasstemplates = extractTemplates(cl,classname)

                if (len(subclasstemplates)):
                    print("  "*lev,f"{subclassname}<{','.join(subclasstemplates[c] for c in subclasstemplates)}>")
                else:
                    print("  "*lev,f"{subclassname}",)
                getInputs(root,subclassname.replace("::","/"),prefix, subclasstemplates,lev+1)
        elif input['type'] == 'queryclass':
            subclassname, subclasstemplates = extractTemplates(templateReplace(input['class'],templates),classname)
            print("  "*lev, input['string'],input['type'],input['string'])
            if input['string']:
                getInputs(root,subclassname.replace("::","/"), prefix + [input['string']], subclasstemplates, lev+1)
            else:
                getInputs(root,subclassname.replace("::","/"), prefix, subclasstemplates, lev)

        else:
            if 'string' in input:
                print("  "*lev,'.'.join(prefix + [input['string']]),input['type'])
        
getInputs("../src","alamo")


#f = open("./src/alamo.cc")
#scraper.sanitize(f.readlines())

#print(data['hydro']) 
#getInputs("alamo")
#scraper.extract('./src/alamo')
#print(stuff)
#print(schema)
#for s in schema:
#    print(s)
#f = open("input_schema.json","w")
#print(json.dump(schema,f))
#f.close()
