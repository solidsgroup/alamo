import scraper
import json

data = scraper.scrape("./src/")

#print(data['hydro'])
#for i in data['hydro']['inputs']: print(i)
#exit(0)

#print(scraper.extract('./src/hydro'))#.extract("./src/hydro.cc")

#exit(0)
schema = dict()
def getInputs(cl,prefix=[]):
    for input in data[cl]["inputs"]:
        if input['type'] in ["query","queryarr","queryarr_required",
                             "query_default","query_required","query_validate","query_file"]:
            name = '.'.join(prefix + [input['string']])
            print(name)
            entry = {"description": input['doc']}
            schema[name] = entry
        elif input['type'] == "querysubclass":
            namespace = '::'.join(cl.split('::')[:-1])
            getInputs(namespace + '::' + input['class'].split('<')[0],prefix)
        elif input['type'] == "queryclass" and input['class']:
            if input['class'] in data.keys():
                theclass = input['class']
            else:
                namespace = '::'.join(cl.split('::')[:-1])
                theclass = namespace + '::' + input['class'].split('<')[0]
            getInputs(theclass,prefix + [input['string']])
        elif input['type'] in ["select","select_default"]:
            print('.'.join(prefix + [input['string'] , 'type']),
                  [subcl.split('::')[-1].lower() for subcl in input['classes']])
            for subcl in input['classes']:
                getInputs(subcl, prefix + [input['string'], subcl.split('::')[-1].lower()] )
        elif input['type'] in ["select_main"]:
            for subcl in input['classes']:
                getInputs(subcl, prefix)
        elif input['type'] == "queryclass":
            False # implicit substitution - not currently allowing
        



#print(data['hydro'])
print(getInputs("alamo"))
#f = open("input_schema.json","w")
#print(json.dump(schema,f))
#f.close()
