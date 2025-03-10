#!/usr/bin/env python3
import sys,os
me = os.path.dirname(__file__)
import scraper

def count(root):
    global num_tot, num_doc
    num_tot = 0
    num_doc = 0

    docs = scraper.scrape(root)
    for classname in docs:
        inputs = docs[classname]['inputs']
        def writeInput(input,lev,prefix):
            prefix = list(filter(lambda x: x != "", prefix))
            if (input["type"]=="group"):
                for subinput in input["inputs"]:
                    writeInput(subinput,lev+1,prefix + [input["prefix"]])
            elif (input["type"] in ["query","queryarr","query_validate",
                                  "query_default","queryarr_default",
                                  "query_required","queryarr_required",
                                  "query_file", "select", "select_default"]):
                global num_tot, num_doc
                num_tot += 1
                if input["doc"] != "": num_doc +=1
                elif 'file' in input.keys():
                    print(input['file'],':',input['line'],' ',input['string'],' missing documentation')
    
        for input in inputs:
            writeInput(input,0,[])
    print("\n{} of {} inputs documented\n".format(num_doc,num_tot))
    return num_doc, num_tot


docs, total = count(root=me+"/../src/")


if docs < total:
    raise Exception("Documentation is missing from inputs.")

