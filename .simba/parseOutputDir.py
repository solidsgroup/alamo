import os
import re

def parseOutputDir(self,directory):

    if not os.path.isfile(directory+'/metadata'):
        return None

    things = dict()

    f = open(directory+"/metadata")
    for line in f.readlines():
        if line.startswith('#'): continue;
        if '::' in line:
            line = re.sub(r'\([^)]*\)', '',line)
            line = line.replace(" :: ", " = ").replace('[','').replace(',','').replace(']','').replace(' ','')
        if len(line.split(' = ')) != 2: continue;
        col = line.split(' = ')[0]#.replace('.','_')
        val = line.split(' = ')[1].replace('\n','')#.replace('  ','').replace('\n','').replace(';','')
        
        things[col] = val

    if os.path.isfile(directory+"/diff.html"):
        difffile = open(directory+"/diff.html")
        things['DIFF'] = difffile.read()
        difffile.close()

    if os.path.isfile(directory+"/diff.patch"):
        difffile = open(directory+"/diff.patch")
        things['DIFF_PATCH'] = difffile.read()
        difffile.close()

    if os.path.isfile(directory+"/stdout"):
        difffile = open(directory+"/stdout","r")
        things['STDOUT'] = difffile.read()
        difffile.close()

    if os.path.isfile(directory+"/stderr"):
        difffile = open(directory+"/stderr")
        things['STDERR'] = difffile.read()
        difffile.close()

    return things
