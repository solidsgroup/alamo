import os
import re
import hashlib

def parse(directory):

    if not os.path.isfile(directory+'/metadata'):
        return None

    things = dict()

    things['DIR'] = os.path.abspath(directory)

    f = open(directory+"/metadata")
    for line in f.readlines():
        if line.startswith('#'): continue;
        if '::' in line:
            line = re.sub(r'\([^)]*\)', '',line)
            line = line.replace(" :: ", " = ").replace('[','').replace(',','').replace(']','').replace(' ','')
        if len(line.split(' = ')) != 2: continue;
        col = line.split(' = ')[0].replace('.','_')
        val = line.split(' = ')[1].replace('  ','').replace('\n','').replace(';','')
        
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

    if not 'HASH' in things:
        f = open(directory+'/metadata')
        sim_hash = str(hashlib.sha224(f.read().encode('utf-8')).hexdigest())
        f.close()
        things['HASH'] = sim_hash        

    return things
