import os
import re
import hashlib

def getHash(self,directory):

    if not os.path.isfile(directory+'/metadata'):
        return None

    f = open(directory+"/metadata")
    for line in f.readlines():
        if line.startswith('#'): continue;
        if not 'HASH = ' in line: continue
        val = line.split(' = ')[1].replace('  ','').replace('\n','').replace(';','')
        return val

    f = open(directory+'/metadata')
    sim_hash = str(hashlib.sha224(f.read().encode('utf-8')).hexdigest())
    f.close()
    return sim_hash        
