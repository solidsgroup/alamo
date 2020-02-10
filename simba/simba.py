import os
import re
import sqlite3


class Status:
    runcode = -1
    compare = "NO"
    performance = 0.0

#
# Scan an output directory and bundle all of the stuff in a dictionary
#
def parse(directory):

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

    if os.path.isfile(directory+"/stdout"):
        difffile = open(directory+"/stdout","r")
        things['STDOUT'] = difffile.read()
        difffile.close()

    if os.path.isfile(directory+"/stderr"):
        difffile = open(directory+"/stderr")
        things['STDERR'] = difffile.read()
        difffile.close()

    return things

#
# Scan through an entry's data dictionary and determine datatypes
#
def getTypes(data):
    types = dict()
    for key in data:
        val = data[key]
        if (key == 'HASH'): 
            continue
        if (key == 'DIFF' or key == 'STDOUT' or key == 'STDERR'):
            types[key] = 'BLOB'
        if key not in types:
            try:
                int(val)
                types[key] = 'INTEGER'
            except ValueError:
                False
        if key not in types:
            try:
                float(val)
                types[key] = 'FLOAT'
            except ValueError:
                False
        if key not in types:
            types[key] = 'VARCHAR(1000)'     
    return types 

#
# Create a table if it does not exist, or update the table if it does.
#
def updateTable(cur,tablename,types,mode="results",verbose=True):
    # If the table does not exist, create it
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cur.fetchall()]
    if "__tables__" not in tables:
        cur.execute('CREATE TABLE __tables__ (NAME UNIQUE, Description VARCHAR(8000));')
    if tablename not in tables:
        if mode == "results":
            cur.execute('CREATE TABLE ' + tablename + ' ('
                        'HASH VARCHAR(255) UNIQUE, ' +
                        'DIR,' +
                        'Description VARCHAR(8000),' +
                        'Tags VARCHAR(1000)' +
                        (',' if len(types)>0 else '') +
                        ','.join([key+' '+types[key] for key in sorted(types)]) +
                        ');')
        elif mode == "regtest":
            cur.execute('CREATE TABLE ' + tablename + ' ('
                        'HASH VARCHAR(255) UNIQUE ' +
                        (',' if len(types)>0 else '') +
                        ','.join([key+' '+types[key] for key in sorted(types)]) +
                        ');')
        if (verbose): print('\033[1;32mADDED TABLE\033[1;0m: ' + tablename)
    else:
        if (verbose): print('\033[1;34mUSING TABLE\033[1;0m: ' + tablename)

    cur.execute('SELECT * FROM __tables__ WHERE NAME = \"' + tablename + '\";')
    if (len(cur.fetchall()) == 0): cur.execute('INSERT INTO __tables__ (NAME, Description) VALUES (\"' + tablename + '\", \"Description\");')

    # If the table exists, but new columns have been added, modify the table
    # accordingly
    cur.execute("PRAGMA table_info("+tablename+")")
    columns=[a[1] for a in cur.fetchall()]
    for key in types:
        if key not in columns:
            cur.execute('ALTER TABLE ' + tablename + ' ADD ' + key + ' ' + types[key])
            if (verbose): print('\033[1;34mADDED COLUMN\033[1;0m: ' + key + ' to ' + tablename)


def updateRecord(cur,tablename,data,hash=None,verbose=True):
    if not hash: hash = data['HASH']
    new_dir = data['DIR']
    
    cur.execute('SELECT HASH FROM ' + tablename + ' WHERE HASH = ?',(hash,))
    if (len(cur.fetchall()) == 0):
        cur.execute("INSERT INTO {} (HASH) VALUES (?)".format(tablename),(hash,))
        if(verbose): print(u'  \u251C\u2574'+'\033[1;32mInserting\033[1;0m: ' + new_dir)
        old_dir = None
    else:
        cur.execute('SELECT DIR FROM {} WHERE HASH = ?'.format(tablename),(hash,))
        old_dir = cur.fetchall()[0][0]
        
    for col in data:
        cur.execute('UPDATE {} SET {} = ? WHERE HASH = ?'.format(tablename,col),(data[col],hash))
    
    if old_dir:
        if (old_dir == new_dir):
            if (verbose): print(u'  \u251C\u2574'+'\033[1;33mUpdating\033[1;0m:  ' + new_dir + ' ( record already exists )')
        else:
            if (verbose): print(u'  \u251C\u2574'+'\033[1;36mMoving\033[1;0m:    ' + old_dir + ' --> ' + new_dir)

    if 'DIFF' in data:
        cur.execute("UPDATE " + tablename + " SET DIFF = ? WHERE HASH = ?",(data['DIFF'],hash))

def updateRegTestTable(cur,tablename,verbose=False):
    types = dict()
    types['RUN'] = 'VARCHAR(255)'
    types['TEST_NAME'] = 'VARCHAR(255)'
    types['RUNCODE'] = 'INT'
    types['COMPARE'] = 'VARCHAR(16)'
    types['PERFORMANCE'] = 'FLOAT'
    types['BENCHMARK_HASH'] = 'VARCHAR(24)'
    updateTable(cur,tablename,types,mode="regtest",verbose=False)

def updateRegTestRecord(cur,tablename,hash,run,test_name,status,benchmark_hash):
    cur.execute("INSERT INTO {} (HASH,RUN,TEST_NAME,RUNCODE,COMPARE,PERFORMANCE,BENCHMARK_HASH) VALUES (?,?,?,?,?,?,?)".format(tablename),
                (hash,run,test_name,status.runcode,status.compare,float(status.performance),benchmark_hash))