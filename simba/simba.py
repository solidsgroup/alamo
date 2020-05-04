import os
import re
import sqlite3


class Status:
    runcode = -1
    compare = "NO"
    diff_stdout = ""
    runtime = 0.0
    bm_runtime = 0.0

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

#
# Scan through an entry's data dictionary and determine datatypes
#
def getTypes(data):
    types = dict()
    for key in data:
        val = data[key]
        if (key == 'HASH'): 
            continue
        if (key == 'DIFF' or key == 'DIFF_PATCH' or key == 'STDOUT' or key == 'STDERR'):
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
                        'Description VARCHAR(8000),' +
                        'Tags VARCHAR(1000)' +
                        (',' if len(types)>0 else '') +
                        ','.join([key+' '+types[key] for key in sorted(types)]) +
                        ');')
        elif mode == "regtest":
            cur.execute('CREATE TABLE regtest ('
                        'HASH VARCHAR(255) UNIQUE ' +
                        (',' if len(types)>0 else '') +
                        ','.join([key+' '+types[key] for key in sorted(types)]) +
                        ');')
        elif mode == "regtest_runs":
            cur.execute('CREATE TABLE regtest_runs ('
                        'RUN VARCHAR(255) UNIQUE ' +
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
    sqlstring = "PRAGMA table_info("+tablename+")"
    print(sqlstring)
    cur.execute(sqlstring)
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

def updateRegTestTable(cur,verbose=False):
    types = dict()
    types['RUN'] = 'VARCHAR(255)'
    types['DIR'] = 'VARCHAR(255)'
    types['TEST_NAME'] = 'VARCHAR(255)'
    types['RUNCODE'] = 'INT'
    types['COMPARE'] = 'VARCHAR(16)'
    types['DIFF_STDOUT'] = 'BLOB'
    types['RUNTIME'] = 'FLOAT'
    types['BM_RUNTIME'] = 'FLOAT'
    types['BENCHMARK_HASH'] = 'VARCHAR(24)'
    types['BENCHMARK_RUN'] = 'VARCHAR(64)'
    updateTable(cur,"regtest",types,mode="regtest",verbose=False)
    types = dict()
    types['COMPILECODE'] = 'INT'
    types['STDIO']       = 'BLOB'
    updateTable(cur,"regtest_runs",types,mode="regtest_runs",verbose=False)


def updateRegTestRecord(cur,hash,run,test_name,status,benchmark_hash,benchmark_run,rt_dir):
    cur.execute("INSERT INTO regtest (HASH,DIR,RUN,TEST_NAME,RUNCODE,COMPARE,DIFF_STDOUT,RUNTIME,BM_RUNTIME,BENCHMARK_HASH,BENCHMARK_RUN) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (hash,rt_dir,run,test_name,status.runcode,status.compare,status.diff_stdout,float(status.runtime),float(status.bm_runtime),benchmark_hash,benchmark_run))

def updateRegTestRun(cur,run,compilecode,stdio):
    cur.execute('SELECT RUN FROM regtest_runs WHERE RUN = ?',(run,))
    if (len(cur.fetchall()) == 0):
        cur.execute("INSERT INTO regtest_runs (RUN,COMPILECODE,STDIO) VALUES (?,?,?)",
                    (run,compilecode,stdio))
    else:
        cur.execute("UPDATE regtest_runs SET COMPILECODE=?, STDIO=? WHERE RUN=?",(compilecode,stdio,run))
    

