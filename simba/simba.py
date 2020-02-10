import os
import re
import sqlite3



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
        if (key == 'DIFF'):
            types['DIFF'] = 'BLOB'
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
def updateTable(cur,tablename,types):
    # If the table does not exist, create it
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cur.fetchall()]
    if "__tables__" not in tables:
        cur.execute('CREATE TABLE __tables__ ('
                    'NAME UNIQUE,' +
                    'Description VARCHAR(8000));')
    if tablename not in tables:
        cur.execute('CREATE TABLE ' + tablename + ' ('
                    'HASH VARCHAR(255) UNIQUE, ' +
                    'DIR,' +
                    'Description VARCHAR(8000),' +
                    'Tags VARCHAR(1000)' +
                    (',' if len(types)>0 else '') +
                    ','.join([key+' '+types[key] for key in sorted(types)]) +
                    ');')
        print('\033[1;32mADDED TABLE\033[1;0m: ' + tablename)
    else:
        print('\033[1;34mUSING TABLE\033[1;0m: ' + tablename)

    cur.execute('SELECT * FROM __tables__ WHERE NAME = \"' + tablename + '\";')
    if (len(cur.fetchall()) == 0): cur.execute('INSERT INTO __tables__ (NAME, Description) VALUES (\"' + tablename + '\", \"Description\");')

    # If the table exists, but new columns have been added, modify the table
    # accordingly
    cur.execute("PRAGMA table_info("+tablename+")")
    columns=[a[1] for a in cur.fetchall()]
    print(types)
    print(columns)
    for key in types:
        if key not in columns:
            cur.execute('ALTER TABLE ' + tablename + ' ADD ' + key + ' ' + types[key])
            print('\033[1;34mADDED COLUMN\033[1;0m: ' + key + ' to ' + tablename)


def updateRecord(cur,tablename,data,hash=None):
    if not hash: hash = data['HASH']
    new_dir = data['DIR']
    
    cur.execute('SELECT HASH FROM ' + tablename + ' WHERE HASH = ?',(hash,))
    if (len(cur.fetchall()) == 0):
        cur.execute("INSERT INTO {} (HASH) VALUES (?)".format(tablename),(hash,))
        print(u'  \u251C\u2574'+'\033[1;32mInserting\033[1;0m: ' + new_dir)
        old_dir = None
    else:
        cur.execute('SELECT DIR FROM {} WHERE HASH = ?'.format(tablename),(hash,))
        old_dir = cur.fetchall()[0][0]
        
    for col in data:
        cur.execute('UPDATE {} SET {} = ? WHERE HASH = ?'.format(tablename,col),(data[col],hash))
    
    # Are we merely updating and changing the home location of a record?
    # If so, let's let the user know.
    if old_dir:
        if (old_dir == new_dir):
            print(u'  \u251C\u2574'+'\033[1;33mUpdating\033[1;0m:  ' + new_dir + ' ( record already exists )')
        else:
            print(u'  \u251C\u2574'+'\033[1;36mMoving\033[1;0m:    ' + old_dir + ' --> ' + new_dir)

    #else:
    #    cur.execute('INSERT INTO ' + tablename + ' (' + ','.join(data.keys()) + ') ' +
    #                'VALUES (' + ','.join(data.values()) + ')')
    #    print(u'  \u251C\u2574'+'\033[1;32mInserting\033[1;0m: ' + new_dir)
    
    if 'DIFF' in data:
        cur.execute("UPDATE " + tablename + " SET DIFF = ? WHERE HASH = ?",(data['DIFF'],hash))
