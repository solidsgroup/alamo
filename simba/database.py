#!/usr/bin/env python3
import argparse
import sqlite3
import hashlib
import os
import re


parser = argparse.ArgumentParser(description='Sift through outputs');
parser.add_argument('directories', nargs='*', help='List of directories containing ALAMO output');
parser.add_argument('-d','--database', default='results.db', help='Name of database');
parser.add_argument('-r','--remove', nargs='*', help='Tables to remove');
parser.add_argument('-t','--table', default='simulation_data', help='Table name in database');

args=parser.parse_args();

db = sqlite3.connect(args.database if args.database.endswith('.db') else args.database+'.db')
db.text_factory = str
cur= db.cursor()
types = dict()

def parse(filename):
    f = open(filename)

    things = dict()

    for line in f.readlines():
        if line.startswith('#'): continue;
        if '::' in line:
            line = re.sub(r'\([^)]*\)', '',line)
            line = line.replace(" :: ", " = ").replace('[','').replace(',','').replace(']','').replace(' ','')
        if len(line.split(' = ')) != 2: continue;
        col = line.split(' = ')[0].replace('.','_')
        val = line.split(' = ')[1].replace('  ','').replace('\n','').replace(';','')
        
        things[col] = val

    return things
    

#
# Scan metadata files to determine columns
#
for directory in args.directories:
    if not os.path.isfile(directory+'/metadata'):
        continue

    data = parse(directory+'/metadata')

    for key in data:
        val = data[key]
        if (key == 'HASH'): continue
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

#
# If the table does not exist, create it
#
cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = [r[0] for r in cur.fetchall()]
if "__tables__" not in tables:
    cur.execute('CREATE TABLE __tables__ ('
#                'HASH VARCHAR(255) UNIQUE, ' +
                'NAME UNIQUE,' +
                'Description VARCHAR(8000));')
if args.table not in tables:
    cur.execute('CREATE TABLE ' + args.table + ' ('
                'HASH VARCHAR(255) UNIQUE, ' +
                'DIR,' +
                'Description VARCHAR(8000),' +
                'Tags VARCHAR(1000)' +
                (',' if len(types)>0 else '') +
                ','.join([key+' '+types[key] for key in sorted(types)]) +
                ');')
    print('\033[1;32mADDED TABLE\033[1;0m: ' + args.table)
else:
    print('\033[1;34mUSING TABLE\033[1;0m: ' + args.table)

#cur.execute('DROP TABLE __tables__')
cur.execute('SELECT * FROM __tables__ WHERE NAME = \"' + args.table + '\";')
if (len(cur.fetchall()) == 0): cur.execute('INSERT INTO __tables__ (NAME, Description) VALUES (\"' + args.table + '\", \"Description\");')

#
# If there are tables to delete, delete them
#
if args.remove:
    for tab in list(args.remove):
        cur.execute('DROP TABLE ' + tab)

#
# If the table exists, but new columns have been added, modify the table
# accordingly
#
cur.execute("PRAGMA table_info("+args.table+")")
columns=[a[1] for a in cur.fetchall()]
for key in types:
    if key not in columns:
        cur.execute('ALTER TABLE ' + args.table + ' ADD ' + key + ' ' + types[key])
        print('\033[1;34mADDED COLUMN\033[1;0m: ' + key + ' to ' + args.table)

#
# Scan each metadata file and add an entry to the table, skipping any
# records that already exist.
#
for directory in args.directories:
    if not os.path.isfile(directory+'/metadata'):
        print(u'  \u251C\u2574\033[1;31mSkipping\033[1;0m:  ' + directory + ' (no metadata) ')
        continue

    data = parse(directory+'/metadata')
    data['DIR'] = os.path.abspath(directory);
    if 'HASH' in data:
        sim_hash = data['HASH']
    else:
        f = open(directory+'/metadata')
        sim_hash = str(hashlib.sha224(f.read().encode('utf-8')).hexdigest())
        f.close()
        data['HASH'] = sim_hash
    
    
    cur.execute('SELECT HASH FROM ' + args.table + ' WHERE HASH = ?',(sim_hash,))
    if (len(cur.fetchall()) != 0):
        cur.execute('SELECT DIR FROM ' + args.table + ' WHERE HASH = ?',(sim_hash,))
        old_dir = cur.fetchall()[0][0]
        new_dir = os.path.abspath(directory)
        cur.execute('UPDATE ' + args.table +
                    ' SET ' + ','.join([col+' = ' + '"' + data[col] + '"' for col in data]) +
                    ' WHERE HASH = ' + '"' + sim_hash + '"')
        if (old_dir == new_dir):
            print(u'  \u251C\u2574'+'\033[1;33mUpdating\033[1;0m:  ' + directory + ' ( record already exists )')
        else:
            print(u'  \u251C\u2574'+'\033[1;36mMoving\033[1;0m:    ' + old_dir + ' --> ' + new_dir)
    else:
        cur.execute('INSERT INTO ' + args.table + ' (' + ','.join([c for c in data]) + ') ' +
                    'VALUES (' + ','.join(['"'+data[c]+'"' for c in data]) + ')')
        print(u'  \u251C\u2574'+'\033[1;32mInserting\033[1;0m: ' + directory)
print(u'  \u2514\u2574' + 'Done')


db.commit()
db.close()
