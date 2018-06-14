#!/usr/bin/python
import argparse
import sqlite3
import hashlib
import os

parser = argparse.ArgumentParser(description='Sift through outputs');
parser.add_argument('directories', nargs='*', help='List of directories containing ALAMO output');
parser.add_argument('-d','--database', default='results.db', help='Name of database');
parser.add_argument('-t','--table', default='simulation_data', help='Table name in database');

args=parser.parse_args();

db = sqlite3.connect(args.database if args.database.endswith('.db') else args.database+'.db')
cur= db.cursor()

types = dict()

for directory in args.directories:
    f = open(directory+'/metadata')
    for line in f.readlines():
        if line.startswith('#'): continue;
        if '::' in line: continue; # TODO: need to replace this to be more robust!
        if len(line.split(' = ')) != 2: continue;
        key = line.split(' = ')[0].replace('.','_')
        val = line.split(' = ')[1].replace('  ','')
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

cur.execute('CREATE TABLE IF NOT EXISTS ' + args.table + '('
            'HASH VARCHAR(255) UNIQUE, ' +
            'ID VARCHAR(255) UNIQUE,' +
            'Description VARCHAR(8000),' +
            'Tags VARCHAR(1000),'
            +','.join([key+' '+types[key] for key in types]) + ');')


for directory in args.directories:
    sim_hash = hashlib.sha224(directory).hexdigest()
    f = open(directory+'/metadata')
    cols = ['HASH','ID']
    vals = ['"'+sim_hash+'"', '"'+os.path.abspath(directory)+'"']
    for line in f.readlines():
        if line.startswith('#'): continue;
        if '::' in line: continue; # TODO: need to replace this to be more robust!
        if len(line.split(' = ')) != 2: continue;
        cols.append(line.split(' = ')[0].replace('.','_'))
        vals.append('"' + line.split(' = ')[1].replace('  ','').replace('\n','') + '"')

    cur.execute('INSERT OR IGNORE INTO ' + args.table + ' (' + ','.join(cols) + ') ' +
                'VALUES (' + ','.join(vals) + ')')

    print('Added/updated entry ' + directory)
db.commit()
db.close()
