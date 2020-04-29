#!/usr/bin/env python3
import argparse
import sqlite3
import hashlib
import os
import re
import simba


parser = argparse.ArgumentParser(description='Sift through outputs')
parser.add_argument('directories', nargs='*', help='List of directories containing ALAMO output')
parser.add_argument('-d','--database', default='results.db', help='Name of database')
parser.add_argument('-r','--remove', nargs='*', help='Tables to remove')
parser.add_argument('-t','--table', default='simulation_data', help='Table name in database')

args=parser.parse_args()

db = sqlite3.connect(args.database if args.database.endswith('.db') else args.database+'.db')
db.text_factory = str
cur= db.cursor()
types = dict()

#
# Scan metadata files to determine columns
#
for directory in args.directories:
    if not os.path.isfile(directory+'/metadata'):
        continue

    data = simba.parse(directory)
    types.update(simba.getTypes(data))

#
# Update/create the chosen table
#
simba.updateTable(cur,args.table,types)

#
# If there are tables to delete, delete them
#
if args.remove:
    for tab in list(args.remove):
        cur.execute('DROP TABLE ' + tab)


#
# Scan each metadata file and add an entry to the table, skipping any
# records that already exist.
#
for directory in args.directories:
    if not os.path.isfile(directory+'/metadata'):
        print(u'  \u251C\u2574\033[1;31mSkipping\033[1;0m:  ' + directory + ' (no metadata) ')
        continue

    data = simba.parse(directory)
    data['DIR'] = os.path.abspath(directory)
    if 'HASH' in data:
        sim_hash = data['HASH']
    else:
        f = open(directory+'/metadata')
        sim_hash = str(hashlib.sha224(f.read().encode('utf-8')).hexdigest())
        f.close()
        data['HASH'] = sim_hash
    
    simba.updateRecord(cur,args.table,data,sim_hash)

        
print(u'  \u2514\u2574' + 'Done')


db.commit()
db.close()
