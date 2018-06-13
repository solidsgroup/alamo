#!/usr/bin/python
import os
import glob
import fnmatch
import sqlite3
from flask import Flask, render_template, send_file

print("==================================")
print("SIMBA: SIMulation Browser Analysis")
print("==================================")



script_directory = os.path.realpath(__file__)

app = Flask(__name__)

@app.route("/")
def go():
    db = sqlite3.connect('results.db')
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT * FROM " + tables[0] )
    data = cur.fetchall()

    cur.execute("PRAGMA table_info("+tables[0]+")")
    columns=[a[1] for a in cur.fetchall()]

    db.commit()
    db.close()
    
    return render_template('template.html',
                           table_name=tables[0],
                           table=data,
                           columns=columns)

imgfiles = []

def find_images(path):
    global imgfiles
    img_fmts = ['.jpg', '.jpeg', '.png', '.pdf']
    imgfiles = []
    for fmt in img_fmts: imgfiles += glob.glob(path+'/*'+fmt)

@app.route('/img/<number>')
def serve_image(number):
    global imgfiles
    return send_file(imgfiles[int(number)])

@app.route('/table/<table>/entry/<entry>')
def database(table,entry):
    global imgfiles

    db = sqlite3.connect('results.db')
    db.text_factory = str
    cur= db.cursor()

    cur.execute("PRAGMA table_info("+table+")")
    columns=[a[1] for a in cur.fetchall()]

    cur.execute("SELECT * FROM " + table + " WHERE HASH='" + entry + "'")
    data = cur.fetchall()[0]

    find_images(data[1])

    return render_template('detail.html',
                           columns=columns,
                           data=data,
                           imgfiles=imgfiles)

if __name__ == '__main__':
    app.run(debug=True,use_reloader=False)
