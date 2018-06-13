#!/usr/bin/python
import os
import sqlite3
from flask import Flask, render_template

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

    return render_template('template.html', my_string="Wheeeee!", tables=tables, my_list=tables, table=data, columns=columns)

@app.route('/entry/<code>')
def database(code):
    return "Hello world, " + code

if __name__ == '__main__':
    app.run(debug=True,use_reloader=False)
