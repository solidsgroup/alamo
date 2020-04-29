#!/usr/bin/env python3
import os
import glob
import fnmatch
import sqlite3
import argparse
import getpass
from functools import wraps
from flask import Flask, request, render_template, send_file, redirect, Response
from flaskext.markdown import Markdown
import datetime

print("====================================")
print("SIMBA: SIMulation Browser Analysis")
print("====================================")

parser = argparse.ArgumentParser(description='Start a webserver to brows database entries');
parser.add_argument('-i','--ip', default='127.0.0.1', help='IP address of server (default: localhost)');
parser.add_argument('-p','--port', default='5000', help='Port (default: 5000)');
parser.add_argument('-d','--database',default='results.db',help='Name of database to read from')
parser.add_argument('-s','--safe',dest='safe',action='store_true',help='Safe mode - disallow permanent record deletion')
parser.add_argument('-f','--fast',dest='fast',action='store_true',help='Fast mode - fewer features for working with large datasets')
parser.add_argument('--pwd',default=False,action='store_true')
args=parser.parse_args()

pwd = None
usr = None
if args.pwd:
    usr = getpass.getuser()
    pwd = getpass.getpass()

def check_auth(username, password):
    """This function is called to check if a username /
    password combination is valid.
    """
    return username == usr and password == pwd

def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response(
    'Could not verify your access level for that URL.\n'
    'You have to login with proper credentials', 401,
    {'WWW-Authenticate': 'Basic realm="Login Required"'})

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if (pwd and usr):
            if not auth or not check_auth(auth.username, auth.password):
                return authenticate()
        return f(*args, **kwargs)
    return decorated    

def format_datetime(value, format='medium'):
    #t = parser.parse(value)
    #datetime.
    if not value: return value
    
    for fmt in ["%a %b%d %H:%M:%S %Y", "%a %b %d %H:%M:%S %Y"]:
        try:
            dt = datetime.datetime.strptime(str(value),fmt)
            return(dt.strftime("%Y-%m-%d %H:%M:%S (%a)"))
        except ValueError:
            pass
    print("Date parsing failed for string " + value + "")
    return value


if not args.safe and not args.ip == '127.0.0.1' or args.ip == 'localhost':
    print("=============  WARNING =============")
    print("It appears that you are starting    ")
    print("SIMBA on a public server NOT in SAFE")
    print("mode. This could allow malicious    ")
    print("users to alter your records. It is  ")
    print("strongly recommended that you run   ")
    print("with the --safe or -s flags enabled!")
    print("====================================")

script_directory = os.path.realpath(__file__)

app = Flask(__name__)
Markdown(app)

app.jinja_env.filters['datetime'] = format_datetime

@app.route("/", methods=['GET','POST'])
@requires_auth
def root():
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()

    if request.method == 'POST':
        if request.form.get('action')=="delete-table" and not args.safe:
            print("deleting table " + request.form.get('table-name'))
            cur.execute("DROP TABLE " + request.form.get('table-name'))
        if request.form.get('action')=="rename-table" and not args.safe:
            print(request.form.get('table-name-old'))
            print(request.form.get('table-name-new'))
            cur.execute("ALTER TABLE " + str(request.form.get('table-name-old')) + " RENAME TO " + str(request.form.get('table-name-new')))
            

    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cur.fetchall()]
    if "__tables__" in tables: tables.remove("__tables__")


    if len(tables) > 0:
        return redirect('/table/'+tables[0])
    return render_template('root.html',
                           tables=tables)

@app.route("/table/<table>", methods=['GET','POST'])
@requires_auth
def table(table):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()

    if request.method == 'POST':
        if request.form.get('table-description') and not args.safe:
            print(request.form.get('table-description'))
            cur.execute("UPDATE __tables__ SET Description = ? WHERE NAME = ?;", (request.form.get('table-description'), table))
        items = request.form.items()
        for f in items:
            print(f)
            if str(f[0]).startswith('hash_'):
                hash = str(f[0]).replace('hash_','')
                dir = str(f[1])
                print("DELETING ",hash)
                cur.execute("DELETE FROM " + table + " WHERE HASH = ?;",(hash,))
                print("DELETING ",dir)
                os.system('rm -rf ' + dir)



    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in sorted(cur.fetchall())]
    if "__tables__" in tables: tables.remove("__tables__")

    if not table: table_name = tables[0]
    else: table_name = table

    if request.method == 'POST':
        if request.form.get('action')=="delete-entry-only" and not args.safe:
            cur.execute("DELETE FROM " + table + " WHERE HASH = ?;",(request.form.get('entry-hash'),))
        if request.form.get('action')=='delete-everything' and not args.safe:
            cur.execute("SELECT DIR FROM " + table + " WHERE HASH = ?",(request.form.get('entry-hash'),))
            os.system('rm -rf ' + cur.fetchall()[0][0])
            cur.execute("DELETE FROM " + table + " WHERE HASH = ?;",(request.form.get('entry-hash'),))
    

    cur.execute("PRAGMA table_info("+table_name+")")
    columns=[a[1] for a in cur.fetchall()]

    cur.execute("SELECT * FROM " + table_name )
    rawdata = cur.fetchall()

    data = []
    for d in rawdata: data.append(dict(zip(columns,d)))


    cur.execute("SELECT Description FROM __tables__ WHERE Name = \"" + table_name  + "\"")
    desc = list(cur.fetchall()[0])[0]


    status = []
    if ("Status" in columns):
        status = [d["Status"] for d in data]

    db.commit()
    db.close()
    
    if table==None or table not in tables: table = tables[0];

    numfiles = []
    if not args.fast:
        for d in data:
            find_images(d['DIR'])
            numfiles.append(len(imgfiles))

    columns.insert(0,columns.pop(columns.index('DIR')))
    columns.insert(1,columns.pop(columns.index('Description')))
    columns.insert(1,columns.pop(columns.index('Tags')))


    return render_template('template.html',
                           tables=tables,
                           table_name=table,
                           table_description=desc,
                           data=data,
                           status=status,
                           numfiles=numfiles,
                           columns=columns)
imgfiles = []

def find_images(path):
    print("Path is ",path)
    global imgfiles
    img_fmts = ['.jpg', '.jpeg', '.png', '.gif','.svg']
    imgfiles = []
    for fmt in img_fmts: imgfiles += glob.glob(path+'/*'+fmt)
    imgfiles.sort()

def find_tarballs(path):
    global tarballfiles
    img_fmts = ['.tar.gz']
    tarballfiles = []
    for fmt in img_fmts: tarballfiles += glob.glob(path+'/*'+fmt)
    tarballfiles.sort()

def find_thermo(path):
    global thermofile
    if os.path.isfile(path+'/thermo.dat'): thermofile = path+"/thermo.dat"
    else: thermofile = None

@app.route('/img/<number>')
@requires_auth
def serve_image(number):
    global imgfiles
    return send_file(imgfiles[int(number)],cache_timeout=-1)

@app.route('/metadata/')
@requires_auth
def serve_metadata():
    global metadatafile
    response = send_file(metadatafile,cache_timeout=-1,as_attachment=True)
    response.headers["x-filename"] = "metadata"
    response.headers["Access-Control-Expose-Headers"] = 'x-filename'
    return response


@app.route('/thermo/')
@requires_auth
def serve_thermo():
    global thermofile
    response = send_file(thermofile,cache_timeout=-1,as_attachment=True)
    response.headers["x-filename"] = "thermo.dat"
    response.headers["Access-Control-Expose-Headers"] = 'x-filename'
    return response

@app.route('/tarball/<filename>/<number>')
@requires_auth
def serve_tarball(filename,number):
    print (filename)
    global tarballfiles
    response = send_file(tarballfiles[int(number)],cache_timeout=-1,as_attachment=True)
    response.headers["x-filename"] = filename
    response.headers["Access-Control-Expose-Headers"] = 'x-filename'
    return response

@app.route('/table/<table>/entry/<entry>', methods=['GET','POST'])
@requires_auth
def table_entry(table,entry):

    global imgfiles
    global metadatafile
    global thermofile

    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    
    if request.method == 'POST':
        if request.form.get('description'):
            cur.execute("UPDATE " + table + " SET Description = ? WHERE HASH = ?;",
                        (request.form.get('description'), entry));
        if request.form.get('tags'):
            cur.execute("UPDATE " + table + " SET Tags = ? WHERE HASH = ?;",
                        (request.form.get('tags'), entry));

    cur.execute("PRAGMA table_info("+table+")")
    columns=[a[1] for a in cur.fetchall()]

    cur.execute("SELECT * FROM " + table + " WHERE HASH='" + entry + "'")
    d = cur.fetchall()[0]

    data = dict(zip(columns,d))

    find_images(data['DIR'])
    find_tarballs(data['DIR'])
    metadatafile=data['DIR']+"/metadata"
    find_thermo(data['DIR'])

    db.commit()
    db.close()
    
    columns.insert(0,columns.pop(columns.index('DIR')))
    columns.insert(1,columns.pop(columns.index('Description')))
    columns.insert(1,columns.pop(columns.index('Tags')))
    
    return render_template('detail.html',
                           table=table,
                           entry=entry,
                           columns=columns,
                           data=data,
                           thermofile=thermofile,
                           imgfiles=[os.path.split(im)[1] for im in imgfiles],
                           tarballfiles=[os.path.split(tb)[1] for tb in tarballfiles])
                           
@app.route('/table/<table>/entry/<entry>/stdout', methods=['GET','POST'])
@requires_auth
def table_entry_stdout(table,entry):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT STDOUT FROM {} WHERE HASH = ?".format(table),(entry,))
    return cur.fetchall()[0][0]

@app.route('/table/<table>/entry/<entry>/diff', methods=['GET','POST'])
@requires_auth
def table_entry_diff(table,entry):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT DIFF FROM {} WHERE HASH = ?".format(table),(entry,))
    return cur.fetchall()[0][0]

@app.route('/table/<table>/entry/<entry>/diff.patch')
@requires_auth
def table_entry_diff_patch(table,entry):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT DIFF_PATCH FROM {} WHERE HASH = ?".format(table),(entry,))
    return Response(cur.fetchall()[0][0],content_type='File')
    #return cur.fetchall()[0][0]

@app.route('/table/<table>/entry/<entry>/stderr', methods=['GET','POST'])
@requires_auth
def table_entry_stderr(table,entry):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT STDERR FROM {} WHERE HASH = ?".format(table),(entry,))
    return cur.fetchall()[0][0]
                           
@app.route('/regtest/<regtest>/<run>/stdout', methods=['GET','POST'])
@requires_auth
def regtest_run_stdout(regtest,run):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()
    cur.execute("SELECT STDIO FROM regtest_runs WHERE RUN = ?",(run,))
    ret = cur.fetchall()
    if len(ret) > 0:
        if len(ret[0]) > 0:
            return ret[0][0]
    else: return "<h2>None</h2>"


@app.route("/regtest/<regtest>", methods=['GET','POST'])
@requires_auth
def regtest(regtest):
    db = sqlite3.connect(args.database)
    db.text_factory = str
    cur= db.cursor()

    cur.execute("SELECT DISTINCT TEST_NAME FROM {}".format(regtest))
    test_names = [tn[0] for tn in cur.fetchall()]

    cur.execute("SELECT RUN,COMPILECODE FROM regtest_runs ORDER BY RUN DESC".format(regtest))
    runs = cur.fetchall()#sorted([tn[0] for tn in cur.fetchall()],reverse=True)
        
    cur.execute("PRAGMA table_info({})".format(regtest))
    columns=[a[1] for a in cur.fetchall()]

    cur.execute("SELECT * FROM {}".format(regtest))
    rawdata = cur.fetchall()


    data = []
    for d in rawdata: data.append(dict(zip(columns,d)))

    return render_template('regtest.html',
                            runs=runs,
                            tests=test_names,
                            data=data,
                            columns=columns)
                           
                           

if __name__ == '__main__':
    app.run(debug=True,
            use_reloader=False,
            host=args.ip,
            port=int(args.port))
