import sys
import os
sys.path.append(os.path.abspath('../../scripts/'))
import glob
import recurse


class HTMLPrinter:
    f = None
    intable = False
    myprefix = "     "
    def __init__(self,exe):
        self.f = open(f"InputIndex_{exe}.rst","w",encoding='utf-8')
        print(exe,file=self.f)
        print("--------------------------",file=self.f)
        print("",file=self.f)
        print("",file=self.f)
        print(".. raw:: html",file=self.f)
        print("",file=self.f)
    def __del__(self):
        self.f.close()
    def starttable(self):
        if not self.intable:
            print(self.myprefix,"<table class='api-inputs-table'>", file=self.f)
            print(self.myprefix,f"<thead><tr>", file=self.f)
            print(self.myprefix,f"<th>Input name</th>", file=self.f)
            print(self.myprefix,f"<th>Description</th>", file=self.f)
            print(self.myprefix,f"</tr></thead>", file=self.f)
        self.intable=True
    def endtable(self):
        if self.intable:
            print(self.myprefix,"</table><br/><br/>",file=self.f)
        self.intable=False
    def printinput(self,input,prefix,lev):
        if input['string']:
            name = f'.'.join(prefix + [input['string']])
            print(self.myprefix,"  "*lev,f"<tr>", file=self.f)
            print(self.myprefix,"  "*lev,f"  <td style='padding-left: {10*lev}px'>", file=self.f)
            print(self.myprefix,"  "*lev,f"    {name}", file=self.f)
            print(self.myprefix,"  "*lev,f"  </td>", file=self.f)
            print(self.myprefix,"  "*lev,f"  <td>", file=self.f)
            print(self.myprefix,"  "*lev,f"    {input['doc']}", file=self.f)
            print(self.myprefix,"  "*lev,f"  </td>", file=self.f)
            print(self.myprefix,"  "*lev,f"</tr>", file=self.f)
    def printconditional(self,inputname,inputvalue,lev):
        print(self.myprefix,"  "*lev,f"<tr>", file=self.f)
        print(self.myprefix,"  "*lev,f"  <td style='padding-left: {10*lev}px' colspan=2>", file=self.f)
        print(self.myprefix,"  "*lev,f"     <b> if </b>", file=self.f)
        print(self.myprefix,"  "*lev,f"     {inputname}", file=self.f)
        print(self.myprefix,"  "*lev,f"     = {inputvalue}", file=self.f)
        print(self.myprefix,"  "*lev,f"  </td>", file=self.f)
        print(self.myprefix,"  "*lev,f"</tr>", file=self.f)
    def printtablename(self,inputname,lev):
        sanitizedname = inputname.replace('<','&lt;').replace('>','&gt;')
        print(self.myprefix,"  "*lev,f"<h3> {sanitizedname}  </h3>", file=self.f)

for exe in ["alamo","mechanics","hydro","sfi","thermoelastic","topop"]:
    try:
        printer = HTMLPrinter(exe)
        recurse.recurse("../../src/",exe,printer=printer)
        del printer
    except Exception as e:
        del printer
        print("Problem writing file for ",exe)
        print(e)


