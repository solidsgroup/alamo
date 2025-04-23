from glob import glob


src2url = {}

try:
	doxsourcefiles = sorted(glob("build/doxygen/*source.html"))
	for doxsourcefile in doxsourcefiles:
	    print(doxsourcefile)
	    with open(doxsourcefile) as f:
	        for line in f.readlines():
	            if r"<title>" in line:
	                line = line.replace(r"<title>Alamo: ","")
	                line = line.replace(r" Source File</title>","")
	                line = line.replace("\n","")
	                src2url[line] = doxsourcefile
	                continue
except Exception as e:
    print(e)
print(src2url)


