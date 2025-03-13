import sys
import os
sys.path.append(os.path.abspath('../../scripts/'))
import glob
import recurse

f = open("InputIndex.html","w")


for line in sorted(glob.glob("../../src/*.cc")):
    print(">>>",line)
    cc = line.replace("../../src/","").replace(".cc","")
    print(">>>",cc)
    print(f"<h1>{cc}</h1>",file=f)
    recurse.recurse("../../src/","alamo",f=f)

f.close()


