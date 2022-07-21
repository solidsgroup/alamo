import os
import glob

for f in glob.glob("output*"):
    if not os.path.isfile("{}/nodeoutput.visit".format(f)):
        print(f," is empty")
        os.system("rm -rf {}".format(f))
                   

