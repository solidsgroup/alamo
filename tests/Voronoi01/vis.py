import argparse

parser = argparse.ArgumentParser(description='Start a webserver to brows database entries');
parser.add_argument('directories', nargs='*', help='List of directories containing ALAMO output');
args=parser.parse_args();
print(args.directories)
for d in args.directories:
    print("Working on directory " + d)
    DeleteAllPlots()
    OpenDatabase(d+"/output.visit", 0)
    AddPlot("Pseudocolor", "Eta001", 1, 1)
    AddPlot("Mesh", "Mesh", 1, 1)
    DrawPlots()
    TimeSliderPreviousState()
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 0
    SaveWindowAtts.outputDirectory = d
    SaveWindowAtts.fileName = "visit"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
exit(0)
