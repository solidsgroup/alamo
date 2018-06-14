# Visit 2.13.0 log file
ScriptVersion = "2.13.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
SetActivePlots()
DeleteActivePlots()
OpenDatabase("output/output.visit", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.00111985207/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.00111985207"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.00395226479/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.00395226479"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.01753020287/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.01753020287"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.03576922417/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.03576922417"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.03660058975/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.03660058975"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.06631302834/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.06631302834"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.06679081917/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.06679081917"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.06885242462/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.06885242462"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.07678103447/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.07678103447"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.07860279083/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.07860279083"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08126592636/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08126592636"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08346748352/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08346748352"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08403873444/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08403873444"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08755111694/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08755111694"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08787131310/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08787131310"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.08919143677/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.08919143677"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.09110426903/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.09110426903"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.09401273727/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.09401273727"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.09616994858/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.09616994858"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.09862065315/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.09862065315"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10053396225/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10053396225"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10074973106/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10074973106"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10091710091/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10091710091"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10485672951/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10485672951"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10844516754/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10844516754"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.10870289803/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.10870289803"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.11164474487/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.11164474487"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.11395835876/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.11395835876"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.12316799164/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.12316799164"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.12597012520/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.12597012520"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.13473391533/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.13473391533"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.14094591141/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.14094591141"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.14987516403/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.14987516403"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.16405677795/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.16405677795"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.16424560547/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.16424560547"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.16503715515/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.16503715515"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.16637539864/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.16637539864"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.16916584969/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.16916584969"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.17586612701/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.17586612701"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.17758989334/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.17758989334"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.18354296684/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.18354296684"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.18471193314/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.18471193314"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.18913912773/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.18913912773"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.20695734024/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.20695734024"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.20762872696/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.20762872696"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.22675728798/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.22675728798"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.23051238060/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.23051238060"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.23457622528/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.23457622528"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.23953843117/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.23953843117"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.25177335739/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.25177335739"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.25697278976/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.25697278976"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.27137351036/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.27137351036"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.27303266525/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.27303266525"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.28048253059/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.28048253059"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.28785395622/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.28785395622"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.28873920441/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.28873920441"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.29806518555/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.29806518555"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.30255174637/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.30255174637"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.30410504341/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.30410504341"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.30487036705/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.30487036705"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.31051993370/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.31051993370"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.31095004082/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.31095004082"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.32996630669/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.32996630669"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.33450508118/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.33450508118"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.33487081528/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.33487081528"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.34812974930/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.34812974930"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.34970140457/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.34970140457"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.36136317253/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.36136317253"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.36588644981/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.36588644981"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.37053894997/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.37053894997"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.37120056152/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.37120056152"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.37454700470/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.37454700470"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.37496829033/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.37496829033"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.38452482224/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.38452482224"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.39749050140/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.39749050140"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.39906167984/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.39906167984"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42121243477/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42121243477"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42123627663/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42123627663"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42280220985/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42280220985"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42469072342/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42469072342"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42567610741/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42567610741"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.42730355263/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.42730355263"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.43358898163/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.43358898163"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.43368387222/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.43368387222"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.43675160408/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.43675160408"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.44686365128/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.44686365128"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.44700431824/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.44700431824"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.44712066650/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.44712066650"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.44964003563/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.44964003563"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.47158122063/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.47158122063"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.49078273773/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.49078273773"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.49622058868/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.49622058868"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.49828600883/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.49828600883"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.51029920578/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.51029920578"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.52233219147/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.52233219147"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.52489662170/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.52489662170"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.52649116516/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.52649116516"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.53445911407/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.53445911407"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.53476953506/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.53476953506"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.53539466858/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.53539466858"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.54520153999/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.54520153999"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.54818105698/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.54818105698"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.54943442345/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.54943442345"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.54998278618/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.54998278618"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.55299234390/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.55299234390"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.55378103256/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.55378103256"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.55734372139/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.55734372139"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.56414556503/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.56414556503"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.56872248650/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.56872248650"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.57159900665/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.57159900665"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots()
DeleteActivePlots()
OpenDatabase("output.old.57614469528/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.57614469528"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.58278083801/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.58278083801"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.58391523361/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.58391523361"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.60699772835/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.60699772835"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.61032056808/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.61032056808"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.61614322662/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.61614322662"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62041568756/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62041568756"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62175416946/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62175416946"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62577509880/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62577509880"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62680125237/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62680125237"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62684893608/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62684893608"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.62688493729/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.62688493729"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.63363289833/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.63363289833"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.64436030388/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.64436030388"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.64905524254/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.64905524254"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.65138745308/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.65138745308"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.66375184059/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.66375184059"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.66738557816/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.66738557816"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.67322921753/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.67322921753"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.67684221268/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.67684221268"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.68038439751/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.68038439751"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.68394041061/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.68394041061"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.69156384468/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.69156384468"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.69669842720/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.69669842720"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.69760298729/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.69760298729"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.70088005066/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.70088005066"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.70178127289/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.70178127289"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.70268487930/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.70268487930"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.71367859840/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.71367859840"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.71559786797/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.71559786797"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.72152972221/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.72152972221"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.72179341316/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.72179341316"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.72618937492/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.72618937492"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73087120056/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73087120056"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73348784447/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73348784447"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73557472229/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73557472229"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73575997353/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73575997353"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73668718338/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73668718338"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73861956596/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73861956596"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.73941254616/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.73941254616"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.74150371552/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.74150371552"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.74229502678/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.74229502678"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.74466943741/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.74466943741"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.75484228134/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.75484228134"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.75658345222/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.75658345222"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.76499176025/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.76499176025"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.76504015923/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.76504015923"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.76642012596/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.76642012596"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.76777505875/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.76777505875"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.77009034157/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.77009034157"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.77735972404/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.77735972404"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.78349804878/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.78349804878"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.78809499741/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.78809499741"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.78926968575/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.78926968575"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.78941559792/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.78941559792"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.80025458336/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.80025458336"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.80205655098/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.80205655098"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.80463528633/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.80463528633"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.80608129501/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.80608129501"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.80889153481/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.80889153481"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.82767581940/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.82767581940"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.83703851700/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.83703851700"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.83733773232/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.83733773232"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.83801651001/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.83801651001"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.84492921829/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.84492921829"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.84630489349/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.84630489349"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.85292577744/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.85292577744"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.85481953621/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.85481953621"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.85641431808/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.85641431808"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.85773706436/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.85773706436"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.86719560623/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.86719560623"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.88202500343/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.88202500343"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.89074444771/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.89074444771"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.89151859283/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.89151859283"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.89901494980/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.89901494980"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.90384697914/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.90384697914"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.91760706902/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.91760706902"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.93431949615/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.93431949615"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.93595528603/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.93595528603"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.94814324379/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.94814324379"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.94992494583/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.94992494583"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.95615077019/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.95615077019"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.95809674263/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.95809674263"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.97038245201/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.97038245201"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.97136998177/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.97136998177"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.97348880768/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.97348880768"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.98047614098/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.98047614098"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.98392963409/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.98392963409"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.98795938492/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.98795938492"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
SetActivePlots((0, 1))
DeleteActivePlots()
OpenDatabase("output.old.99632906914/output.visit", 0)
AddPlot("Pseudocolor", "Eta001", 1, 1)
AddPlot("Mesh", "Mesh", 1, 1)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "output.old.99632906914"
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
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
