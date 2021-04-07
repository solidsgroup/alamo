import sys
import os

print(sys.argv)
filename = sys.argv[-2]
imgname = sys.argv[-1]
print("filename = " + filename)
print("imgname  = " +imgname)

if not os.path.isfile(filename+"/nodeoutput.visit"): exit()

OpenDatabase("localhost:{}/nodeoutput.visit".format(filename), 0)
AddPlot("Pseudocolor", "crack", 1, 1)
DefineScalarExpression("operators/ConnectedComponents/Mesh", "cell_constant(<Mesh>, 0.)")
DefineCurveExpression("operators/DataBinning/1D/Mesh", "cell_constant(<Mesh>, 0)")
DefineScalarExpression("operators/DataBinning/2D/Mesh", "cell_constant(<Mesh>, 0)")
DefineScalarExpression("operators/DataBinning/3D/Mesh", "cell_constant(<Mesh>, 0)")
DefineScalarExpression("operators/Flux/Mesh", "cell_constant(<Mesh>, 0.)")
DefineCurveExpression("operators/Lineout/crack", "cell_constant(<crack>, 0.)")
DefineCurveExpression("operators/Lineout/disp001", "cell_constant(<disp001>, 0.)")
DefineCurveExpression("operators/Lineout/disp002", "cell_constant(<disp002>, 0.)")
DefineCurveExpression("operators/Lineout/materials001", "cell_constant(<materials001>, 0.)")
DefineCurveExpression("operators/Lineout/materials002", "cell_constant(<materials002>, 0.)")
DefineCurveExpression("operators/Lineout/stress001", "cell_constant(<stress001>, 0.)")
DefineCurveExpression("operators/Lineout/stress002", "cell_constant(<stress002>, 0.)")
DefineCurveExpression("operators/Lineout/stress003", "cell_constant(<stress003>, 0.)")
DefineCurveExpression("operators/Lineout/stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/ModelFit/distance", "point_constant(<Mesh>, 0)")
DefineScalarExpression("operators/ModelFit/model", "point_constant(<Mesh>, 0)")
DefineScalarExpression("operators/StatisticalTrends/Mean/crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Mean/stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Residuals/stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Slope/stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Std. Dev./stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Sum/stress004", "cell_constant(<stress004>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/crack", "cell_constant(<crack>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/disp001", "cell_constant(<disp001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/disp002", "cell_constant(<disp002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/materials001", "cell_constant(<materials001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/materials002", "cell_constant(<materials002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/stress001", "cell_constant(<stress001>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/stress002", "cell_constant(<stress002>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/stress003", "cell_constant(<stress003>, 0.)")
DefineScalarExpression("operators/StatisticalTrends/Variance/stress004", "cell_constant(<stress004>, 0.)")
DefineVectorExpression("operators/SurfaceNormal/Mesh", "cell_constant(<Mesh>, 0.)")
DefineScalarExpression("bndry", "materials001*materials002*crack")
AddPlot("Pseudocolor", "bndry", 1, 1)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.useBelowMinColor = 0
PseudocolorAtts.belowMinColor = (0, 0, 0, 255)
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.useAboveMaxColor = 0
PseudocolorAtts.aboveMaxColor = (0, 0, 0, 255)
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "xray"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.Ramp  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)
DrawPlots()
TimeSliderPreviousState()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = filename
SaveWindowAtts.fileName = filename+"/"+imgname
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate, LZW
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.pixelData = 1
SaveWindowAtts.advancedMultiWindowSave = 0
SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
SaveWindowAtts.subWindowAtts.win1.layer = 0
SaveWindowAtts.subWindowAtts.win1.transparency = 0
SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
SaveWindowAtts.subWindowAtts.win2.position = (0, 0)
SaveWindowAtts.subWindowAtts.win2.size = (128, 128)
SaveWindowAtts.subWindowAtts.win2.layer = 0
SaveWindowAtts.subWindowAtts.win2.transparency = 0
SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
SaveWindowAtts.subWindowAtts.win3.size = (128, 128)
SaveWindowAtts.subWindowAtts.win3.layer = 0
SaveWindowAtts.subWindowAtts.win3.transparency = 0
SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
SaveWindowAtts.subWindowAtts.win4.position = (0, 0)
SaveWindowAtts.subWindowAtts.win4.size = (128, 128)
SaveWindowAtts.subWindowAtts.win4.layer = 0
SaveWindowAtts.subWindowAtts.win4.transparency = 0
SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
SaveWindowAtts.subWindowAtts.win5.layer = 0
SaveWindowAtts.subWindowAtts.win5.transparency = 0
SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
SaveWindowAtts.subWindowAtts.win6.layer = 0
SaveWindowAtts.subWindowAtts.win6.transparency = 0
SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
SaveWindowAtts.subWindowAtts.win7.layer = 0
SaveWindowAtts.subWindowAtts.win7.transparency = 0
SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
SaveWindowAtts.subWindowAtts.win8.layer = 0
SaveWindowAtts.subWindowAtts.win8.transparency = 0
SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
SaveWindowAtts.subWindowAtts.win9.layer = 0
SaveWindowAtts.subWindowAtts.win9.transparency = 0
SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
SaveWindowAtts.subWindowAtts.win10.layer = 0
SaveWindowAtts.subWindowAtts.win10.transparency = 0
SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
SaveWindowAtts.subWindowAtts.win11.layer = 0
SaveWindowAtts.subWindowAtts.win11.transparency = 0
SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
SaveWindowAtts.subWindowAtts.win12.layer = 0
SaveWindowAtts.subWindowAtts.win12.transparency = 0
SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
SaveWindowAtts.subWindowAtts.win13.layer = 0
SaveWindowAtts.subWindowAtts.win13.transparency = 0
SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
SaveWindowAtts.subWindowAtts.win14.layer = 0
SaveWindowAtts.subWindowAtts.win14.transparency = 0
SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
SaveWindowAtts.subWindowAtts.win15.layer = 0
SaveWindowAtts.subWindowAtts.win15.transparency = 0
SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
SaveWindowAtts.subWindowAtts.win16.layer = 0
SaveWindowAtts.subWindowAtts.win16.transparency = 0
SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
SaveWindowAtts.opts.types = ()
SaveWindowAtts.opts.help = ""
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
exit()
