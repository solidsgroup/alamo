import alamo

alamo.Util.Initialize()

mytest = alamo.Test.Operator.Elastic()

mytest.setBounds([1.0,1.0,1.0])

for lev in [1,2,3]:
    print("Levels of refinement: " + str(lev))
    mytest.Define([32,32,32],lev,1,mytest.Grid.XYZ)
    failed = mytest.TrigTest(0,0,1,"")
    if failed: print("Test failed")
    else: print ("Test passed")

alamo.Util.Finalize()
