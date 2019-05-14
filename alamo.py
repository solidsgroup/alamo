import alamo
import sys

alamo.Initialize()
mytest = alamo.Elastic()

#mytest.Define(16,3)
#print(alamo.Grid.values)
mytest.Define([32,32,32],3,1,mytest.Grid.XYZ)

#mytest.setBounds(1.0,1.0,1.0)
mytest.setBounds([1.0,1.0,1.0])


ret = mytest.TrigTest(3,0,1,"testoutput")
print("Return value = " + str(ret))
