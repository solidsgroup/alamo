import alamo
import sys

print(sys.argv)
alamo.Initialize()
mytest = alamo.Elastic()

mytest.Define(16,2)

mytest.setBounds(1.0,1.0,1.0)

mytest.TrigTest(3,0,1,"testoutput")
