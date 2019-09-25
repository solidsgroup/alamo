import alamo
import numpy as np
import random
import pylab

sin = np.sin
cos = np.cos

alamo.Util.Initialize()

model = alamo.Model.Interface.GB.SH()
model.Define(0.0,1.0,1.0)


pi = np.pi

def getRandomUnitVector():
    phi = random.random() * 2.0*pi
    theta = random.random() * pi
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])

N = getRandomUnitVector();
B = getRandomUnitVector();
B = B - N*(B.dot(N));
B = B / np.sqrt(B.dot(B));

alpha = np.linspace(0,2*pi,10000)
dalpha = alpha[1]-alpha[0]

W = []
DW = []
DDW = []
for i in range(0,len(alpha)):
    n =  N*cos(alpha[i]) + B*sin(alpha[i])
    b = -N*sin(alpha[i]) + B*cos(alpha[i])
    W.append(model.W(n))
    DW.append(model.DW(n,b))
    DDW.append(model.DDW(n,b))

DWn = []
for i in range(0,len(alpha)):
    DWn.append((W[i]-W[i-1])/dalpha)

DDWn = []
for i in range(0,len(alpha)):
    DDWn.append((DW[i]-DW[i-1])/dalpha)


pylab.plot(alpha,W,label="W")
pylab.plot(alpha,DW,label="DW - implemented",linewidth=4)
pylab.plot(alpha,DWn,label="DW - numerical")
pylab.plot(alpha,DDW,label="DDW - implemented",linewidth=4)
pylab.plot(alpha,DDWn,label="DDW - numerical")

#DW = [np.array(model.DW(n*cos(a) + b*sin(a))) for a in alpha]
#pylab.plot(alpha,DW)



pylab.legend()
pylab.show()


alamo.Util.Finalize()
