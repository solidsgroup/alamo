import alamo
from mayavi import mlab
import numpy as np

alamo.Util.Initialize()

model = alamo.Model.Interface.GB.SH()

model.Define(0.0,1.0,1.0)

def plot(func,offx=0,offy=0):
    r = 1.0
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:201j]

    y = r * sin(theta) * cos(phi)
    z = r * sin(theta) * sin(phi)
    x = r * cos(theta)


    def Wwrapper(a,b,c): return func([a,b,c])
    W = np.vectorize(Wwrapper)
    s = W(x,y,z)
    print(np.amax(s),np.amin(s))
    mlab.mesh(x-2*offx, y-2*offy, z, scalars=s, colormap='jet')#,vmin=-1.0,vmax=1.0)


alpha = 0.000000001
dn0 = np.array([alpha,0,0])
dn1 = np.array([0,alpha,0])
dn2 = np.array([0,0,alpha])

def numericDerivative0(n): return (model.W(n+dn0) - model.W(n-dn0))/2.0/alpha
def numericDerivative1(n): return (model.W(n+dn1) - model.W(n-dn1))/2.0/alpha
def numericDerivative2(n): return (model.W(n+dn2) - model.W(n-dn2))/2.0/alpha

def exactDerivative0(n): return model.DW(n,np.array([1,0,0]))
def exactDerivative1(n): return model.DW(n,np.array([0,1,0]))
def exactDerivative2(n): return model.DW(n,np.array([0,0,1]))

def diff0(n): return abs(numericDerivative0(n) - exactDerivative0(n))
def diff1(n): return abs(numericDerivative1(n) - exactDerivative1(n))
def diff2(n): return abs(numericDerivative2(n) - exactDerivative2(n))


mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.clf()

plot(model.W)
plot(numericDerivative0,0,1)
plot(numericDerivative1,1,1)
plot(numericDerivative2,2,1)
plot(exactDerivative0,0,2)
plot(exactDerivative1,1,2)
plot(exactDerivative2,2,2)
plot(diff0,0,3)
plot(diff1,1,3)
plot(diff2,2,3)

mlab.show()

alamo.Util.Finalize()
