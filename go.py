import alamo
from mayavi import mlab
import numpy as np


alamo.Util.Initialize()
q4 = alamo.FiniteElement.Q4()
q4.Define(1.0,1.0)


phi0 = np.vectorize(lambda x,y: q4.Phi(0,[x,y]))
phi1 = np.vectorize(lambda x,y: q4.Phi(1,[x,y]))
phi2 = np.vectorize(lambda x,y: q4.Phi(2,[x,y]))
phi3 = np.vectorize(lambda x,y: q4.Phi(3,[x,y]))


X, Y = np.meshgrid(np.linspace(0,1,100),np.linspace(0,1,100))

P0 = phi0(X,Y)
P1 = phi1(X,Y)
P2 = phi2(X,Y)
P3 = phi3(X,Y)






mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.mesh(X,       Y,       P0, scalars=P0, colormap='jet'); mlab.axes()
mlab.mesh(X,       Y + 2.0, P1, scalars=P1, colormap='jet'); mlab.axes()
mlab.mesh(X + 2.0, Y ,      P2, scalars=P2, colormap='jet'); mlab.axes()
mlab.mesh(X + 2.0, Y + 2.0, P3, scalars=P3, colormap='jet'); mlab.axes()

mlab.show()

alamo.Util.Finalize()
