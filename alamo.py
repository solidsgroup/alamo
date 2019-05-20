import alamo
from mayavi import mlab
import numpy as np
from scipy.special import sph_harm

alamo.Util.Initialize()

model = alamo.Model.Interface.GB.SH()

model.Define(0.0,1.0,1.0)



# Create a sphere
r = 1.0
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:201j, 0:2 * pi:201j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.clf()

print(model.DW([1,2,3]))

def myfunc(a,b,c): return model.W([a,b,c])

W = np.vectorize(myfunc)

s = W(x,y,z)

mlab.mesh(x, y, z, scalars=s, colormap='jet')
mlab.show()


alamo.Util.Finalize()
