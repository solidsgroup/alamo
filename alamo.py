import alamo
from mayavi import mlab
import numpy as np
from scipy.special import sph_harm

alamo.Util.Initialize()

model = alamo.Model.Interface.GB.SH()

model.Define(0.0,1.0,1.0)

def plotW(func):
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

    def Wwrapper(a,b,c): return func([a,b,c])
    W = np.vectorize(Wwrapper)
    s = W(x,y,z)
    mlab.mesh(x, y, z, scalars=s, colormap='jet')
    mlab.show()

def checkDerivative():
    r = 1.0
    pi = np.pi
    cos = np.cos
    sin = np.sin

    PHI = np.linspace(0,pi,51)
    THETA = np.linspace(0,2*pi,101)

    for phi in PHI:
        for theta in THETA:
            n = np.array([sin(theta) * cos(phi),
                          sin(theta) * sin(phi),
                          cos(theta)])
            
            alpha = 0.0000001

            dn0 = np.array([alpha,0,0])
            dn1 = np.array([0,alpha,0])
            dn2 = np.array([0,0,alpha])

            dw_numeric = np.array([(model.W(n+dn0) - model.W(n-dn0))/2.0/alpha,
                                   (model.W(n+dn1) - model.W(n-dn1))/2.0/alpha,
                                   (model.W(n+dn2) - model.W(n-dn2))/2.0/alpha])
            dw_exact = model.DW(n)
            print(dw_numeric[0],dw_exact[0])

            #dwnum = [(model.W(theta+alpha,phi) - model.W(theta-alpha,phi))/2./alpha,
            #         (model.W(theta,phi+alpha) - model.W(theta,phi-alpha))/2./alpha]
            #print(model.DW(theta,phi),dwnum)
            

            #print(n,model.W(n),model.W(phi,theta))


#plotW(model.W)

checkDerivative()


alamo.Util.Finalize()
