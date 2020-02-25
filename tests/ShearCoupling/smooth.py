#!/usr/bin/env python3
import numpy, scipy, pylab



#raw = numpy.loadtxt("energy_fcc_sigma5.dat")
raw = numpy.loadtxt("test2.dat")

theta = raw[:,0]
dtheta = theta[1]-theta[0]

w = raw[:,1]

sigma = 5.0
moll = dtheta*(1./sigma/numpy.sqrt(2.*numpy.pi)) * numpy.exp(-0.5*(theta/sigma)**2)

ws = numpy.convolve(w,moll,mode='same')

pylab.plot((theta),w)
#pylab.plot(theta,moll)
pylab.plot((theta),ws)


pylab.tight_layout()
pylab.show()

#f = open("energy_fcc_sigma5_smooth.dat","w")
#for _t, _w in zip(theta,ws):
#    f.write(str(numpy.radians(2*_t))+" "+str(_w-0.5)+"\n")
#f.close()
#
f = open("test2_smooth.dat","w")
for _t, _w in zip(theta,ws):
    f.write(str(numpy.radians(_t+90))+" "+str(2*_w-0.65)+"\n")
f.close()
