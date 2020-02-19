#!/usr/bin/env python3
import numpy, scipy, pylab



raw = numpy.loadtxt("energy_fcc_sigma5.dat")

theta = raw[:,0]

w = raw[:,1]

sigma = 3.0
moll = 0.25*(1./sigma/numpy.sqrt(2.*numpy.pi)) * numpy.exp(-0.5*(theta/sigma)**2)

ws = numpy.convolve(w,moll,mode='same')

pylab.plot(numpy.radians(theta),w)
#pylab.plot(theta,moll)
pylab.plot(numpy.radians(theta),ws)


pylab.tight_layout()
pylab.show()

f = open("energy_fcc_sigma5_smooth.dat","w")
for _t, _w in zip(theta,ws):
    f.write(str(numpy.radians(2*_t))+" "+str(_w-0.5)+"\n")
f.close()
