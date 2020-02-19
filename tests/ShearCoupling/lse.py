import numpy, pylab

mode = "abssin"

#sigma0 = 0.075
sigma0 = 0.01
sigma1 = 0.07
theta0 = 0

theta = numpy.linspace(-1.1*numpy.pi, 1.1*numpy.pi, 400)

if mode == "sin":
    energy = sigma0 + 0.5*sigma1*(1.0 - numpy.cos(4.0*theta - numpy.radians(theta0)))
if mode == "abssin":
    sin = numpy.sin(2*(theta-numpy.radians(theta0)))
    eps = 0.15
    energy = sigma0 + sigma1 * (eps*numpy.log(numpy.exp(sin/eps) + numpy.exp(-sin/eps)))


f = open("energy_abssin_v2.dat",'w')
for t,e in zip(theta,energy):
    f.write(str((t))+" "+str(e)+'\n')

pylab.plot(theta,energy,marker='o')
pylab.show()
