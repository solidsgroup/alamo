import alamo
import pylab

alamo.include("Model/Chamber/Ballistic.H")

model = alamo.Model.Chamber.Ballistic()


model.At = 1.0
model.T0 = 500
model.R = 0.000287
model.gamma = 1.4
model.T0 = 300
model.pressure = 4.0


dt = 0.01
mdot = 1.0

t = [0.0]
p = [model.pressure]
volume = [0.01]

while (t[-1] < 10.0):
    t.append(t[-1]+dt)
    p.append(model.Advance(dt, mdot, volume[-1], p[-1]))
    print(p[-1])



pylab.plot(t,p)
pylab.show()
