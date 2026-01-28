import alamo
import pylab

alamo.include("Model/Chamber/Ballistic.H")
alamo.include("Unit/Unit.H")

alamo.Unit.setLengthUnit("m")
alamo.Unit.setTimeUnit("s")

alamo.Util.Initialize()

model = alamo.Model.Chamber.Ballistic()



model.At = alamo.Unit.Parse("0.00007_m^2").normalized_value()
model.R = alamo.Unit.Parse("287_J/kg/K").normalized_value()
model.gamma = 1.25
model.T0 = alamo.Unit.Parse("1500_K").normalized_value()
model.pressure = alamo.Unit.Parse("2.0_MPa").normalized_value()


dt = 0.001
mdot = alamo.Unit.Parse("0.26_kg/s").normalized_value()

t = [0.0]
p = [model.pressure]
volume = [alamo.Unit.Parse("0.0254088_m^3").normalized_value()]


while (t[-1] < 5.0):
    t.append(t[-1]+dt)
    # Advance(Set::Scalar dt, Set::Scalar mdot_in, Set::Scalar volume, Set::Scalar pressure)
    # return {new_pressure, current_dpdt};
    p.append(model.Advance(dt, mdot, volume[-1], p[-1])[0])
    print(p[-1])


pylab.plot(t,p)
pylab.show()
