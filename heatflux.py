import alamo
import pylab
import numpy

alamo.include("Model/Chamber/Ballistic.H")
alamo.include("Unit/Unit.H")
alamo.include("Model/Propellant/Homogenize.H")

alamo.Unit.setLengthUnit("m")
alamo.Unit.setTimeUnit("s")

alamo.Util.Initialize()



model = alamo.Model.Propellant.Homogenize()


model.dispersion1       = alamo.Unit.Parse("0.1").normalized_value()
model.dispersion2       = alamo.Unit.Parse("1.0").normalized_value()
model.dispersion3       = alamo.Unit.Parse("1.0").normalized_value()
model.h1                = alamo.Unit.Parse("10.67_W/m^2").normalized_value()
model.h2                = alamo.Unit.Parse("9.94e6_W/m^2").normalized_value()
model.P_reference       = alamo.Unit.Parse("2_MPa").normalized_value()
model.rho_prop          = alamo.Unit.Parse("1344.0_kg/m^3").normalized_value() # density
model.k_prop            = alamo.Unit.Parse("0.1366_W/m/K").normalized_value() # thermal conductivity
model.cp_prop           = alamo.Unit.Parse("1476.0_J/kg/K").normalized_value() # specific heat
model.m_prop            = alamo.Unit.Parse("0.8e5_cm/s").normalized_value() # pre-exponential factor for Arrhenius Law
model.E_prop            = alamo.Unit.Parse("16000.0_K").normalized_value() # activation energy for Arrhenius Law
model.mob_prop          = True # bool - whether to include pressure in arrhenius law
model.pressure_exponent = alamo.Unit.Parse("0.5").normalized_value()

model.k4 = 1.0

mdot = numpy.linspace(0,10)
qdot = [model.get_qdot(m,5.0) for m in mdot]

pylab.plot(mdot,qdot)
pylab.show()


# model = alamo.Model.Chamber.Ballistic()
# 
# 
# 
# model.At = alamo.Unit.Parse("0.00007_m^2").normalized_value()
# model.R = alamo.Unit.Parse("287_J/kg/K").normalized_value()
# model.gamma = 1.25
# model.T0 = alamo.Unit.Parse("1500_K").normalized_value()
# model.pressure = alamo.Unit.Parse("2.0_MPa").normalized_value()
# 
# 
# dt = 0.001
# mdot = alamo.Unit.Parse("0.26_kg/s").normalized_value()
# 
# t = [0.0]
# p = [model.pressure]
# volume = [alamo.Unit.Parse("0.0254088_m^3").normalized_value()]
# 
# 
# while (t[-1] < 5.0):
#     t.append(t[-1]+dt)
#     # Advance(Set::Scalar dt, Set::Scalar mdot_in, Set::Scalar volume, Set::Scalar pressure)
#     # return {new_pressure, current_dpdt};
#     p.append(model.Advance(dt, mdot, volume[-1], p[-1])[0])
#     print(p[-1])
# 
# 
# pylab.plot(t,p)
# pylab.show()
