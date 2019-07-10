#import alamo
#import numpy as np
#import matplotlib.pyplot as plt

#alamo.Util.Initialize()

#model = alamo.Model.Solid.CrystalPlastic.CrystalPlastic()
#es = np.ndarray(shape=(3,3))
#final = 2/1e-7
#dt = 1e-5
#model.Setdt(dt)
#a = 1e-5
#stressx = []
#strainx = []
#for i in range(int(110)):
#	es[0][0] = a * i
#	model.SetEs(es)
#	model.UpdateSigma()
#	sigma = model.GetSigma()
#	for j in range(1000000): #int(.001/dt*i)
#		#model.UpdateSigma()
#		model.AdvanceEsp()
#		#model.UpdateSigma()
#	if i%10 == 0:
#		print(sigma)
#		strainx.append(es[0][0])
#		stressx.append(sigma[0][0])
#	#model.reset()
#esp = model.GetEsp()
#print(esp)
#plt.plot(strainx, stressx)
#plt.xlabel('Strain')
#plt.ylabel('Stress')
#plt.title('Crystal plastic stress-strain curve')
#plt.show()

#alamo.Util.Finalize()

import alamo
import numpy as np
import matplotlib.pyplot as plt

alamo.Util.Initialize()

model = alamo.Model.Solid.CrystalPlastic.CrystalPlastic()
model.Randomize()
es = np.ndarray(shape=(3,3))
sigma = np.ndarray(shape=(3,3))
T = 3.0
dt = 1e-6
model.Setdt(dt)
counter = 0
t = 0

stressx = []
strainx = []

while(t<T):
	es[0][0] = t
	sigma = model.UpdateSigma(es)
	model.Update(es,sigma,dt)
	t += dt
	if(counter % 100 ==0):
		strainx.append(es[0][0])
		stressx.append(sigma[0][0])
	counter= counter + 1

plt.plot(strainx, stressx)
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('Crystal plastic stress-strain curve')
plt.show()
alamo.Util.Finalize()
