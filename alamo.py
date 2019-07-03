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
es = np.ndarray(shape=(3,3))
final = 2/1e-7
model.Setdt(1e-4)
dt = 1e-3
stressx = []
strainx = []
for i in range(int(9e3)):
	es[0][0] = dt * i
	model.SetEs(es)
	
	for j in range(int(1e4)):
		model.UpdateSigma()
		model.AdvanceEsp()
	if i%1 == 0:
		sigma = model.GetSigma()
		print(sigma)
		strainx.append(es[0][0])
		stressx.append(sigma[0][0])
esp = model.GetEsp()
print(esp)
plt.plot(strainx, stressx)
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('Crystal plastic stress-strain curve')
plt.show()

alamo.Util.Finalize()
