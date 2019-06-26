import alamo
import numpy as np
import matplotlib.pyplot as plt

alamo.Util.Initialize()

model = alamo.Model.Solid.CrystalPlastic.CrystalPlastic()
es = np.ndarray(shape=(3,3))
final = 2/1e-7
model.Setdt(1e-5)
dt = 1e-4
stressx = []
strainx = []
for i in range(110):
	es[0][0] = dt * i
	model.SetEs(es)
	model.UpdateSigma()
	for j in range(int(1/1e-5*i)):
		model.AdvanceEsp()
	if i%10 == 0:
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
