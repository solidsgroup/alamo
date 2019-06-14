import alamo
import numpy as np
import matplotlib.pyplot as plt

# Initialize Alamo
alamo.Util.Initialize()

# Define the Linear Elastic Isotropic Model
# The C++ implementation is defined in 
#     src/Model/Solid/LinearElastic/Isotropic
# and the python bindings are defined in
#     py/Model/Solid/LinearElastic/LinearElastic.cpy
model = alamo.Model.Solid.LinearElastic.Isotropic()
lamb = 1.0 # Lame parameter
mu = 1.0   # Shear modulus
model.Define(lamb,mu)

# Define a strain tensor.
# NOTE: you MUST define the strain tensor this way
# (i.e. eps = np.ndarray...)
# or the library will break
eps = np.ndarray(shape=(3,3))
eps *= 0
eps[0][0] = 0
stressx = []
strainx = []
for i in range(10):
	eps[0][0] = eps[0][0] + 0.01
	#eps[1][2] = 0
	print(eps)
	strainx.append(eps[0][0])
	# Compute the strain energy as defined by the model
	energy = model.W(eps)

	# Compute the stress tensor as defined by the model
	sigma = model.DW(eps)
	print(sigma[0][0])
	stressx.append(sigma[0][0])

plt.plot(strainx, stressx)
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('Linear elastic isotropic stress-strain curve')
plt.show()
# Finalize Alamo
alamo.Util.Finalize()
