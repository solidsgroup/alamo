import alamo
import numpy as np

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
eps[0][0] = 0.1

# Compute the strain energy as defined by the model
energy = model.W(eps)
print(energy)

# Compute the stress tensor as defined by the model
sigma = model.DW(eps)
print(sigma)

# Finalize Alamo
alamo.Util.Finalize()
