import alamo
import numpy

alamo.include("Numeric/Stencil.H")
dx = numpy.array([1.0, 1.0], dtype=numpy.float64)
ret = alamo.Numeric.Gradient_Diagonal[alamo.Set.Matrix3](1,dx)


for i in range(0,2):
    for j in range(0,2):
        for k in range(0,2):
            print(ret(i,j,k))


