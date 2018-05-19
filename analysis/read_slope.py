import os
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pylab import imshow,show,plot,pcolor
from contours_frequency import frequency
from contours_frequency import analysis

vector_slope = []



for i in range(148,179):
    filename = "examples/perturbed_bar/PerturbedBar" + str(i) + "/data_analysis/av_slope.txt"
    slope= abs(float(open(filename).readlines()[0]))
    vector_slope.append(slope)

print(vector_slope)

np.savetxt("slope_148_178",vector_slope)
