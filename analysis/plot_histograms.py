import os
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pylab import imshow,show,plot,pcolor

filename_1 = "examples/perturbed_bar/PerturbedBar41/"
filename_2 = "examples/perturbed_bar/PerturbedBar42/"
filename_3 = "examples/perturbed_bar/PerturbedBar43/"
filename_4 = "examples/perturbed_bar/PerturbedBar44/"
filename_5 = "examples/perturbed_bar/PerturbedBar45/"
filename_6 = "examples/perturbed_bar/PerturbedBar47/"

x_1 = np.loadtxt(filename_1 + "data_analysis/final_density_x.txt")
x_2 = np.loadtxt(filename_2 + "data_analysis/final_density_x.txt")
x_3 = np.loadtxt(filename_3 + "data_analysis/final_density_x.txt")
x_4 = np.loadtxt(filename_4 + "data_analysis/final_density_x.txt")
x_5 = np.loadtxt(filename_5 + "data_analysis/final_density_x.txt")
x_6 = np.loadtxt(filename_6 + "data_analysis/final_density_x.txt")


y_1 = np.loadtxt(filename_1 + "data_analysis/final_density.txt")
y_2 = np.loadtxt(filename_2 + "data_analysis/final_density.txt")
y_3 = np.loadtxt(filename_3 + "data_analysis/final_density.txt")
y_4 = np.loadtxt(filename_4 + "data_analysis/final_density.txt")
y_5 = np.loadtxt(filename_5 + "data_analysis/final_density.txt")
y_6 = np.loadtxt(filename_6 + "data_analysis/final_density.txt")

plt.figure(figsize = (14.5,6.5))
plt.figure(1)

ax = plt.subplot(111)
ax.plot(x_1,y_1, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.001\ l_{gb}= 0.08$' )
ax.plot(x_2,y_2, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0012\ l_{gb}= 0.08$')
ax.plot(x_3,y_3, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0014\ l_{gb}= 0.08$')
ax.plot(x_4,y_4, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0016\ l_{gb}= 0.08$')
ax.plot(x_5,y_5, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0020\ l_{gb}= 0.08$')
ax.plot(x_6,y_6, label = r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.001\ l_{gb}= 0.05$')
plt.xlim([-1.5,1.5])
plt.ylim([0,2])
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("Histograms_41_45-47.pdf")
