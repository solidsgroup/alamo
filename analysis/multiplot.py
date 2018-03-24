import os
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pylab import imshow,show,plot,pcolor
from wand.image import Image
filename = "examples/perturbed_bar/PerturbedBar"
Image_1 = mpimg.imread(filename + "40/Frequency_plot.png")
Image_2 = mpimg.imread(filename + "41/Frequency_plot.png")
Image_3 = mpimg.imread(filename + "42/Frequency_plot.png")
Image_4 = mpimg.imread(filename + "43/Frequency_plot.png")


plt.figure(figsize=(9,6.5))
plt.subplot(2,2,1)
plt.axis('off')
plt.imshow(Image_1)
plt.title(r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0008\ l_{gb}= 0.08$')

plt.subplot(2,2,2)
plt.imshow(Image_2)
plt.axis('off')
plt.title(r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.001\ l_{gb}= 0.08$')

plt.subplot(2,2,3)
plt.imshow(Image_3)
plt.axis('off')
plt.title(r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0012\ l_{gb}= 0.08$')

plt.subplot(2,2,4)
plt.imshow(Image_4)
plt.axis('off')
plt.title(r'$\sigma_1 = 2.05\ \sigma_0 = 0.205\ \beta = 0.0014\ l_{gb}= 0.08$')


plt.savefig("Multiplot_histogram_40-43__betas_plot.pdf")

plt.show()
