import os
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pylab import imshow,show,plot,pcolor
from wand.image import Image
filename = "Images_paper/"
Image_1 = mpimg.imread(filename + "theta_45_low_f.png")
Image_2 = mpimg.imread(filename + "theta_45_high_f.png")
Image_3 = mpimg.imread(filename + "theta_30_low_f.png")
Image_4 = mpimg.imread(filename + "theta_30_high_f.png")




plt.figure(figsize=(15,6.5))
plt.subplot(2,2,1)
plt.axis('off')
plt.imshow(Image_1)
plt.title(r'(a) $\beta= 0.0008$')

plt.subplot(2,2,2)
plt.imshow(Image_2)
plt.axis('off')
plt.title(r'(b) $\beta = 0.0001$')

plt.subplot(2,2,3)
plt.imshow(Image_3)
plt.axis('off')
plt.title(r'(c) $\beta = 0.0008$')

plt.subplot(2,2,4)
plt.imshow(Image_4)
plt.axis('off')
plt.title(r'(d) $\beta = 0.0001$')


plt.savefig("Multiplot_verstility.pdf")

plt.show()
