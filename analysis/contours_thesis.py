import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread
import scipy.fftpack
from pylab import imshow,show,plot
from skimage.color import rgb2gray
from skimage import measure
from skimage import data, img_as_float
from scipy.interpolate import interp1d
#from skimage.segmentation import active_contour
import matplotlib.pyplot as plt
from PIL import Image
import os.path


headerfile= "Images_paper/analysis_gray0000.png"
filename = headerfile
im = imread(filename)
img = rgb2gray(im)
contours = measure.find_contours(img, 0.5)
Fs=512.0
Ts=1.0/Fs;
t=np.arange(0,1,Ts)
y = interp1d(contours[0][:,1]/max(contours[0][:,1]), contours[0][:,0])
y=y(t)-np.mean(y(t))
y=-y
#y=y/max(y);
n=len(y)
#print(n)
y_slope =2.0*y/(208.0)
y= y*2.0/208.0
x=np.linspace(0,10,512)
dx = 10.0/512.0;
dy= np.gradient(y_slope,dx);

fourier = np.fft.rfft(y,512)
k=np.arange(n)
T = n/Fs
frq= k/T
plt.figure(1)
plt.plot(x,y)
plt.ylim([-2,2])
plt.axis('scaled')
plt.yticks([-0.2 ,0.2])
plt.savefig('Capture_image.png')
#print(n)
#print(mean_slope)
frq=frq[range(n/2)]
fourier=fourier[range(n/2)]
#frequency Time of maximum and standard deviation
fourier = abs(fourier)

plt.figure(2)
plt.plot(frq,fourier)
plt.xlim([0 ,32])
plt.xlabel(r'Frequency ($\omega$)')
plt.ylabel('Intensity')
plt.savefig('Analysis_frequency.png')
show()
