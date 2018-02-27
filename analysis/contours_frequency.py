"""
===============
Contour finding
===============

We use a marching squares method to find constant valued contours in an image.
In ``skimage.measure.find_contours``, array values are linearly interpolated
to provide better precision of the output contours. Contours which intersect
the image edge are open; all others are closed.

The `marching squares algorithm
<http://www.essi.fr/~lingrand/MarchingCubes/algo.html>`__ is a special case of
the marching cubes algorithm (Lorensen, William and Harvey E. Cline. Marching
Cubes: A High Resolution 3D Surface Construction Algorithm. Computer Graphics
(SIGGRAPH 87 Proceedings) 21(4) July 1987, p. 163-170).

"""
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

def frequency(filename):
        im = imread(filename)
        img = rgb2gray(im)
        contours = measure.find_contours(img, 0.5)
        Fs=128.0
        Ts=1.0/Fs;
        t=np.arange(0,1,Ts)
        y = interp1d(contours[0][:,1]/max(contours[0][:,1]), contours[0][:,0])
        y=y(t)-np.mean(y(t))
        y=y/max(y);
        n=len(y)
        fourier = np.fft.rfft(y,128)
        k=np.arange(n)
        T = n/Fs
        frq= k/T
        #print(frq)
        #print(k)
        frq=frq[range(n/2)]
        fourier=fourier[range(n/2)]
        #frequency Time of maximum and standard deviation
        fourier = abs(fourier)

        #array = np.array(frq , fourier)
        
        return fourier,frq


def analysis(fourier,frq):
	max_freq_i = np.argmax(fourier)
	max_frequency = frq[max_freq_i]

        return max_frequency


#fourier,frq = frequency("../doc_prov/results/out/0008.00.png")
#standar_deviation = np.std(fourier)
#max_freq_i = np.argmax(fourier)
#max_frequency = frq[max_freq_i]
#print(standar_deviation)
#print(max_frequency)
#print(max_freq_i)
#np.save("fourier", fourier)
#np.savetxt("fourier",fourier)

