import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import imshow,show,plot,pcolor
from contours_frequency import frequency
from contours_frequency import analysis

filename = "/home/josep/wheeler/output7.visit";
output_base = "./out7"
i=50
header_images= []
time_vector = []
for f in open(filename).readlines():
    if (i == 50) :
    	headerfile=os.path.dirname(filename)+"/"+f.replace("\n","");
        timestep=float(open(headerfile).readlines()[7])
        timestep_formatted = '{num:4.02f}'.format(num=timestep).zfill(7)
    	outputfile=os.path.basename(os.path.dirname(f.replace("\n","")));
        print(timestep_formatted)
        i= 1     
        header_images.append(timestep_formatted)
        #print(header_images)
        time_vector.append(timestep)
    else :
	i=i+1

np.savetxt("header_images.txt",header_images, fmt="%s")



filename = "/home/josep/overleaf/doc_prov/results/header_images.txt"
frequency_vector= []
fourier_matrix=[]
for f in open(filename).readlines():
    headerfile=os.path.dirname(filename)+"/out7/"+f.replace("\n","")+".png";
    fourier, frq = frequency(headerfile)
    max_frequency = analysis(fourier,frq)
    print(max_frequency)
    frequency_vector.append(max_frequency)
    fourier_matrix.append(fourier)
pcolor(np.transpose(np.array(fourier_matrix)))
#show()
#np.savetxt("freq_matrix",fourier_matrix)
plot(time_vector,frequency_vector,color='white')
show()


