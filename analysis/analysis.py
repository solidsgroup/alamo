import os
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pylab import imshow,show,plot,pcolor
from contours_frequency import frequency
from contours_frequency import analysis

filename = "examples/perturbed_bar/PerturbedBar163/"
output_base =filename + "/Images"
i=1
header_images= []
time_vector = []
read_file = filename+"output/output.visit"
for f in open(read_file).readlines():
    if (i == 1) :
        headerfile=os.path.dirname(read_file)+"/"+f.replace("\n","");
        timestep=float(open(headerfile).readlines()[7])
        timestep_formatted = '{num:4.02f}'.format(num=timestep).zfill(7)
        outputfile=os.path.basename(os.path.dirname(f.replace("\n","")));
        print(timestep_formatted)
        i = 1    
        header_images.append(timestep_formatted)
        #print(header_images)
        time_vector.append(timestep)
    else:
        i=i+1

header_save = filename+"Images/header_images.txt"
np.savetxt(header_save ,header_images, fmt="%s")
    
    
    
header_txt = filename + "Images/header_images.txt"
max_frequency_vector= []
fourier_matrix=[]
density_matrix = []
av_slope = []
max_freq = []
i=0
for f in open(header_txt).readlines():
    headerfile=os.path.dirname(header_txt)+"/"+f.replace("\n","")+".png";
    fourier, frq,dy = frequency(headerfile)
    max_frequency = analysis(fourier,frq)
    print(max_frequency)
    max_frequency_vector.append(max_frequency)
    fourier_matrix.append(fourier)
    density = gaussian_kde(dy)
    xs= np.linspace(-10,10,1000)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    slope_density = density(xs)
    density_matrix.append(slope_density)
    plt.figure(1)
    num_lines = float(sum(1 for line in open(header_txt)))
    grey_intensity = str(0.9-0.9*float((i/num_lines)))
    i=i+1
        
slope_density = slope_density[:]
max_slope_density = max(slope_density)
max_slope = xs[slope_density.argmax()]
av_slope.append(max_slope)
print(av_slope)
max_frequency =  analysis(fourier,frq)
print(max_frequency)
max_freq.append(max_frequency)
np.savetxt(filename + "data_analysis/final_density.txt",density(xs))
np.savetxt(filename + "data_analysis/final_density_x.txt",xs)
np.savetxt(filename + "data_analysis/av_slope.txt",av_slope)
np.savetxt(filename + "data_analysis/final_freq.txt",max_freq)

plt.plot(xs,density(xs), color=grey_intensity )
plt.xlim([-1.2,1.2])
plt.ylim([0,2])
plt.xlabel('Slope')
plt.ylabel('Probability density')
        
np.savetxt(filename+ "data_analysis/max_freq_vector.txt",max_frequency_vector)
plt.savefig(filename + "slope_analysis.png")
plt.savefig(filename + "slope_analysis.pdf")


plt.figure(2)
pcolor(np.transpose(np.array(fourier_matrix)))
plt.ylim([0,40])
plt.xlabel('t')
h =plt.ylabel(r'$\omega$')
h.set_rotation(0)
plt.savefig(filename + "Frequency_plot.pdf")
plt.savefig(filename + "Frequency_plot.png")

np.savetxt(filename + "data_analysis/freq_matrix",fourier_matrix)
np.savetxt(filename + "data_analysis/density_dy",density_matrix)
    #plot(time_vector,frequency_vector,color='white')
    

#show()
    
