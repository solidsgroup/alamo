#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:29:09 2022

@author: mmeierdo
"""
import numpy as np 
from matplotlib import pyplot as plt
import math 
from scipy import integrate
import getpass 

user = getpass.getuser()

def exact(x,t,Ti = 300,P = 4.0, k = 1*0.4186e0, cp = 1297.9, rho = 1950, div=1):
    alpha = k / cp / rho
   
    q = 1e6  
    q = k * q
    T = []
    H = 0.002#x[-1]
    for xx in x:
        F0 = alpha * t / H**2
        qsi = xx / H
        th = Ti + (q * H / k) * (np.sqrt(4*F0/ np.pi) * np.exp(- ( qsi**2 / 4 / F0) ) - qsi*math.erfc(qsi / np.sqrt(4*F0)))
        T.append(th)
      
    return T, x

def exact_(x,t,Ti = 300,P = 4.0, k = 1*0.4186e0, cp = 1297.9, rho = 1950, div=1):
    alpha = k / cp / rho
   
    q = 1.2e6  
    q = k * q
    T = []
    H = x[-1] - x[0]
    for xx in x:
        F0 = alpha * t / H**2
        qsi = (xx - x[0] )/ H
        th = Ti + (q * H / k) * (np.sqrt(4*F0/ np.pi) * np.exp(- ( qsi**2 / 4 / F0) ) - qsi*math.erfc(qsi / np.sqrt(4*F0)))
        T.append(th)
      
    return T, x

def alamo(name, time, ini=0, Pressure = 4.0):
    r = np.loadtxt(name, skiprows=0)
    x1 = r[ini:, 0] 
    rr = r[:,1]    
    T, x = exact(x=x1, t=time, P=Pressure)   
    T2, x2 = exact(x=x1, t = time, P = Pressure)   
    return T, rr, x
    

def check(T, r, x, x1):
    y1 = [i**2 for i in r[ini:,1]]
    y2 = [i*j for i,j in zip(T,r[ini:,1])]
    y3 = [i**2 for i in T]

    exact_top = integrate.cumtrapz(y3, x)
    alamo_top = integrate.cumtrapz(y1, x1)
    und = integrate.cumtrapz(y2, x1)

    check = alamo_top[-1] / und[-1]
    check1 = exact_top[-1] / und[-1]
    print("check: " + str(check))
    print("check1: " + str(check1))

    check1 = 0.98
    
    return check, check1

k = 0.4186 # W / m / K
Pressure = 4.
alp=1 

ini = 0
time = [0.001, .1, .2, .3, .4, .5]
time = [i - 0.00 for i in time]
ots = [ "output1", "output2","output3"]
file = ["/heat0000.curve",
         "/heat0001.curve",
         "/heat0002.curve",
         "/heat0003.curve",
         "/heat0004.curve",
         "/heat0005.curve",
         ]

files = [["./"+j+i for i in file] for j in ots]

cmap = plt.get_cmap('coolwarm')
c = ['orangered','darkorange','gold','yellowgreen','seagreen', "lightseagreen"] 
fig1, ((ax, bx), (cx, dx)) = plt.subplots(2, 2, figsize=(14,12))
col = ['#5572df', '#5572df', '#81a5fb', '#a0bfff', '#e9785e', '#b40426']

T,r,x = alamo("./output0.0002/heat0000.curve",0.5)
ax.plot(x, r, label="Alamo Solution", ls = 'none', marker='o',markeredgecolor='#b40426',markerfacecolor='none', markersize=5)
T,r,x = alamo("./output0.0001/heat0000.curve",0.5)
bx.plot(x, r, label="Alamo Solution", ls = 'none', marker='o',markeredgecolor='#b40426',markerfacecolor='none', markersize=5)
T,r,x = alamo("./output0.00002/heat0000.curve",0.5)
cx.plot(x, r, label="Alamo Solution", ls = 'none', marker='o',markeredgecolor='#b40426',markerfacecolor='none', markersize=5)
T,r,x = alamo("./output0.00001/heat0006.curve",0.5)
dx.plot(x, r, label="Alamo Solution", ls = 'none', marker='o',markeredgecolor='#b40426',markerfacecolor='none', markersize=5)

        
T, x = exact_(np.linspace(0.0020,0.004, 1000), 0.5)
ax.plot(x,T, label = "Exact Solution", c="black", linewidth=2.5)  
bx.plot(x,T, label = "Exact Solution", c="black", linewidth=2.5)  
cx.plot(x,T, label = "Exact Solution", c="black", linewidth=2.5)  
dx.plot(x,T, label = "Exact Solution", c="black", linewidth=2.5)  

ax.set_xlim(0.001, 0.003)
bx.set_xlim(0.001, 0.003)
cx.set_xlim(0.001, 0.003)
dx.set_xlim(0.001, 0.003)

ax.grid()
bx.grid()
cx.grid()
dx.grid()

ax.set_ylabel("Temperature [K]")
bx.set_ylabel("Temperature [K]")
cx.set_ylabel("Temperature [K]")
dx.set_ylabel("Temperature [K]")

ax.set_xlabel("Distance [m]")
bx.set_xlabel("Distance [m]")
cx.set_xlabel("Distance [m]")
dx.set_xlabel("Distance [m]")

ax.set_title("$\epsilon$=0.0002")
bx.set_title("$\epsilon$=0.0001")
cx.set_title("$\epsilon$=0.00002")
dx.set_title("$\epsilon$=0.00001")


ax.legend(loc="upper left")
bx.legend(loc="upper left")
cx.legend(loc="upper left")
dx.legend(loc="upper left")

ax2 = ax.twinx()
bx2 = bx.twinx()
cx2 = cx.twinx()
dx2 = dx.twinx()

phi = []
eps = 0.0002
cx = 0.002
xphi =np.linspace(0,0.004,10000)
phi = [0.5+0.5*np.tanh(((i-cx))/eps) for i in xphi]
ax2.plot(xphi,phi, ls = "-", label="$\eta$", c="grey")
phi = []; eps = 0.0001
phi = [0.5+0.5*np.tanh(((i-cx))/eps) for i in xphi]
bx2.plot(xphi,phi, ls = "-", label="$\eta$", c="grey")
phi = []; eps = 0.00002
phi = [0.5+0.5*np.tanh(((i-cx))/eps) for i in xphi]
cx2.plot(xphi,phi, ls = "-", label="$\eta$", c="grey")
phi = []; eps = 0.00001
phi = [0.5+0.5*np.tanh(((i-cx))/eps) for i in xphi]
dx2.plot(xphi,phi, ls = "-", label="$\eta$", c="grey")

ax2.legend(loc="upper right")
bx2.legend(loc="upper right")
cx2.legend(loc="upper right")
dx2.legend(loc="upper right")

ax2.set_ylabel("$\eta$ [-]")
bx2.set_ylabel("$\eta$ [-]")
cx2.set_ylabel("$\eta$ [-]")
dx2.set_ylabel("$\eta$ [-]")

fig1.tight_layout()
fig1.savefig(fname="heatequation.pdf", dpi=600, format="pdf")
