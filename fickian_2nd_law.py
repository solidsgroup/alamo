import numpy as np
import matplotlib.pyplot as plt
from math import erf

# Parameters
L = 1.0
D = 2.074e-5
D = 0.2074
C = [1,0]
times = [1, 10, 50, 200, 500]
times = [.01, .05, .15, .50]
N_terms = 300  # number of Fourier terms in series

# Spatial grid
x = np.linspace(0, L, 30)

def f(x):
    # Initial condition
    s = np.zeros(len(x))
    s[x<=L/2] = 1.0
    eps = 1e-2
    for i in range(len(x)):
        s[i] = 0.5 * (1.0 - erf((x[i] - L/2.0)/eps))
    return s

def h(x):
    # Steady state solution
    return C[0] + (C[1] - C[0]) * x / L

def v0(x):
    return f(x) - h(x)

def C_xt(x, t, L, D, N_terms=200):
    # steady-state part
    C_val = h(x)
    # transient Fourier series
    for m in range(1, N_terms + 1):
        B = 2.0/L * np.trapz(np.sin(m*np.pi*x/L)*v0(x), x)
        term =  B * np.sin(m*np.pi*x/L) * np.exp(-D*(m*np.pi/L)**2.0*t)
        C_val += term
    return C_val

# Plot
plt.figure(figsize=(8,6))
plt.plot(x,f(x), 'k.', markersize=10)
dat = np.loadtxt(f'fickian{0:.2f}.curve', skiprows=1)
plt.plot(dat[:,0], dat[:,1], label=f"t={0:.2f}")
for t in times:
    C_vals = C_xt(x, t, L, D, N_terms)
    plt.plot(x, C_vals, 'k.', markersize=10)
    dat = np.loadtxt(f'fickian{t:.2f}.curve', skiprows=1)
    plt.plot(dat[:,0], dat[:,1], label=f"t={t:.2f}")
plt.plot([],[],'k.', markersize=10, label='Analytic')
plt.plot(x,h(x), '--', color='tab:grey', label=r"t=$\infty$")

plt.xlabel("x")
plt.ylabel("C(x,t)")
plt.title(r"$\dfrac{\partial C}{\partial t} = D_{AB} \dfrac{\partial^2 C}{\partial x^2}$" + "\n" + r"$D_{AB} = $" + f"{D}")
plt.legend()
plt.gca().set_xticks(np.linspace(0,1,11))
plt.gca().set_yticks(np.linspace(0,1,11))
plt.grid(True)
plt.show()

