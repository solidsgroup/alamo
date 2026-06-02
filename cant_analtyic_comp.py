#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Path to your file
filename = "disp_y_curve.curve"

# Read data, skipping the header line
data = np.loadtxt(filename, skiprows=2)

def analytic_solution(x):
    a = 0.00000039
    l = 10.75
    return -a*x**2*(x**2-4*l*x+6*l**2)

x_data = data[:, 0]
y_data = data[:, 1]

y_analytic = analytic_solution(x_data)
error = (y_data - y_analytic)/(y_analytic)*100

l2_error = np.sqrt(np.mean(error**2))
max_error = np.max(np.abs(error))

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(8, 8),
    sharex=True,
    constrained_layout=True
)

# Solution comparison
ax1.plot(x_data, y_data, '--', label='Simulation')
ax1.plot(x_data, y_analytic, '-', label='Analytic')
ax1.set_ylabel('disp_y')
ax1.grid(True)
ax1.legend()

# Pointwise error
ax2.plot(x_data, error, '-.', color='red')
ax2.axhline(0.0, color='black', linestyle='--', linewidth=1)
ax2.set_xlabel('x')
ax2.set_ylabel('Error')
ax2.set_title(
    f'Pointwise Error\n'
    f'L2 = {l2_error:.3e}, Max = {max_error:.3e}'
)
ax2.grid(True)

plt.show()