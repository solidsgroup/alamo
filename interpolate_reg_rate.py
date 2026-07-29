#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import interp1d

# Original data
pressure = np.array([
    2.76, 2.85, 3.54, 4.14, 4.23, 5.52, 5.61, 6.99,
    8.27, 8.37, 10.34, 10.44, 11.72, 12.41, 13.10, 13.79
])

value = np.array([
    3.231940e-01, 3.350700e-01, 4.040823e-01, 4.610793e-01,
    4.695766e-01, 5.867199e-01, 5.950946e-01, 7.215921e-01,
    8.342342e-01, 8.426099e-01, 1.003589e+00, 1.010982e+00,
    1.105394e+00, 1.153453e+00, 1.200085e+00, 1.244354e+00
])

# Pressures to evaluate
new_pressure = np.arange(0, 7, 1)

# Linear interpolation (with extrapolation outside the data range)
interp_func = interp1d(
    pressure,
    value,
    kind='linear',
    fill_value='extrapolate'
)

new_value = interp_func(new_pressure)

# Print results
print("Pressure    Interpolated Value")
for p, v in zip(new_pressure, new_value):
    print(f"{p:8.2f}    {v:.6e}")
