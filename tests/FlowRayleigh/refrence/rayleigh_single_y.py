"""
===============================================================================
Rayleigh Flow (Stokes' First Problem) Analysis Script - SLIDING LEFT PLATE
===============================================================================

This script extracts velocity profiles from AMReX Rayleigh flow simulations,
compares them with the analytical transient solution, and generates visualization plots.

CONFIGURATION: Left plate at x=x_min SLIDES with velocity U in the +y direction.
Y-velocity diffuses in the +x direction (away from the left plate).

Analytical Solution:
    v_y(x,t) = U * erfc((x - x_min) / sqrt(4*nu*t))
    
    where:
    - U: plate velocity in y-direction (m/s)
    - nu: kinematic viscosity (m^2/s) = mu/rho
    - x: distance from plate (m)
    - t: time (s)
    - erfc: complementary error function

Features:
- Extracts 1D velocity profile along x-direction sampling velocityy
- Compares numerical solution with analytical Rayleigh solution
- Creates multi-panel plot showing transient evolution (up to 6 timesteps)
- Customizable plotting parameters (fonts, titles, line weights)
- Outputs plots as both .eps and .png formats
- Optional GIF generation showing velocity profile evolution
- Error metrics (L2 and L-infinity norms) for each timestep
- Final timestep comparison plot
- Semilog error plot

Required inputs:
- U: Plate velocity in y-direction (m/s)
- mu: Dynamic viscosity (Pa-s)
- rho: Fluid density (kg/m^3)
- AMReX output directory path
===============================================================================
"""

import yt
import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image
import glob
import re
from scipy.special import erfc

# Suppress yt's verbose output
yt.funcs.mylog.setLevel(40)

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

# Physical parameters (matching your input file)
U = 1.0                          # Plate velocity in Y-direction (m/s)
mu = 2.5                         # Dynamic viscosity (Pa-s)
rho = 100.0                      # Fluid density (kg/m^3)
nu = mu / rho                    # Kinematic viscosity (m^2/s)

# File paths
amrex_output_dir = r'..\..\..\bin\tests\FlowRayleigh\FlowRayleighy'

# Plotting customization
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 14
FONT_SIZE_LEGEND = 11
FONT_SIZE_TICK = 11
LINE_WIDTH_EXACT = 2.5
LINE_WIDTH_NUMERICAL = 2.0
PLOT_TITLE = "Rayleigh Flow (Stokes' First Problem) - Sliding Left Plate"

# Multi-panel plot settings
MAX_PANELS = 6                   # Maximum number of timesteps to show
PANEL_ROWS = 2
PANEL_COLS = 3
PANEL_FIGSIZE = (15, 10)

# Output settings
output_folder = './Images'
create_gif = False                # Set to True to generate GIF animation
gif_duration = 300               # Duration per frame in milliseconds
gif_filename = 'rayleigh_evolution.gif'

# Sampling location (y-coordinate where we sample the x-profile)
# Set to 0.0 or domain center depending on your setup
y_sample_location = 0.0          # Sample at y=0

# ============================================================================
# ANALYTICAL RAYLEIGH SOLUTION
# ============================================================================

def analytical_rayleigh_velocity(x_dist, t, U, nu):
    """
    Analytical solution for Rayleigh flow (Stokes' first problem).
    
    A vertical plate at x=0 slides in the y-direction with velocity U.
    The y-velocity diffuses into the semi-infinite fluid domain (x>0).
    
    Formula: v_y(x,t) = U * erfc(x / sqrt(4*nu*t))
    
    Parameters:
    -----------
    x_dist : array-like
        Distance from the plate (m), must be >= 0
    t : float
        Time since impulsive start (s), must be > 0
    U : float
        Plate velocity in y-direction (m/s)
    nu : float
        Kinematic viscosity (m^2/s)
    
    Returns:
    --------
    v_y : array-like
        Y-velocity at each x position (m/s)
    """
    if t <= 0:
        # At t=0, velocity is U at plate and 0 everywhere else
        return np.where(x_dist == 0, U, 0.0)
    
    # Similarity variable
    eta = x_dist / np.sqrt(4.0 * nu * t)
    
    # Rayleigh solution
    return U * erfc(eta)

# ============================================================================
# HELPER FUNCTION TO EXTRACT TIMESTEP NUMBER
# ============================================================================

def extract_timestep_number(filename):
    """
    Extract the timestep number from a plot file name.
    Handles formats like: 00000cell, 00100cell, plt00000, etc.
    """
    match = re.search(r'(\d+)', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return 0

# ============================================================================
# SETUP
# ============================================================================

print("=" * 70)
print("RAYLEIGH FLOW (STOKES' FIRST PROBLEM) - SLIDING LEFT PLATE")
print("=" * 70)
print(f"\nPhysical Parameters:")
print(f"  Plate Velocity (U):       {U:.4f} m/s (Y-DIRECTION)")
print(f"  Density (rho):            {rho:.2e} kg/m^3")
print(f"  Dynamic Viscosity (mu):   {mu:.2e} Pa-s")
print(f"  Kinematic Viscosity (nu): {nu:.2e} m^2/s")
print(f"\nSampling Configuration:")
print(f"  Direction: X-axis (horizontal from left plate)")
print(f"  Velocity Component: velocityy (Y-velocity)")
print(f"  Y-location: {y_sample_location} m")

# Create output directory if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"\nCreated output directory: {output_folder}")

# ============================================================================
# FIND PLOT FILES
# ============================================================================

print("\n" + "=" * 70)
print("SCANNING FOR AMReX OUTPUT FILES")
print("=" * 70)

plot_files = []
for item in os.listdir(amrex_output_dir):
    item_path = os.path.join(amrex_output_dir, item)
    if os.path.isdir(item_path) and item.endswith('cell'):
        plot_files.append(item_path)

if not plot_files:
    print(f"ERROR: No *cell directories found in {amrex_output_dir}")
    print("Looking for directories like: 00000cell, 00100cell, etc.")
    exit(1)

# Sort plot files by timestep number
plot_files.sort(key=extract_timestep_number)

print(f"\nFound {len(plot_files)} plot files")
print(f"First: {os.path.basename(plot_files[0])} (timestep {extract_timestep_number(plot_files[0])})")
print(f"Last:  {os.path.basename(plot_files[-1])} (timestep {extract_timestep_number(plot_files[-1])})")

# ============================================================================
# SELECT SUBSET OF TIMESTEPS FOR MULTI-PANEL PLOT
# ============================================================================

num_timesteps = len(plot_files)
num_panels = min(MAX_PANELS, num_timesteps)

# Select evenly-spaced timesteps
if num_timesteps <= MAX_PANELS:
    selected_indices = list(range(num_timesteps))
else:
    selected_indices = np.linspace(0, num_timesteps - 1, num_panels, dtype=int)

selected_plot_files = [plot_files[i] for i in selected_indices]

print("\n" + "=" * 70)
print(f"SELECTED {num_panels} TIMESTEPS FOR ANALYSIS")
print("=" * 70)
for i, plot_file in enumerate(selected_plot_files):
    print(f"  Panel {i+1}: {os.path.basename(plot_file)} (timestep {extract_timestep_number(plot_file)})")

# ============================================================================
# EXTRACT DATA AND CREATE MULTI-PANEL PLOT
# ============================================================================

print("\n" + "=" * 70)
print("EXTRACTING DATA AND CREATING MULTI-PANEL PLOT")
print("=" * 70)

fig, axes = plt.subplots(PANEL_ROWS, PANEL_COLS, figsize=PANEL_FIGSIZE)
axes = axes.flatten()  # Flatten to 1D array for easy indexing

# Storage for error metrics
error_data = []

for panel_idx, plot_file in enumerate(selected_plot_files):
    print(f"\nProcessing panel {panel_idx + 1}/{num_panels}...")
    print(f"  File: {os.path.basename(plot_file)}")
    
    # Load dataset
    ds = yt.load(plot_file)
    simulation_time = float(ds.current_time)
    print(f"  Time: {simulation_time:.6e} s")
    
    # Get domain bounds
    x_min = float(ds.domain_left_edge[0])
    x_max = float(ds.domain_right_edge[0])
    y_min = float(ds.domain_left_edge[1])
    y_max = float(ds.domain_right_edge[1])
    
    # Determine y-sampling location (use specified or domain center)
    if y_sample_location is None:
        y_sample = (y_min + y_max) / 2.0
    else:
        y_sample = y_sample_location
    
    print(f"  Sampling at y = {y_sample:.6e} m")
    
    # Create a ray along x-direction at y=y_sample, z=0
    ray_start = ds.arr([x_min, y_sample, 0.0], 'code_length')
    ray_end = ds.arr([x_max, y_sample, 0.0], 'code_length')
    ray = ds.ray(ray_start, ray_end)
    
    # Sort by x coordinate
    sort_indices = np.argsort(ray['x'])
    x_numerical = np.array(ray['x'][sort_indices])
    
    # *** KEY CHANGE: Extract velocityy instead of velocityx ***
    velocity_numerical = np.array(ray['velocityy'][sort_indices])
    
    # Calculate distance from left plate
    x_dist_numerical = x_numerical - x_min
    
    print(f"  Extracted {len(x_numerical)} points")
    print(f"  X range: [{x_min:.6e}, {x_max:.6e}] m")
    print(f"  Velocity component: velocityy")
    
    # Compute analytical solution
    x_dist_analytical = np.linspace(0, x_max - x_min, 1000)
    velocity_analytical = analytical_rayleigh_velocity(x_dist_analytical, simulation_time, U, nu)
    
    # Convert back to absolute x coordinates for plotting
    x_analytical = x_dist_analytical + x_min
    
    # Plot on current panel
    ax = axes[panel_idx]
    ax.plot(x_analytical, velocity_analytical, 'b-', 
            linewidth=LINE_WIDTH_EXACT, label='Analytical', zorder=1)
    ax.plot(x_numerical, velocity_numerical, 'r--', 
            linewidth=LINE_WIDTH_NUMERICAL, label='Numerical', 
            zorder=2, alpha=0.8)
    
    # Formatting
    ax.set_xlabel('Position x (m)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Y-Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
    ax.set_title(f't = {simulation_time:.4e} s', 
                 fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([0, U * 1.1])
    
    # Compute error metrics
    velocity_analytical_interp = np.interp(x_numerical, x_analytical, velocity_analytical)
    velocity_error = velocity_numerical - velocity_analytical_interp
    l2_error = np.linalg.norm(velocity_error) / np.sqrt(len(velocity_error))
    linf_error = np.max(np.abs(velocity_error))
    
    error_data.append({
        'time': simulation_time,
        'l2_error': l2_error,
        'linf_error': linf_error
    })
    
    print(f"  L2 error: {l2_error:.6e}, L-inf error: {linf_error:.6e}")

# Hide unused panels if num_panels < MAX_PANELS
for i in range(num_panels, len(axes)):
    axes[i].axis('off')

# Overall title
fig.suptitle(PLOT_TITLE, fontsize=FONT_SIZE_TITLE + 2, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])

# Save multi-panel plot
multipanel_filename = os.path.join(output_folder, 'rayleigh_transient_multipanel')
plt.savefig(multipanel_filename + '.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(multipanel_filename + '.eps', format='eps', bbox_inches='tight')

print("\n" + "=" * 70)
print("MULTI-PANEL PLOT SAVED")
print("=" * 70)
print(f"  {multipanel_filename}.png")
print(f"  {multipanel_filename}.eps")

plt.show()

# ============================================================================
# CREATE FINAL TIMESTEP COMPARISON PLOT
# ============================================================================

print("\n" + "=" * 70)
print("CREATING FINAL TIMESTEP COMPARISON PLOT")
print("=" * 70)

last_plot = plot_files[-1]
ds_final = yt.load(last_plot)
simulation_time_final = float(ds_final.current_time)

x_min_final = float(ds_final.domain_left_edge[0])
x_max_final = float(ds_final.domain_right_edge[0])
y_min_final = float(ds_final.domain_left_edge[1])
y_max_final = float(ds_final.domain_right_edge[1])

# Determine y-sampling location
if y_sample_location is None:
    y_sample_final = (y_min_final + y_max_final) / 2.0
else:
    y_sample_final = y_sample_location

ray_start_final = ds_final.arr([x_min_final, y_sample_final, 0.0], 'code_length')
ray_end_final = ds_final.arr([x_max_final, y_sample_final, 0.0], 'code_length')
ray_final = ds_final.ray(ray_start_final, ray_end_final)

sort_indices_final = np.argsort(ray_final['x'])
x_numerical_final = np.array(ray_final['x'][sort_indices_final])

# *** KEY CHANGE: Extract velocityy instead of velocityx ***
velocity_numerical_final = np.array(ray_final['velocityy'][sort_indices_final])

x_dist_numerical_final = x_numerical_final - x_min_final

x_dist_analytical_final = np.linspace(0, x_max_final - x_min_final, 1000)
velocity_analytical_final = analytical_rayleigh_velocity(x_dist_analytical_final, simulation_time_final, U, nu)
x_analytical_final = x_dist_analytical_final + x_min_final

fig_final, ax_final = plt.subplots(figsize=(10, 8))

ax_final.plot(x_analytical_final, velocity_analytical_final, 'b-', 
              linewidth=LINE_WIDTH_EXACT, label='Analytical', zorder=1)
ax_final.plot(x_numerical_final, velocity_numerical_final, 'r--', 
              linewidth=LINE_WIDTH_NUMERICAL, label='Numerical', 
              zorder=2, alpha=0.8)

ax_final.set_xlabel('Position x (m)', fontsize=FONT_SIZE_LABEL)
ax_final.set_ylabel('Y-Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
ax_final.set_title(f"{PLOT_TITLE}\nFinal Time: t = {simulation_time_final:.6e} s", 
                   fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax_final.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
ax_final.grid(True, alpha=0.3)
ax_final.tick_params(labelsize=FONT_SIZE_TICK)
ax_final.set_xlim([x_min_final, x_max_final])
ax_final.set_ylim([0, U * 1.1])

plt.tight_layout()

final_filename = os.path.join(output_folder, 'rayleigh_final_comparison')
plt.savefig(final_filename + '.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(final_filename + '.eps', format='eps', bbox_inches='tight')

print(f"  {final_filename}.png")
print(f"  {final_filename}.eps")

plt.show()

# ============================================================================
# COMPUTE AND DISPLAY ERROR METRICS
# ============================================================================

print("\n" + "=" * 70)
print("ERROR METRICS SUMMARY")
print("=" * 70)
print(f"\n{'Time (s)':<15} {'L2 Error':<15} {'L-inf Error':<15}")
print("-" * 45)
for error_entry in error_data:
    print(f"{error_entry['time']:<15.6e} {error_entry['l2_error']:<15.6e} {error_entry['linf_error']:<15.6e}")

# Final timestep error
velocity_analytical_final_interp = np.interp(x_numerical_final, 
                                              x_analytical_final, 
                                              velocity_analytical_final)
velocity_error_final = velocity_numerical_final - velocity_analytical_final_interp
l2_error_final = np.linalg.norm(velocity_error_final) / np.sqrt(len(velocity_error_final))
linf_error_final = np.max(np.abs(velocity_error_final))

print("\n" + "=" * 70)
print("FINAL TIMESTEP ERROR METRICS")
print("=" * 70)
print(f"  Time:          {simulation_time_final:.6e} s")
print(f"  L2 error:      {l2_error_final:.6e}")
print(f"  L-inf error:   {linf_error_final:.6e}")

# ============================================================================
# CREATE SEMILOG ERROR PLOT FOR FINAL TIMESTEP
# ============================================================================

print("\n" + "=" * 70)
print("CREATING SEMILOG ERROR PLOT")
print("=" * 70)

abs_velocity_error_final = np.abs(velocity_error_final)
epsilon = 1e-16
abs_error_safe = abs_velocity_error_final + epsilon

fig_error, ax_error = plt.subplots(figsize=(10, 8))

ax_error.semilogy(x_numerical_final, abs_error_safe, 'k-', linewidth=LINE_WIDTH_NUMERICAL)

ax_error.set_xlabel('Position x (m)', fontsize=FONT_SIZE_LABEL)
ax_error.set_ylabel('Absolute Error |vy - v_true| (m/s)', fontsize=FONT_SIZE_LABEL)
ax_error.set_title(f'Rayleigh Flow: Absolute Error\nt = {simulation_time_final:.6e} s', 
                   fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax_error.grid(True, alpha=0.3, which='both')
ax_error.tick_params(labelsize=FONT_SIZE_TICK)

plt.tight_layout()

error_filename = os.path.join(output_folder, 'rayleigh_error_semilog')
plt.savefig(error_filename + '.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(error_filename + '.eps', format='eps', bbox_inches='tight')

print(f"  {error_filename}.png")
print(f"  {error_filename}.eps")

plt.show()

# ============================================================================
# GENERATE GIF ANIMATION (OPTIONAL)
# ============================================================================

if create_gif:
    print("\n" + "=" * 70)
    print("GENERATING GIF ANIMATION")
    print("=" * 70)
    
    temp_frame_files = []
    
    # Determine which timesteps to include
    max_frames = 100
    if len(plot_files) > max_frames:
        step = len(plot_files) // max_frames
        gif_plot_files = plot_files[::step]
    else:
        gif_plot_files = plot_files
    
    print(f"Processing {len(gif_plot_files)} timesteps for GIF...")
    
    for i, plot_file in enumerate(gif_plot_files):
        # Load dataset
        ds_temp = yt.load(plot_file)
        time_temp = float(ds_temp.current_time)
        
        # Extract data
        x_min_temp = float(ds_temp.domain_left_edge[0])
        x_max_temp = float(ds_temp.domain_right_edge[0])
        y_min_temp = float(ds_temp.domain_left_edge[1])
        y_max_temp = float(ds_temp.domain_right_edge[1])
        
        # Determine y-sampling location
        if y_sample_location is None:
            y_sample_temp = (y_min_temp + y_max_temp) / 2.0
        else:
            y_sample_temp = y_sample_location
        
        ray_start_temp = ds_temp.arr([x_min_temp, y_sample_temp, 0.0], 'code_length')
        ray_end_temp = ds_temp.arr([x_max_temp, y_sample_temp, 0.0], 'code_length')
        ray_temp = ds_temp.ray(ray_start_temp, ray_end_temp)
        
        sort_indices_temp = np.argsort(ray_temp['x'])
        x_temp = np.array(ray_temp['x'][sort_indices_temp])
        
        # *** KEY CHANGE: Extract velocityy instead of velocityx ***
        velocity_temp = np.array(ray_temp['velocityy'][sort_indices_temp])
        
        x_dist_temp = x_temp - x_min_temp
        
        # Analytical solution
        x_dist_analytical_temp = np.linspace(0, x_max_temp - x_min_temp, 1000)
        velocity_analytical_temp = analytical_rayleigh_velocity(x_dist_analytical_temp, time_temp, U, nu)
        x_analytical_temp = x_dist_analytical_temp + x_min_temp
        
        # Create frame
        fig_temp, ax_temp = plt.subplots(figsize=(10, 8))
        
        ax_temp.plot(x_analytical_temp, velocity_analytical_temp, 'b-', 
                    linewidth=LINE_WIDTH_EXACT, label='Analytical', zorder=1)
        ax_temp.plot(x_temp, velocity_temp, 'r--', 
                    linewidth=LINE_WIDTH_NUMERICAL, label='Numerical', 
                    zorder=2, alpha=0.8)
        
        ax_temp.set_xlabel('Position x (m)', fontsize=FONT_SIZE_LABEL)
        ax_temp.set_ylabel('Y-Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
        ax_temp.set_title(f"{PLOT_TITLE}\nt = {time_temp:.6e} s", 
                         fontsize=FONT_SIZE_TITLE, fontweight='bold')
        ax_temp.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
        ax_temp.grid(True, alpha=0.3)
        ax_temp.tick_params(labelsize=FONT_SIZE_TICK)
        
        # Consistent axis limits
        ax_temp.set_xlim([x_min_temp, x_max_temp])
        ax_temp.set_ylim([0, U * 1.1])
        
        plt.tight_layout()
        
        # Save temporary frame
        temp_filename = os.path.join(output_folder, f'temp_frame_{i:04d}.png')
        plt.savefig(temp_filename, format='png', dpi=150, bbox_inches='tight')
        temp_frame_files.append(temp_filename)
        plt.close(fig_temp)
        
        if (i + 1) % 10 == 0 or i == len(gif_plot_files) - 1:
            print(f"  Processed {i + 1}/{len(gif_plot_files)} frames")
    
    # Create GIF from frames
    print("\nAssembling GIF...")
    images = [Image.open(frame) for frame in temp_frame_files]
    gif_path = os.path.join(output_folder, gif_filename)
    images[0].save(gif_path, save_all=True, append_images=images[1:], 
                   duration=gif_duration, loop=0)
    
    # Clean up temporary files
    for temp_file in temp_frame_files:
        os.remove(temp_file)
    
    print(f"\nGIF saved: {gif_path}")
    print(f"  Frames: {len(images)}")
    print(f"  Duration per frame: {gif_duration} ms")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nFinal timestep: {simulation_time_final:.6e} s")
print(f"Output directory: {output_folder}")
print(f"\nFiles generated:")
print(f"  - rayleigh_transient_multipanel.png")
print(f"  - rayleigh_transient_multipanel.eps")
print(f"  - rayleigh_final_comparison.png")
print(f"  - rayleigh_final_comparison.eps")
print(f"  - rayleigh_error_semilog.png")
print(f"  - rayleigh_error_semilog.eps")
if create_gif:
    print(f"  - {gif_filename}")
print("\n" + "=" * 70)
