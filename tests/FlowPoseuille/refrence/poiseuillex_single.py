"""
Poiseuille Flow Analysis Script - Pressure Driven IN X DIRECTION

This script extracts velocity profiles from AMReX Poiseuille flow simulations,
compares them with the analytical pressure-driven solution, and generates visualization plots.

Features:
- Extracts 1D velocity profile along y-direction at x=0
- Compares numerical solution with analytical Poiseuille flow (parabolic profile)
- Uses pressure gradient formulation: u_x(y) = -(dP/dx) * (1/2mu) * y(H-y)
- Customizable plotting parameters (fonts, titles, line weights)
- Outputs plots as both .eps and .png formats
- Optional GIF generation showing velocity profile evolution
- Error metrics (L2 and L-infinity norms)

Required inputs:
- gamma: Ratio of specific heats
- density: Fluid density (kg/m^3)
- mu: Dynamic viscosity (Pa-s)
- dpdx: Pressure gradient (Pa/m)
- Channel dimensions from geometry.prob_lo and geometry.prob_hi
"""

import yt
import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image
import glob
import re

# Suppress yt's verbose output
yt.funcs.mylog.setLevel(40)

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================
# Physical parameters (matching your input file)
gamma = 1.4                      # Ratio of specific heats
mu = 0.1                         # Dynamic viscosity (Pa-s)
dpdx = -1.0                      # Pressure gradient (Pa/m)
density = 1.0                  # Fluid density (kg/m^3)

# Channel Dimensions (from your geometry settings)
y_min_sim = -0.5                 # geometry.prob_lo (y-component)
y_max_sim = 0.5                  # geometry.prob_hi (y-component)
H = y_max_sim - y_min_sim        # Total height (1.0 m)

# File paths
amrex_output_dir = r'..\..\..\bin\tests\FlowPoiseuille\FlowPoiseuillex'

# Plotting customization
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_TICK = 11
LINE_WIDTH_EXACT = 2.5
LINE_WIDTH_NUMERICAL = 2.0
PLOT_TITLE = 'Plane Poiseuille Flow: Pressure Driven'

# Output settings
output_folder = './Images'
create_gif = False               # Set to True to generate GIF animation
gif_duration = 200               # Duration per frame in milliseconds
gif_filename = 'poiseuille_evolution.gif'

# ============================================================================
# ANALYTICAL POISEUILLE FLOW SOLUTION (PRESSURE DRIVEN)
# ============================================================================
def analytical_poiseuille_velocity(y_sim, dpdx, mu, H):
    """
    Analytical solution for pressure-driven plane Poiseuille flow.
    
    Formula: u_x(y) = -(dP/dx) * (1/2mu) * y_wall * (H - y_wall)
    
    Where y_wall is the distance from the bottom wall (0 at bottom, H at top).
    
    Parameters:
    -----------
    y_sim : array-like
        Vertical positions in simulation coordinates (e.g., -0.5 to 0.5)
    dpdx : float
        Pressure gradient (Pa/m) - negative for flow in positive x direction
    mu : float
        Dynamic viscosity (Pa-s)
    H : float
        Total channel height (m)
    
    Returns:
    --------
    u : array-like
        Velocity at each y position (m/s)
    """
    # Shift simulation coordinates to wall coordinates (0 to H)
    y_wall = y_sim - y_min_sim
    
    # Apply the pressure-driven Poiseuille formula
    return (-dpdx / (2.0 * mu)) * y_wall * (H - y_wall)

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
print("PLANE POISEUILLE FLOW ANALYSIS - PRESSURE DRIVEN")
print("=" * 70)
print(f"\nPhysical Parameters:")
print(f"  Gamma:              {gamma}")
print(f"  Density:            {density:.2e} kg/m^3")
print(f"  Viscosity (mu):     {mu:.2e} Pa-s")
print(f"  Pressure Gradient:  {dpdx:.2e} Pa/m")
print(f"  Channel Height (H): {H:.4f} m")
print(f"  y_min (sim):        {y_min_sim:.4f} m")
print(f"  y_max (sim):        {y_max_sim:.4f} m")

# Calculate expected maximum velocity (at centerline, y_wall = H/2)
u_max_analytical = (-dpdx / (2.0 * mu)) * (H/2) * (H - H/2)
print(f"  Max Velocity (analytical): {u_max_analytical:.4f} m/s")

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
# EXTRACT DATA FROM LAST TIMESTEP
# ============================================================================
last_plot = plot_files[-1]
print("\n" + "=" * 70)
print("EXTRACTING DATA FROM LAST TIMESTEP")
print("=" * 70)
print(f"Using: {os.path.basename(last_plot)}")
print(f"Timestep number: {extract_timestep_number(last_plot)}")

# Load the dataset
ds = yt.load(last_plot)
simulation_time = float(ds.current_time)
print(f"\nSimulation time: {simulation_time:.6e} s")
print(f"Domain: {ds.domain_left_edge} to {ds.domain_right_edge}")

# Print available fields
print(f"\nAvailable fields:")
for field in ds.field_list:
    print(f"  {field}")

# Create a ray along y-direction at x=0, z=0
print(f"\nExtracting 1D slice at x=0, z=0...")
y_min = float(ds.domain_left_edge[1])
y_max = float(ds.domain_right_edge[1])

# Define the ray from bottom to top at x=0, z=0
ray_start = ds.arr([0.0, y_min, 0.0], 'code_length')
ray_end = ds.arr([0.0, y_max, 0.0], 'code_length')
ray = ds.ray(ray_start, ray_end)

# Sort by y coordinate
sort_indices = np.argsort(ray['y'])
y_numerical = np.array(ray['y'][sort_indices])
velocity_numerical = np.array(ray['velocityx'][sort_indices])

print(f"Extracted {len(y_numerical)} points")

# Print data ranges
print(f"\nNumerical data ranges:")
print(f"  y: [{np.min(y_numerical):.6f}, {np.max(y_numerical):.6f}]")
print(f"  velocityx: [{np.min(velocity_numerical):.6e}, {np.max(velocity_numerical):.6e}]")

# ============================================================================
# COMPUTE ANALYTICAL SOLUTION
# ============================================================================
print("\n" + "=" * 70)
print("COMPUTING ANALYTICAL SOLUTION")
print("=" * 70)

# Create fine grid for analytical solution
y_analytical = np.linspace(y_min, y_max, 1000)
velocity_analytical = analytical_poiseuille_velocity(y_analytical, dpdx, mu, H)

print(f"Analytical solution computed on {len(y_analytical)} points")
print(f"  y (sim coords): [{np.min(y_analytical):.6f}, {np.max(y_analytical):.6f}]")
print(f"  velocity: [{np.min(velocity_analytical):.6e}, {np.max(velocity_analytical):.6e}]")
print(f"  Max velocity (at y=0): {np.max(velocity_analytical):.6e} m/s")

# ============================================================================
# CREATE COMPARISON PLOT
# ============================================================================
print("\n" + "=" * 70)
print("CREATING COMPARISON PLOT")
print("=" * 70)

fig, ax = plt.subplots(figsize=(10, 8))

# Plot analytical solution
ax.plot(velocity_analytical, y_analytical, 'b-', 
        linewidth=LINE_WIDTH_EXACT, label='Analytical (Pressure Driven)', zorder=1)

# Plot numerical solution
ax.plot(velocity_numerical, y_numerical, 'r--', 
        linewidth=LINE_WIDTH_NUMERICAL, label='Numerical Solution', 
        zorder=2, alpha=0.8)

# Formatting
ax.set_xlabel('Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
ax.set_ylabel('Height (m)', fontsize=FONT_SIZE_LABEL)
ax.set_title(f'{PLOT_TITLE}\ndP/dx = {dpdx} Pa/m, mu = {mu} Pa-s\nt = {simulation_time:.6e} s', 
             fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=FONT_SIZE_TICK)

plt.tight_layout()

# Save plots
output_filename = os.path.join(output_folder, 'poiseuille_comparison')
plt.savefig(output_filename + '.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(output_filename + '.eps', format='eps', bbox_inches='tight')

print(f"\nPlots saved:")
print(f"  {output_filename}.png")
print(f"  {output_filename}.eps")

plt.show()

# ============================================================================
# COMPUTE ERROR METRICS
# ============================================================================
print("\n" + "=" * 70)
print("ERROR METRICS")
print("=" * 70)

# Interpolate analytical solution to numerical grid points
velocity_analytical_interp = np.interp(y_numerical, y_analytical, velocity_analytical)

# Calculate errors
velocity_error = velocity_numerical - velocity_analytical_interp
abs_velocity_error = np.abs(velocity_error)
l2_error = np.linalg.norm(velocity_error) / np.sqrt(len(velocity_error))
linf_error = np.max(np.abs(velocity_error))
relative_l2_error = l2_error / (np.linalg.norm(velocity_analytical_interp) / np.sqrt(len(velocity_analytical_interp)))

print(f"\nVelocity Errors:")
print(f"  L2 error:          {l2_error:.6e}")
print(f"  L-infinity error:  {linf_error:.6e}")
print(f"  Relative L2 error: {relative_l2_error:.6e} ({relative_l2_error*100:.4f}%)")

# ============================================================================
# CREATE SEMILOG ERROR PLOT
# ============================================================================
print("\n" + "=" * 70)
print("CREATING SEMILOG ERROR PLOT")
print("=" * 70)

fig, ax = plt.subplots(figsize=(10, 8))

# Semilog plot of absolute error
epsilon = 1e-16
abs_error_safe = abs_velocity_error + epsilon

ax.semilogy(y_numerical, abs_error_safe, 'k-', linewidth=LINE_WIDTH_NUMERICAL)

ax.set_xlabel('Height (m)', fontsize=FONT_SIZE_LABEL)
ax.set_ylabel('Absolute Error |vx - v_true| (m/s)', fontsize=FONT_SIZE_LABEL)
ax.set_title(f'Poiseuille Flow: Absolute Error\nt = {simulation_time:.6e} s', 
             fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
ax.tick_params(labelsize=FONT_SIZE_TICK)

plt.tight_layout()

# Save error plot
error_filename = os.path.join(output_folder, 'poiseuille_error_semilog')
plt.savefig(error_filename + '.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(error_filename + '.eps', format='eps', bbox_inches='tight')

print(f"\nError plot saved:")
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
    
    gif_frames = []
    temp_frame_files = []
    
    # Determine which timesteps to include
    max_frames = 100
    if len(plot_files) > max_frames:
        step = len(plot_files) // max_frames
        selected_plots = plot_files[::step]
    else:
        selected_plots = plot_files
    
    print(f"Processing {len(selected_plots)} timesteps for GIF...")
    
    for i, plot_file in enumerate(selected_plots):
        # Load dataset
        ds_temp = yt.load(plot_file)
        time_temp = float(ds_temp.current_time)
        
        # Extract data
        ray_start_temp = ds_temp.arr([0.0, y_min, 0.0], 'code_length')
        ray_end_temp = ds_temp.arr([0.0, y_max, 0.0], 'code_length')
        ray_temp = ds_temp.ray(ray_start_temp, ray_end_temp)
        
        sort_indices_temp = np.argsort(ray_temp['y'])
        y_temp = np.array(ray_temp['y'][sort_indices_temp])
        velocity_temp = np.array(ray_temp['velocityx'][sort_indices_temp])
        
        # Create frame
        fig_temp, ax_temp = plt.subplots(figsize=(10, 8))
        
        ax_temp.plot(velocity_analytical, y_analytical, 'b-', 
                    linewidth=LINE_WIDTH_EXACT, label='Analytical', zorder=1)
        ax_temp.plot(velocity_temp, y_temp, 'r--', 
                    linewidth=LINE_WIDTH_NUMERICAL, label='Numerical', 
                    zorder=2, alpha=0.8)
        
        ax_temp.set_xlabel('Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
        ax_temp.set_ylabel('Height (m)', fontsize=FONT_SIZE_LABEL)
        ax_temp.set_title(f'{PLOT_TITLE}\nt = {time_temp:.6e} s', 
                         fontsize=FONT_SIZE_TITLE, fontweight='bold')
        ax_temp.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
        ax_temp.grid(True, alpha=0.3)
        ax_temp.tick_params(labelsize=FONT_SIZE_TICK)
        
        # Set consistent axis limits
        max_vel = np.max(velocity_analytical)
        ax_temp.set_xlim([0, max_vel * 1.1])
        ax_temp.set_ylim([y_min, y_max])
        
        plt.tight_layout()
        
        # Save temporary frame
        temp_filename = os.path.join(output_folder, f'temp_frame_{i:04d}.png')
        plt.savefig(temp_filename, format='png', dpi=150, bbox_inches='tight')
        temp_frame_files.append(temp_filename)
        plt.close(fig_temp)
        
        if (i + 1) % 10 == 0 or i == len(selected_plots) - 1:
            print(f"  Processed {i + 1}/{len(selected_plots)} frames")
    
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
print(f"\nFinal timestep: {simulation_time:.6e} s")
print(f"Output directory: {output_folder}")
print(f"\nFiles generated:")
print(f"  - poiseuille_comparison.png")
print(f"  - poiseuille_comparison.eps")
print(f"  - poiseuille_error_semilog.png")
print(f"  - poiseuille_error_semilog.eps")
if create_gif:
    print(f"  - {gif_filename}")
print("\n" + "=" * 70)
