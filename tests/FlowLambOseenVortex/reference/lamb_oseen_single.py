"""
===============================================================================
LAMB-OSEEN VORTEX COMPREHENSIVE ANALYSIS SCRIPT
===============================================================================

PURPOSE:
    Perform full diagnostic of Lamb-Oseen vortex simulations, comparing 
    Domain KE, velocity profiles (U and V), and 3D error surfaces against 
    analytical solutions. Includes effective viscosity (mu_eff) tracking.

FEATURES:
    - 9 standalone diagnostic plots (.png and .eps)
    - Evolution GIFs for velocity profiles and 3D error surfaces
    - Static color mapping for 3D error surfaces across all timesteps

INPUTS:
    - AMReX plot files from Lamb-Oseen simulation
    - Physical parameters: rho, mu, Gamma, rc0, xc, yc

OUTPUTS:
    - 9 diagnostic plots 
    - 6 evolution GIFs 

===============================================================================
"""

import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
import glob
import re
from PIL import Image

# Suppress yt verbose output
yt.funcs.mylog.setLevel(40)

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

# Physical parameters (matching Lamb-Oseen input file)
gamma = 1.4
mu = 0.01                        # Dynamic viscosity (Pa-s)
density = 1.0                    # Fluid density (kg/m^3)
nu = mu / density                # Kinematic viscosity (m^2/s)

# Lamb-Oseen Specifics
Gamma = 1.0                      # Circulation
rc0 = 0.5                        # Initial core radius
xc, yc = 5.0, 5.0                # Center of the vortex
L = 10.0                         # Domain length (m)

# File paths
amrex_output_dir = r'..\..\..\bin\tests\FlowLambOseenVortex\FlowLambOseenVortex'

# Plotting customization
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_TICK = 11
LINE_WIDTH_EXACT = 2.5
LINE_WIDTH_NUMERICAL = 2.0

# Output settings
output_folder = './LambOseen_Analysis'
create_gifs = False               # Set to True to generate evolution GIFs
gif_duration = 200               # Duration per frame in milliseconds

# Sampling location for velocity profiles (through the center)
x_sample = xc                 

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ============================================================================
# ANALYTICAL SOLUTIONS
# ============================================================================

def get_lamb_oseen_vtheta(r, t, Gamma, rc0, nu):
    """Tangential velocity: v_theta = (Gamma / (2*pi*r)) * (1 - exp(-r^2 / (rc0^2 + 4*nu*t)))"""
    if t < 0: 
        t = 0
    r_core_sq = rc0**2 + 4.0 * nu * t
    # Avoid division by zero at r=0
    return (Gamma / (2.0 * np.pi * (r + 1e-12))) * (1.0 - np.exp(-(r**2) / r_core_sq))

def analytical_lamb_oseen_u(x, y, t, Gamma, rc0, nu, xc, yc):
    """U-velocity: u = -v_theta * (y - yc) / r"""
    r = np.sqrt((x - xc)**2 + (y - yc)**2)
    v_theta = get_lamb_oseen_vtheta(r, t, Gamma, rc0, nu)
    return -v_theta * (y - yc) / (r + 1e-12)

def analytical_lamb_oseen_v(x, y, t, Gamma, rc0, nu, xc, yc):
    """V-velocity: v = v_theta * (x - xc) / r"""
    r = np.sqrt((x - xc)**2 + (y - yc)**2)
    v_theta = get_lamb_oseen_vtheta(r, t, Gamma, rc0, nu)
    return v_theta * (x - xc) / (r + 1e-12)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def extract_timestep_number(filename):
    """Extract timestep number from plot file name"""
    match = re.search(r'(\d+)', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return 0

# ============================================================================
# FIND AND SORT PLOT FILES
# ============================================================================

print("=" * 70)
print("LAMB-OSEEN VORTEX COMPREHENSIVE ANALYSIS")
print("=" * 70)

plot_files = []
for item in os.listdir(amrex_output_dir):
    item_path = os.path.join(amrex_output_dir, item)
    if os.path.isdir(item_path) and item.endswith('cell'):
        plot_files.append(item_path)

if not plot_files:
    print(f"ERROR: No plot files found in {amrex_output_dir}")
    exit(1)

plot_files.sort(key=extract_timestep_number)
print(f"\nFound {len(plot_files)} plot files")

# ============================================================================
# EXTRACT DATA FROM ALL TIMESTEPS
# ============================================================================

print("\n" + "=" * 70)
print("EXTRACTING DATA FROM ALL TIMESTEPS")
print("=" * 70)

times = []
kinetic_energies = []
all_u_data = []
all_v_data = []
all_y_coords = []
analytical_kes = []

for i, plot_file in enumerate(plot_files):
    ds = yt.load(plot_file)
    t = float(ds.current_time)
    
    # Extract all data for KE calculation
    ad = ds.all_data()
    rho = np.array(ad['density'])
    vx = np.array(ad['velocityx'])
    vy = np.array(ad['velocityy'])
    
    # Calculate domain kinetic energy
    KE = 0.5 * np.sum(rho * (vx**2 + vy**2))
    
    times.append(t)
    kinetic_energies.append(KE)
    
    # Calculate analytical KE by integrating analytical profile over domain
    x_coords = np.array(ad['x'])
    y_coords = np.array(ad['y'])
    u_ana = analytical_lamb_oseen_u(x_coords, y_coords, t, Gamma, rc0, nu, xc, yc)
    v_ana = analytical_lamb_oseen_v(x_coords, y_coords, t, Gamma, rc0, nu, xc, yc)
    analytical_kes.append(0.5 * np.sum(density * (u_ana**2 + v_ana**2)))
    
    # Extract 1D profile at x_sample
    y_min = float(ds.domain_left_edge[1])
    y_max = float(ds.domain_right_edge[1])
    
    ray_start = ds.arr([x_sample, y_min, 0.0], 'code_length')
    ray_end = ds.arr([x_sample, y_max, 0.0], 'code_length')
    ray = ds.ray(ray_start, ray_end)
    
    sort_indices = np.argsort(ray['y'])
    y_coords_1d = np.array(ray['y'][sort_indices])
    u_profile = np.array(ray['velocityx'][sort_indices])
    v_profile = np.array(ray['velocityy'][sort_indices])
    
    all_u_data.append(u_profile)
    all_v_data.append(v_profile)
    all_y_coords.append(y_coords_1d)
    
    if (i + 1) % 10 == 0 or i == len(plot_files) - 1:
        print(f"  Processed {i + 1}/{len(plot_files)} timesteps")

times = np.array(times)
kinetic_energies = np.array(kinetic_energies)
analytical_kes = np.array(analytical_kes)

# ============================================================================
# CALCULATE EFFECTIVE VISCOSITY
# ============================================================================

print("\n" + "=" * 70)
print("CALCULATING EFFECTIVE VISCOSITY")
print("=" * 70)

# Calculate mu_eff from KE decay rate
# For Lamb-Oseen, use 1/rc as characteristic wavenumber
k_char = 1.0 / rc0

ln_ke = np.log(kinetic_energies + 1e-16)
d_ln_ke_dt = np.gradient(ln_ke, times)
nu_eff_time = -d_ln_ke_dt / (2.0 * k_char**2)
mu_eff_time = nu_eff_time * density

if len(times) > 10:
    fit_start = 5
    log_KE = np.log(kinetic_energies[fit_start:] + 1e-16)
    times_fit = times[fit_start:]
    
    coeffs = np.polyfit(times_fit, log_KE, 1)
    slope_measured = coeffs[0]
    
    nu_effective = -slope_measured / (2.0 * k_char**2)
    mu_effective = nu_effective * density
    
    print(f"\nViscosity Analysis:")
    print(f"  Expected mu:    {mu:.6e} Pa-s")
    print(f"  Effective mu:   {mu_effective:.6e} Pa-s")
    print(f"  Ratio:          {mu_effective/mu:.6f}")

# ============================================================================
# SELECT 6 TIMESTEPS FOR 2x3 GRIDS
# ============================================================================

num_panels = 6
if len(plot_files) <= num_panels:
    selected_indices = list(range(len(plot_files)))
else:
    selected_indices = np.linspace(0, len(plot_files) - 1, num_panels, dtype=int)

selected_times = times[selected_indices]
print(f"\nSelected {num_panels} timesteps for 2x3 grids:")
for i, idx in enumerate(selected_indices):
    print(f"  Panel {i+1}: t = {times[idx]:.6e} s")

# ============================================================================
# PLOT 1: DOMAIN KE VS ANALYTICAL KE (LINEAR-LINEAR)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 1: KE COMPARISON (LINEAR-LINEAR)")
print("=" * 70)

fig1, ax1 = plt.subplots(figsize=(10, 8))
ax1.plot(times, analytical_kes, 'b-', linewidth=LINE_WIDTH_EXACT, 
         label='Analytical KE', zorder=1)
ax1.plot(times, kinetic_energies, 'ro', markersize=6, 
         label='Domain KE', alpha=0.7, zorder=2)
ax1.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
ax1.set_ylabel('Kinetic Energy', fontsize=FONT_SIZE_LABEL)
ax1.set_title(f'Lamb-Oseen Vortex: KE Comparison\nmu = {mu} Pa-s, Gamma = {Gamma}, rc0 = {rc0}', 
              fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax1.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
ax1.grid(True, alpha=0.3)
ax1.tick_params(labelsize=FONT_SIZE_TICK)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, '01_KE_Comparison.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '01_KE_Comparison.eps'))
print("  Saved: 01_KE_Comparison.png/.eps")
plt.close()

# ============================================================================
# PLOT 2: KE ERROR (SEMILOG)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 2: KE ERROR (SEMILOG)")
print("=" * 70)

fig2, ax2 = plt.subplots(figsize=(10, 8))
ke_error = np.abs(kinetic_energies - analytical_kes)
epsilon = 1e-16
ke_error_safe = ke_error + epsilon

ax2.semilogy(times, ke_error_safe, 'k-', linewidth=LINE_WIDTH_NUMERICAL)
ax2.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
ax2.set_ylabel('|Domain KE - Analytical KE|', fontsize=FONT_SIZE_LABEL)
ax2.set_title('Lamb-Oseen Vortex: KE Error', 
              fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')
ax2.tick_params(labelsize=FONT_SIZE_TICK)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, '02_KE_Error.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '02_KE_Error.eps'))
print("  Saved: 02_KE_Error.png/.eps")
plt.close()

# ============================================================================
# PLOT 3: U VELOCITY COMPARISON (2x3 GRID)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 3: U VELOCITY COMPARISON (2x3 GRID)")
print("=" * 70)

fig3, axes3 = plt.subplots(2, 3, figsize=(15, 10))
axes3 = axes3.flatten()

for i, idx in enumerate(selected_indices):
    t = times[idx]
    y_h = all_y_coords[idx]
    u_num = all_u_data[idx]
    u_ana = analytical_lamb_oseen_u(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
    
    axes3[i].plot(u_ana, y_h, 'b-', linewidth=LINE_WIDTH_EXACT, 
                  label='Analytical', zorder=1)
    axes3[i].plot(u_num, y_h, 'r--', linewidth=LINE_WIDTH_NUMERICAL, 
                  label='Numerical', zorder=2, alpha=0.8)
    
    axes3[i].set_xlabel('U Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
    axes3[i].set_ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
    axes3[i].set_title(f't = {t:.4e} s', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    axes3[i].legend(fontsize=FONT_SIZE_LEGEND, loc='best')
    axes3[i].grid(True, alpha=0.3)
    axes3[i].tick_params(labelsize=FONT_SIZE_TICK)

fig3.suptitle(f'U Velocity Comparison at x = {x_sample:.4f} m', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(output_folder, '03_U_Comparison.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '03_U_Comparison.eps'))
print("  Saved: 03_U_Comparison.png/.eps")
plt.close()

# ============================================================================
# PLOT 4: U VELOCITY ERROR (2x3 GRID)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 4: U VELOCITY ERROR (2x3 GRID)")
print("=" * 70)

fig4, axes4 = plt.subplots(2, 3, figsize=(15, 10))
axes4 = axes4.flatten()

for i, idx in enumerate(selected_indices):
    t = times[idx]
    y_h = all_y_coords[idx]
    u_num = all_u_data[idx]
    u_ana = analytical_lamb_oseen_u(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
    u_error = np.abs(u_num - u_ana)
    
    axes4[i].plot(u_error, y_h, 'k-', linewidth=LINE_WIDTH_NUMERICAL)
    
    axes4[i].set_xlabel('|U - U_analytical| (m/s)', fontsize=FONT_SIZE_LABEL)
    axes4[i].set_ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
    axes4[i].set_title(f't = {t:.4e} s', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    axes4[i].grid(True, alpha=0.3)
    axes4[i].tick_params(labelsize=FONT_SIZE_TICK)

fig4.suptitle(f'U Velocity Error at x = {x_sample:.4f} m', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(output_folder, '04_U_Error.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '04_U_Error.eps'))
print("  Saved: 04_U_Error.png/.eps")
plt.close()

# ============================================================================
# PLOT 5: V VELOCITY COMPARISON (2x3 GRID)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 5: V VELOCITY COMPARISON (2x3 GRID)")
print("=" * 70)

fig5, axes5 = plt.subplots(2, 3, figsize=(15, 10))
axes5 = axes5.flatten()

for i, idx in enumerate(selected_indices):
    t = times[idx]
    y_h = all_y_coords[idx]
    v_num = all_v_data[idx]
    v_ana = analytical_lamb_oseen_v(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
    
    axes5[i].plot(v_ana, y_h, 'b-', linewidth=LINE_WIDTH_EXACT, 
                  label='Analytical', zorder=1)
    axes5[i].plot(v_num, y_h, 'r--', linewidth=LINE_WIDTH_NUMERICAL, 
                  label='Numerical', zorder=2, alpha=0.8)
    
    axes5[i].set_xlabel('V Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
    axes5[i].set_ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
    axes5[i].set_title(f't = {t:.4e} s', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    axes5[i].legend(fontsize=FONT_SIZE_LEGEND, loc='best')
    axes5[i].grid(True, alpha=0.3)
    axes5[i].tick_params(labelsize=FONT_SIZE_TICK)

fig5.suptitle(f'V Velocity Comparison at x = {x_sample:.4f} m', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(output_folder, '05_V_Comparison.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '05_V_Comparison.eps'))
print("  Saved: 05_V_Comparison.png/.eps")
plt.close()

# ============================================================================
# PLOT 6: V VELOCITY ERROR (2x3 GRID)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 6: V VELOCITY ERROR (2x3 GRID)")
print("=" * 70)

fig6, axes6 = plt.subplots(2, 3, figsize=(15, 10))
axes6 = axes6.flatten()

for i, idx in enumerate(selected_indices):
    t = times[idx]
    y_h = all_y_coords[idx]
    v_num = all_v_data[idx]
    v_ana = analytical_lamb_oseen_v(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
    v_error = np.abs(v_num - v_ana)
    
    axes6[i].plot(v_error, y_h, 'k-', linewidth=LINE_WIDTH_NUMERICAL)
    
    axes6[i].set_xlabel('|V - V_analytical| (m/s)', fontsize=FONT_SIZE_LABEL)
    axes6[i].set_ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
    axes6[i].set_title(f't = {t:.4e} s', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    axes6[i].grid(True, alpha=0.3)
    axes6[i].tick_params(labelsize=FONT_SIZE_TICK)

fig6.suptitle(f'V Velocity Error at x = {x_sample:.4f} m', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(output_folder, '06_V_Error.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '06_V_Error.eps'))
print("  Saved: 06_V_Error.png/.eps")
plt.close()

# ============================================================================
# EXTRACT 2D FIELDS FOR 3D ERROR PLOTS
# ============================================================================

print("\n" + "=" * 70)
print("EXTRACTING 2D FIELDS FOR 3D ERROR PLOTS")
print("=" * 70)

all_u_fields = []
all_v_fields = []
X_grid = None
Y_grid = None

for i, idx in enumerate(selected_indices):
    ds = yt.load(plot_files[idx])
    
    # Create 2D slice at z=0
    slc = ds.slice('z', 0.0)
    frb = slc.to_frb((L, 'code_length'), 128)
    
    u_field = np.array(frb['velocityx'])
    v_field = np.array(frb['velocityy'])
    
    if X_grid is None:
        x_1d = np.linspace(0, L, u_field.shape[1])
        y_1d = np.linspace(0, L, u_field.shape[0])
        X_grid, Y_grid = np.meshgrid(x_1d, y_1d)
    
    all_u_fields.append(u_field)
    all_v_fields.append(v_field)
    
    print(f"  Extracted 2D field for timestep {i+1}/{num_panels}")

# Calculate global min/max for static color mapping
u_errors_all = []
v_errors_all = []

for i, idx in enumerate(selected_indices):
    t = times[idx]
    u_ana_field = analytical_lamb_oseen_u(X_grid, Y_grid, t, Gamma, rc0, nu, xc, yc)
    v_ana_field = analytical_lamb_oseen_v(X_grid, Y_grid, t, Gamma, rc0, nu, xc, yc)
    
    u_error_field = np.abs(all_u_fields[i] - u_ana_field)
    v_error_field = np.abs(all_v_fields[i] - v_ana_field)
    
    u_errors_all.append(u_error_field)
    v_errors_all.append(v_error_field)

u_error_min = min([np.min(e) for e in u_errors_all])
u_error_max = max([np.max(e) for e in u_errors_all])
v_error_min = min([np.min(e) for e in v_errors_all])
v_error_max = max([np.max(e) for e in v_errors_all])

print(f"\nU Error range: [{u_error_min:.6e}, {u_error_max:.6e}]")
print(f"V Error range: [{v_error_min:.6e}, {v_error_max:.6e}]")

# ============================================================================
# PLOT 7: 3D U ERROR SURFACES (2x3 GRID, SEMILOG Z-AXIS)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 7: 3D U ERROR SURFACES (2x3 GRID)")
print("=" * 70)

fig7 = plt.figure(figsize=(18, 12))

for i in range(num_panels):
    ax = fig7.add_subplot(2, 3, i+1, projection='3d')
    
    # Use log scale for Z
    u_error_log = np.log10(u_errors_all[i] + 1e-16)
    
    surf = ax.plot_surface(X_grid, Y_grid, u_error_log, 
                           cmap=cm.viridis, 
                           vmin=np.log10(u_error_min + 1e-16),
                           vmax=np.log10(u_error_max + 1e-16),
                           linewidth=0, antialiased=True)
    
    ax.set_xlabel('X (m)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Y (m)', fontsize=FONT_SIZE_LABEL)
    ax.set_zlabel('log10(|U Error|)', fontsize=FONT_SIZE_LABEL)
    ax.set_title(f't = {selected_times[i]:.4e} s', fontsize=FONT_SIZE_TITLE)

fig7.suptitle('3D U Velocity Error Surfaces (Semilog Z)', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_folder, '07_U_3D_Error.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '07_U_3D_Error.eps'))
print("  Saved: 07_U_3D_Error.png/.eps")
plt.close()

# ============================================================================
# PLOT 8: 3D V ERROR SURFACES (2x3 GRID, SEMILOG Z-AXIS)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 8: 3D V ERROR SURFACES (2x3 GRID)")
print("=" * 70)

fig8 = plt.figure(figsize=(18, 12))

for i in range(num_panels):
    ax = fig8.add_subplot(2, 3, i+1, projection='3d')
    
    # Use log scale for Z
    v_error_log = np.log10(v_errors_all[i] + 1e-16)
    
    surf = ax.plot_surface(X_grid, Y_grid, v_error_log, 
                           cmap=cm.viridis, 
                           vmin=np.log10(v_error_min + 1e-16),
                           vmax=np.log10(v_error_max + 1e-16),
                           linewidth=0, antialiased=True)
    
    ax.set_xlabel('X (m)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Y (m)', fontsize=FONT_SIZE_LABEL)
    ax.set_zlabel('log10(|V Error|)', fontsize=FONT_SIZE_LABEL)
    ax.set_title(f't = {selected_times[i]:.4e} s', fontsize=FONT_SIZE_TITLE)

fig8.suptitle('3D V Velocity Error Surfaces (Semilog Z)', 
              fontsize=FONT_SIZE_TITLE + 2, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_folder, '08_V_3D_Error.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '08_V_3D_Error.eps'))
print("  Saved: 08_V_3D_Error.png/.eps")
plt.close()

# ============================================================================
# PLOT 9: MU_EFF VS TIME (LINEAR-LINEAR)
# ============================================================================

print("\n" + "=" * 70)
print("CREATING PLOT 9: MU_EFF VS TIME (LINEAR-LINEAR)")
print("=" * 70)

fig9, ax9 = plt.subplots(figsize=(10, 8))
ax9.plot(times, mu_eff_time, 'g-', linewidth=LINE_WIDTH_NUMERICAL, 
         label='Numerical mu_eff')
ax9.axhline(y=mu, color='b', linestyle='--', linewidth=LINE_WIDTH_EXACT, 
            label=f'Analytical mu = {mu} Pa-s')
ax9.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
ax9.set_ylabel('Effective Viscosity mu_eff (Pa-s)', fontsize=FONT_SIZE_LABEL)
ax9.set_title('Effective Viscosity vs Time', 
              fontsize=FONT_SIZE_TITLE, fontweight='bold')
ax9.legend(fontsize=FONT_SIZE_LEGEND, loc='best')
ax9.grid(True, alpha=0.3)
ax9.tick_params(labelsize=FONT_SIZE_TICK)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, '09_Mu_Eff_Time.png'), dpi=300)
plt.savefig(os.path.join(output_folder, '09_Mu_Eff_Time.eps'))
print("  Saved: 09_Mu_Eff_Time.png/.eps")
plt.close()

# ============================================================================
# GIF GENERATION (PLOTS 3-8)
# ============================================================================

if create_gifs:
    print("\n" + "=" * 70)
    print("GENERATING EVOLUTION GIFS")
    print("=" * 70)
    
    # Subsample timesteps for GIF
    max_gif_frames = 50
    if len(plot_files) > max_gif_frames:
        gif_step = len(plot_files) // max_gif_frames
        gif_indices = list(range(0, len(plot_files), gif_step))
    else:
        gif_indices = list(range(len(plot_files)))
    
    print(f"  Using {len(gif_indices)} frames for GIFs")
    
    # GIF 1: U Velocity Comparison Evolution
    print("\n  Creating GIF 1: U Velocity Comparison Evolution...")
    temp_frames = []
    for idx in gif_indices:
        t = times[idx]
        y_h = all_y_coords[idx]
        u_num = all_u_data[idx]
        u_ana = analytical_lamb_oseen_u(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
        
        fig_temp = plt.figure(figsize=(10, 8))
        plt.plot(u_ana, y_h, 'b-', linewidth=LINE_WIDTH_EXACT, label='Analytical')
        plt.plot(u_num, y_h, 'r--', linewidth=LINE_WIDTH_NUMERICAL, label='Numerical')
        plt.xlabel('U Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
        plt.ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
        plt.title(f'U Velocity Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        plt.legend(fontsize=FONT_SIZE_LEGEND)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_u_comp.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_03_U_Comparison.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_03_U_Comparison.gif")
    
    # GIF 2: U Velocity Error Evolution
    print("  Creating GIF 2: U Velocity Error Evolution...")
    temp_frames = []
    for idx in gif_indices:
        t = times[idx]
        y_h = all_y_coords[idx]
        u_num = all_u_data[idx]
        u_ana = analytical_lamb_oseen_u(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
        u_error = np.abs(u_num - u_ana)
        
        fig_temp = plt.figure(figsize=(10, 8))
        plt.plot(u_error, y_h, 'k-', linewidth=LINE_WIDTH_NUMERICAL)
        plt.xlabel('|U - U_analytical| (m/s)', fontsize=FONT_SIZE_LABEL)
        plt.ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
        plt.title(f'U Velocity Error Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_u_err.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_04_U_Error.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_04_U_Error.gif")
    
    # GIF 3: V Velocity Comparison Evolution
    print("  Creating GIF 3: V Velocity Comparison Evolution...")
    temp_frames = []
    for idx in gif_indices:
        t = times[idx]
        y_h = all_y_coords[idx]
        v_num = all_v_data[idx]
        v_ana = analytical_lamb_oseen_v(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
        
        fig_temp = plt.figure(figsize=(10, 8))
        plt.plot(v_ana, y_h, 'b-', linewidth=LINE_WIDTH_EXACT, label='Analytical')
        plt.plot(v_num, y_h, 'r--', linewidth=LINE_WIDTH_NUMERICAL, label='Numerical')
        plt.xlabel('V Velocity (m/s)', fontsize=FONT_SIZE_LABEL)
        plt.ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
        plt.title(f'V Velocity Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        plt.legend(fontsize=FONT_SIZE_LEGEND)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_v_comp.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_05_V_Comparison.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_05_V_Comparison.gif")
    
    # GIF 4: V Velocity Error Evolution
    print("  Creating GIF 4: V Velocity Error Evolution...")
    temp_frames = []
    for idx in gif_indices:
        t = times[idx]
        y_h = all_y_coords[idx]
        v_num = all_v_data[idx]
        v_ana = analytical_lamb_oseen_v(x_sample, y_h, t, Gamma, rc0, nu, xc, yc)
        v_error = np.abs(v_num - v_ana)
        
        fig_temp = plt.figure(figsize=(10, 8))
        plt.plot(v_error, y_h, 'k-', linewidth=LINE_WIDTH_NUMERICAL)
        plt.xlabel('|V - V_analytical| (m/s)', fontsize=FONT_SIZE_LABEL)
        plt.ylabel('Height Y (m)', fontsize=FONT_SIZE_LABEL)
        plt.title(f'V Velocity Error Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_v_err.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_06_V_Error.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_06_V_Error.gif")
    
    # GIF 5: 3D U Error Surface Evolution (with static color mapping)
    print("  Creating GIF 5: 3D U Error Surface Evolution...")
    temp_frames = []
    
    # Pre-extract all 2D fields for GIF
    gif_u_fields = []
    gif_times = []
    for idx in gif_indices:
        ds = yt.load(plot_files[idx])
        slc = ds.slice('z', 0.0)
        frb = slc.to_frb((L, 'code_length'), 128)
        u_field = np.array(frb['velocityx'])
        gif_u_fields.append(u_field)
        gif_times.append(times[idx])
    
    for i, idx in enumerate(gif_indices):
        t = gif_times[i]
        u_ana_field = analytical_lamb_oseen_u(X_grid, Y_grid, t, Gamma, rc0, nu, xc, yc)
        u_error_field = np.abs(gif_u_fields[i] - u_ana_field)
        u_error_log = np.log10(u_error_field + 1e-16)
        
        fig_temp = plt.figure(figsize=(12, 9))
        ax = fig_temp.add_subplot(111, projection='3d')
        
        surf = ax.plot_surface(X_grid, Y_grid, u_error_log, 
                               cmap=cm.viridis,
                               vmin=np.log10(u_error_min + 1e-16),
                               vmax=np.log10(u_error_max + 1e-16),
                               linewidth=0, antialiased=True)
        
        ax.set_xlabel('X (m)', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('Y (m)', fontsize=FONT_SIZE_LABEL)
        ax.set_zlabel('log10(|U Error|)', fontsize=FONT_SIZE_LABEL)
        ax.set_title(f'3D U Error Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        ax.set_zlim([np.log10(u_error_min + 1e-16), np.log10(u_error_max + 1e-16)])
        
        fig_temp.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_u_3d.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_07_U_3D_Error.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_07_U_3D_Error.gif")
    
    # GIF 6: 3D V Error Surface Evolution (with static color mapping)
    print("  Creating GIF 6: 3D V Error Surface Evolution...")
    temp_frames = []
    
    # Pre-extract all 2D fields for GIF
    gif_v_fields = []
    for idx in gif_indices:
        ds = yt.load(plot_files[idx])
        slc = ds.slice('z', 0.0)
        frb = slc.to_frb((L, 'code_length'), 128)
        v_field = np.array(frb['velocityy'])
        gif_v_fields.append(v_field)
    
    for i, idx in enumerate(gif_indices):
        t = gif_times[i]
        v_ana_field = analytical_lamb_oseen_v(X_grid, Y_grid, t, Gamma, rc0, nu, xc, yc)
        v_error_field = np.abs(gif_v_fields[i] - v_ana_field)
        v_error_log = np.log10(v_error_field + 1e-16)
        
        fig_temp = plt.figure(figsize=(12, 9))
        ax = fig_temp.add_subplot(111, projection='3d')
        
        surf = ax.plot_surface(X_grid, Y_grid, v_error_log, 
                               cmap=cm.viridis,
                               vmin=np.log10(v_error_min + 1e-16),
                               vmax=np.log10(v_error_max + 1e-16),
                               linewidth=0, antialiased=True)
        
        ax.set_xlabel('X (m)', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('Y (m)', fontsize=FONT_SIZE_LABEL)
        ax.set_zlabel('log10(|V Error|)', fontsize=FONT_SIZE_LABEL)
        ax.set_title(f'3D V Error Evolution\nt = {t:.6e} s', fontsize=FONT_SIZE_TITLE)
        ax.set_zlim([np.log10(v_error_min + 1e-16), np.log10(v_error_max + 1e-16)])
        
        fig_temp.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
        plt.tight_layout()
        
        temp_file = os.path.join(output_folder, 'temp_v_3d.png')
        plt.savefig(temp_file, dpi=150)
        temp_frames.append(Image.open(temp_file))
        plt.close()
    
    gif_path = os.path.join(output_folder, 'GIF_08_V_3D_Error.gif')
    temp_frames[0].save(gif_path, save_all=True, append_images=temp_frames[1:], 
                        duration=gif_duration, loop=0)
    print(f"    Saved: GIF_08_V_3D_Error.gif")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nOutput directory: {output_folder}")
print(f"\nFiles generated:")
print(f"  - 01_KE_Comparison.png/.eps")
print(f"  - 02_KE_Error.png/.eps")
print(f"  - 03_U_Comparison.png/.eps")
print(f"  - 04_U_Error.png/.eps")
print(f"  - 05_V_Comparison.png/.eps")
print(f"  - 06_V_Error.png/.eps")
print(f"  - 07_U_3D_Error.png/.eps")
print(f"  - 08_V_3D_Error.png/.eps")
print(f"  - 09_Mu_Eff_Time.png/.eps")

if create_gifs:
    print(f"  - GIF_03_U_Comparison.gif")
    print(f"  - GIF_04_U_Error.gif")
    print(f"  - GIF_05_V_Comparison.gif")
    print(f"  - GIF_06_V_Error.gif")
    print(f"  - GIF_07_U_3D_Error.gif")
    print(f"  - GIF_08_V_3D_Error.gif")

print("\n" + "=" * 70)
