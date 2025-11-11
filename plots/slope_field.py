"""
Create slope field visualizations showing particle dynamics as height changes
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
import processing.friends as friend
from scipy.constants import g


def create_slope_field(save_path="slope_field.png"):
    """
    Create slope fields showing:
    1. dv/dh as function of (height, velocity) for different radii
    2. dr/dh as function of (height, radius) for different velocities
    """

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # ===== LEFT PLOT: Velocity vs Height Slope Field =====
    # Grid for height and velocity
    heights = np.linspace(100, 10, 20)  # km
    velocities = np.linspace(0.5, 20, 20)  # km/s

    H, V = np.meshgrid(heights, velocities)

    # Calculate slopes for several different particle radii
    radii_to_plot = [10, 20, 50, 100]  # microns
    colors_r = ['blue', 'green', 'orange', 'red']

    for radius, color in zip(radii_to_plot, colors_r):
        DV_DH = np.zeros_like(H)

        for i in range(H.shape[0]):
            for j in range(H.shape[1]):
                h_km = H[i, j]
                v_kms = V[i, j]

                particle = {'r': radius, 'v': v_kms, 'active': True}
                m = friend.calculate_mass(particle)
                drag = friend.calculate_drag(particle, h_km * 1000)

                v_ms = v_kms * 1000
                if v_ms > 10:
                    dv_dh = ((drag / m) - g) / v_ms * 1000  # (m/s²) / (m/s) * 1000m/km
                    DV_DH[i, j] = dv_dh
                else:
                    DV_DH[i, j] = 0

        # Plot slope field for this radius
        # Normalize slopes for consistent arrow lengths
        DH = np.ones_like(DV_DH)  # Change in h is always 1 (descending)

        # Normalize for visualization
        magnitude = np.sqrt(DH**2 + DV_DH**2)
        max_mag = np.nanmax(magnitude)
        if max_mag > 0:
            scale_factor = 2.0 / max_mag
            DH_norm = DH * scale_factor
            DV_norm = DV_DH * scale_factor
        else:
            DH_norm = DH
            DV_norm = DV_DH

        # Plot with reduced density for clarity
        skip = 1
        ax1.quiver(H[::skip, ::skip], V[::skip, ::skip],
                   -DH_norm[::skip, ::skip], DV_norm[::skip, ::skip],  # negative DH because descending
                   color=color, alpha=0.6, scale=20, width=0.003,
                   label=f'r={radius}μm')

    ax1.set_xlabel('Height (km)', fontsize=13)
    ax1.set_ylabel('Velocity (km/s)', fontsize=13)
    ax1.set_title('Slope Field: dv/dh for Different Particle Sizes', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(10, 100)
    ax1.set_ylim(0, 20)

    # ===== RIGHT PLOT: Radius vs Height Slope Field =====
    # Grid for height and radius
    heights_r = np.linspace(100, 10, 20)  # km
    radii_grid = np.linspace(5, 100, 20)  # microns

    H_r, R = np.meshgrid(heights_r, radii_grid)

    # Calculate slopes for several different velocities
    velocities_to_plot = [5, 10, 15, 20]  # km/s
    colors_v = ['purple', 'blue', 'orange', 'red']

    for velocity, color in zip(velocities_to_plot, colors_v):
        DR_DH = np.zeros_like(H_r)

        for i in range(H_r.shape[0]):
            for j in range(H_r.shape[1]):
                h_km = H_r[i, j]
                r_um = R[i, j]

                if r_um < 1:
                    DR_DH[i, j] = 0
                    continue

                particle = {'r': r_um, 'v': velocity, 'active': True}
                drag = friend.calculate_drag(particle, h_km * 1000)

                # Calculate dr/dh (microns per km)
                dr_dh = friend.calculate_ablation(particle, drag)
                DR_DH[i, j] = dr_dh

        # Plot slope field
        DH_r = np.ones_like(DR_DH)

        # Normalize for visualization
        magnitude_r = np.sqrt(DH_r**2 + DR_DH**2)
        max_mag_r = np.nanmax(np.abs(magnitude_r))
        if max_mag_r > 0:
            scale_factor_r = 2.0 / max_mag_r
            DH_norm_r = DH_r * scale_factor_r
            DR_norm_r = DR_DH * scale_factor_r
        else:
            DH_norm_r = DH_r
            DR_norm_r = DR_DH

        # Plot with reduced density
        skip = 1
        ax2.quiver(H_r[::skip, ::skip], R[::skip, ::skip],
                   -DH_norm_r[::skip, ::skip], DR_norm_r[::skip, ::skip],  # negative because descending
                   color=color, alpha=0.6, scale=20, width=0.003,
                   label=f'v={velocity}km/s')

    ax2.set_xlabel('Height (km)', fontsize=13)
    ax2.set_ylabel('Radius (μm)', fontsize=13)
    ax2.set_title('Slope Field: dr/dh for Different Velocities', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(10, 100)
    ax2.set_ylim(5, 100)

    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f"Slope field saved to {save_path}")
    plt.close()


def overlay_trajectories(save_path="slope_field_with_trajectories.png"):
    """
    Create slope field with actual particle trajectories overlaid
    """

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # ===== LEFT: Velocity slope field with trajectories =====
    heights = np.linspace(100, 10, 15)
    velocities = np.linspace(0.5, 20, 15)
    H, V = np.meshgrid(heights, velocities)

    # Calculate slope field for one representative radius
    radius = 20  # microns
    DV_DH = np.zeros_like(H)

    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            h_km = H[i, j]
            v_kms = V[i, j]

            particle = {'r': radius, 'v': v_kms, 'active': True}
            m = friend.calculate_mass(particle)
            drag = friend.calculate_drag(particle, h_km * 1000)

            v_ms = v_kms * 1000
            if v_ms > 10:
                dv_dh = ((drag / m) - g) / v_ms * 1000
                DV_DH[i, j] = dv_dh
            else:
                DV_DH[i, j] = 0

    # Plot slope field
    DH = np.ones_like(DV_DH)
    magnitude = np.sqrt(DH**2 + DV_DH**2)
    max_mag = np.nanmax(magnitude)
    if max_mag > 0:
        scale = 2.0 / max_mag
        DH_norm = DH * scale
        DV_norm = DV_DH * scale
    else:
        DH_norm = DH
        DV_norm = DV_DH

    ax1.quiver(H, V, -DH_norm, DV_norm, color='gray', alpha=0.4, scale=20, width=0.002)

    # Overlay actual trajectories
    initial_velocities = [5, 10, 15, 20]
    colors = ['blue', 'green', 'orange', 'red']

    for v0, color in zip(initial_velocities, colors):
        h_traj = [100]
        v_traj = [v0]

        h = 100
        v = v0
        r = radius

        while h > 0 and v > 0.01 and r > 0.1:
            particle = {'r': r, 'v': v, 'active': True}
            m = friend.calculate_mass(particle)
            drag = friend.calculate_drag(particle, h * 1000)

            # Update velocity
            v_ms = v * 1000
            if v_ms > 10:
                dv_dh = ((drag / m) - g) / v_ms * 1000
                v = v + dv_dh * (-0.5)
                if v < 0:
                    break

            # Update radius
            dr_dh = friend.calculate_ablation(particle, drag)
            r = r - (dr_dh * (-0.5))

            h = h - 0.5

            h_traj.append(h)
            v_traj.append(v)

        ax1.plot(h_traj, v_traj, color=color, linewidth=2.5, label=f'v₀={v0} km/s', zorder=5)

    ax1.set_xlabel('Height (km)', fontsize=13)
    ax1.set_ylabel('Velocity (km/s)', fontsize=13)
    ax1.set_title(f'Velocity Trajectories (r={radius}μm)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 22)

    # ===== RIGHT: Radius slope field with trajectories =====
    heights_r = np.linspace(100, 10, 15)
    radii_grid = np.linspace(5, 100, 15)
    H_r, R = np.meshgrid(heights_r, radii_grid)

    velocity = 15  # km/s
    DR_DH = np.zeros_like(H_r)

    for i in range(H_r.shape[0]):
        for j in range(H_r.shape[1]):
            h_km = H_r[i, j]
            r_um = R[i, j]

            if r_um < 1:
                continue

            particle = {'r': r_um, 'v': velocity, 'active': True}
            drag = friend.calculate_drag(particle, h_km * 1000)
            dr_dh = friend.calculate_ablation(particle, drag)
            DR_DH[i, j] = dr_dh

    # Plot slope field
    DH_r = np.ones_like(DR_DH)
    magnitude_r = np.sqrt(DH_r**2 + DR_DH**2)
    max_mag_r = np.nanmax(np.abs(magnitude_r))
    if max_mag_r > 0:
        scale_r = 2.0 / max_mag_r
        DH_norm_r = DH_r * scale_r
        DR_norm_r = DR_DH * scale_r
    else:
        DH_norm_r = DH_r
        DR_norm_r = DR_DH

    ax2.quiver(H_r, R, -DH_norm_r, DR_norm_r, color='gray', alpha=0.4, scale=20, width=0.002)

    # Overlay trajectories
    initial_radii = [20, 40, 60, 80]

    for r0, color in zip(initial_radii, colors):
        h_traj = [100]
        r_traj = [r0]

        h = 100
        v = velocity
        r = r0

        while h > 0 and v > 0.01 and r > 0.1:
            particle = {'r': r, 'v': v, 'active': True}
            m = friend.calculate_mass(particle)
            drag = friend.calculate_drag(particle, h * 1000)

            # Update velocity
            v_ms = v * 1000
            if v_ms > 10:
                dv_dh = ((drag / m) - g) / v_ms * 1000
                v = v + dv_dh * (-0.5)
                if v < 0:
                    break

            # Update radius
            dr_dh = friend.calculate_ablation(particle, drag)
            r = r - (dr_dh * (-0.5))

            h = h - 0.5

            h_traj.append(h)
            r_traj.append(r)

        ax2.plot(h_traj, r_traj, color=color, linewidth=2.5, label=f'r₀={r0}μm', zorder=5)

    ax2.set_xlabel('Height (km)', fontsize=13)
    ax2.set_ylabel('Radius (μm)', fontsize=13)
    ax2.set_title(f'Radius Trajectories (v={velocity} km/s)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f"Slope field with trajectories saved to {save_path}")
    plt.close()


if __name__ == "__main__":
    print("Generating slope field...")
    create_slope_field()

    print("Generating slope field with trajectories...")
    overlay_trajectories()

    print("Done!")
