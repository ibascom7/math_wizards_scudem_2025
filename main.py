"""
FUN FUN FUN TIME 
FUN FUN RUN TIME
FUN RUN RUN TIME
RUN RUN RUN TIME
"""

import numpy as np
from scipy.constants import g

from gen_particles import generate_particles as gen_p
import processing.friends as friend
from plots.distribution_graph import plot_distribution as plot_d
from plots.distribution_animation import kawaii_anime as sugoi

def go_do_thing_time():
    acc = 0
    history = {}
    particles = gen_p(10000)
    burned_up_count = 0  
    settled_count = 0    
    plot_d(particles, title="Initial Particles", save_path="initial_plot.png")
    step_size = -0.1
    heights = np.arange(100, 0, step_size)

    # Debug: check first particle at several heights
    debug_particle = particles[0]
    debug_heights = [100]  # Just check 100 km in detail
    print("\n=== DETAILED DRAG DIAGNOSTIC ===")
    for h_km in debug_heights:
        h_m = h_km * 1000

        # Get all intermediate values
        air_density = friend.calculate_air_density(h_m)
        air_viscosity = friend.calculate_air_viscosity(h_m)
        Kn = friend.calculate_knudsen(debug_particle, air_density)

        # Get temperature for thermal velocity
        if h_m > 25000:
            T_celsius = -131.21 + (0.00299 * h_m)
        elif h_m >= 11000:
            T_celsius = -56.46
        else:
            T_celsius = 15.04 - (0.00649 * h_m)
        T_kelvin = T_celsius + 273.15

        # Thermal velocity
        k_boltzmann = 1.380649e-23
        m_air = 4.8e-26
        v_thermal = np.sqrt((8 * k_boltzmann * T_kelvin) / (np.pi * m_air))

        # Particle properties
        r_meters = debug_particle["r"] * 1e-6
        v_ms = debug_particle["v"] * 1000

        # Epstein drag
        drag_epstein = (4.0/3.0) * np.pi * (r_meters ** 2) * air_density * v_ms * v_thermal

        drag_force = friend.calculate_drag(debug_particle, h_m)
        mass = friend.calculate_mass(debug_particle)
        accel = drag_force / mass if mass > 0 else 0

        print(f"Height: {h_km} km ({h_m} m)")
        print(f"  Particle: r={debug_particle['r']:.3f} μm, v={debug_particle['v']:.2f} km/s ({v_ms:.0f} m/s)")
        print(f"  Atmosphere: ρ={air_density:.2e} kg/m³, T={T_kelvin:.1f} K ({T_celsius:.1f}°C)")
        print(f"  Knudsen: Kn={Kn:.1f} (regime: {'Free molecular' if Kn > 10 else 'Transition' if Kn > 0.01 else 'Continuum'})")
        print(f"  Thermal velocity: v_th={v_thermal:.1f} m/s")
        print(f"  Epstein drag: F={(4.0/3.0):.3f} × π × r² × ρ × v × v_th")
        print(f"              = {(4.0/3.0):.3f} × {np.pi:.3f} × {(r_meters**2):.2e} × {air_density:.2e} × {v_ms:.0f} × {v_thermal:.1f}")
        print(f"              = {drag_epstein:.2e} N")
        print(f"  Function returned: {drag_force:.2e} N")
        print(f"  Mass: {mass:.2e} kg")
        print(f"  Acceleration: {accel:.2e} m/s² (={accel/9.81:.0f}g)")
    print("======================\n")

    for height in heights:
        for particle in particles:
            # Skip particles that have already burned up or stopped
            if not particle["active"]:
                continue

            m = friend.calculate_mass(particle)
            # Convert height from km to meters for air density calculation
            height_meters = height * 1000
            drag = friend.calculate_drag(particle, height_meters)

            # Convert velocity to m/s for unit consistency
            v_ms = particle["v"] * 1000  # km/s to m/s

            # Check if particle has slowed to settling velocity
            # Cosmic dust slower than ~10 m/s is essentially settling - mark as done
            if particle["v"] < 0.01:  # Less than 10 m/s
                particle["active"] = False
                settled_count += 1
                continue

            # Avoid divide by zero - only update if velocity is reasonable
            if v_ms > 10:  # Only update if velocity > 10 m/s to avoid numerical instability
                # dv/dh in units of 1/s (dimensionless per unit length)
                dv_dh = ((drag / m) - g) / v_ms
                # Multiply by step_size (in km) gives change in velocity (km/s)
                dv_kms = dv_dh * step_size
                particle["v"] = particle["v"] + dv_kms

                # If velocity goes negative, particle has stopped
                if particle["v"] <= 0:
                    particle["active"] = False
                    settled_count += 1
                    continue

            dr_dh = friend.calculate_ablation(particle, drag)
            # dr_dh is negative (ablation), step_size is negative (descending)
            # We want radius to DECREASE, so subtract the product:
            new_radius = particle["r"] - (dr_dh * step_size)

            # Prevent negative radius (numerical issues)
            if new_radius < 0:
                particle["r"] = 0
                particle["active"] = False
                burned_up_count += 1
            elif new_radius < 0.001:
                # Mark as inactive if too small
                particle["r"] = new_radius
                particle["active"] = False
                burned_up_count += 1
            else:
                particle["r"] = new_radius

        # Filter for active particles with valid radii and velocities
        valid_particles = [p for p in particles if p['active'] and p['r'] > 0.01 and p['v'] >= 0]

        history[height] = {
        'radii': [p['r'] for p in valid_particles],
        'velocities': [p['v'] for p in valid_particles],
        'n_particles': len(valid_particles)
        }
        acc += 1

        # Stop simulation if all particles are inactive
        if len(valid_particles) == 0:
            print(f"\nAll particles accounted for at height {height:.1f} km")
            print("Ending simulation early - no active particles remaining")
            break

        # Debug output every 10 km
        if acc % 100 == 0:
            active_particles = [p for p in particles if p['active']]
            if active_particles:
                avg_v = np.mean([p['v'] for p in active_particles])
                avg_r = np.mean([p['r'] for p in active_particles])
                min_r = min([p['r'] for p in active_particles])
                max_r = max([p['r'] for p in active_particles])
                print(f"h={height:.1f} km | n={len(active_particles)} | avg_v={avg_v:.2f} km/s | avg_r={avg_r:.1f} μm (r: {min_r:.1f}-{max_r:.1f})")

                # Detail on one particle
                if len(active_particles) > 0:
                    p = active_particles[0]
                    m = friend.calculate_mass(p)
                    drag = friend.calculate_drag(p, height * 1000)
                    accel = (drag / m) if m > 0 else 0
                    dr_dh = friend.calculate_ablation(p, drag)
                    dr_per_step = dr_dh * step_size
                    print(f"  → Sample: v={p['v']:.3f} km/s, r={p['r']:.1f} μm, drag={drag:.2e} N, a_drag={accel:.2e} m/s²")
                    print(f"     Ablation: dr/dh={dr_dh:.2e} μm/km, dr_this_step={dr_per_step:.4f} μm")
    final_particles = [p for p in particles if p["active"]]

    if len(final_particles) > 0:
        plot_d(final_particles, title="Final Destination", save_path="final_destination")
        print(f"\n{len(final_particles)} particles survived to ground level")
    else:
        print("\nNo particles survived to ground level - all burned up or settled during descent")

    # Print particle fate statistics
    print(f"\n=== PARTICLE FATE STATISTICS ===")
    print(f"Burned up (ablated): {burned_up_count} particles")
    print(f"Settled (slowed down): {settled_count} particles")
    print(f"Still active: {len(final_particles)} particles")
    print(f"Total: {burned_up_count + settled_count + len(final_particles)} particles")
    print(f"================================\n")

    sugoi(history, "beauty_and_grace.gif")
    


go_do_thing_time()