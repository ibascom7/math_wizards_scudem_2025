"""
FUN FUN FUN TIME 
FUN FUN RUN TIME
FUN RUN RUN TIME
RUN RUN RUN TIME
"""

import numpy as np
from scipy.constants import g

from processing.init_particles import generate_particles as gen_p
import processing.friends as friend
from plots.distribution_graph import plot_distribution as plot_d
from plots.distribution_animation import kawaii_anime as sugoi

def go_do_thing_time():
    acc = 0
    history = {}
    particles = gen_p(100) 
    plot_d(particles, title="Initial Particles", save_path="initial_plot.png")
    step_size = -0.01
    heights = np.arange(100, 0, step_size)

    # Debug: check first particle at several heights
    debug_particle = particles[0]
    debug_heights = [100, 80, 60, 40, 20, 10]
    print("\n=== DRAG DIAGNOSTIC ===")
    for h_km in debug_heights:
        h_m = h_km * 1000
        drag_force = friend.calculate_drag(debug_particle, h_m)
        mass = friend.calculate_mass(debug_particle)
        accel = drag_force / mass if mass > 0 else 0
        print(f"Height: {h_km} km | v={debug_particle['v']:.2f} km/s | r={debug_particle['r']:.3f} μm | Drag={drag_force:.2e} N | a={accel:.2e} m/s²")
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

            # Avoid divide by zero - skip velocity update if already at minimum
            if v_ms > 0.1:  # Only update if velocity > 0.1 m/s
                # dv/dh in units of 1/s (dimensionless per unit length)
                dv_dh = ((drag / m) - g) / v_ms
                # Multiply by step_size (in km) gives change in velocity (km/s)
                # because: (1/s) * (km) = (1/s) * (1000 m) = (1000 m/s) = (1 km/s) works out!
                dv_kms = dv_dh * step_size
                particle["v"] = particle["v"] + dv_kms

                # If velocity drops too low or goes negative, clamp to minimum
                if particle["v"] <= 0:
                    particle["v"] = 0

            dr_dh = friend.calculate_ablation(particle, drag)
            # dr_dh is negative (ablation), step_size is negative (descending)
            # We want radius to DECREASE, so subtract the product:
            particle["r"] = particle["r"] - (dr_dh * step_size)
            if particle["r"] < 0.001:
                particle["active"] = False

        history[height] = {
        'radii': [p['r'] for p in particles if p['active']],
        'velocities': [p['v'] for p in particles if p['active']],
        'n_particles': len([p for p in particles if p['active']])
        }
        acc += 1

        # Debug output every 10 km
        if acc % 1000 == 0:
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
                    print(f"  → Sample: v={p['v']:.3f} km/s, r={p['r']:.1f} μm, drag={drag:.2e} N, a_drag={accel:.2e} m/s², a_grav={g:.2f} m/s²")
    final_particles = [p for p in particles if p["active"]]
    plot_d(final_particles, title="Final Destination", save_path="final_destination")
    sugoi(history, "beauty_and_grace.gif")
    


go_do_thing_time()