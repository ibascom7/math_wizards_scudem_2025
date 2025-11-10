"""
This will create the initial particles array.
Based on whatever initial conditions are given.
"""

import numpy as np
from scipy import stats

def generate_particles(n_particles):
    """
    Generate cosmic dust particles with realistic distributions.

    Radius: Power law distribution (exponent -3.5)
    Velocity: Gamma distribution (peaked at ~14 km/s)

    Returns:
        list of dicts: Each dict contains particle properties
    """

    # ===== RADIUS DISTRIBUTION (Power Law) =====
    # Power law: P(r) ~ r^(-3.5)
    # Range: 0.1 to 100 microns (avoiding exactly 0)
    r_min = 0.1  # microns
    r_max = 100  # microns
    alpha = 3.5  # power law exponent

    # Generate power law distributed radii
    # Using inverse transform sampling: r = (r_min^(1-alpha) + u*(r_max^(1-alpha) - r_min^(1-alpha)))^(1/(1-alpha))
    u = np.random.uniform(0, 1, n_particles)
    radii = (r_min**(1-alpha) + u * (r_max**(1-alpha) - r_min**(1-alpha)))**(1/(1-alpha))

    # ===== VELOCITY DISTRIBUTION (Gamma) =====
    # Gamma distribution parameters chosen to match z-MIF characteristics:
    # - Mode around 14 km/s
    # - Mean around 18 km/s
    # - Range: 0 to 72.5 km/s
    shape = 4.5   # k parameter (mode = (k-1)*scale for k >= 1)
    scale = 4.0   # theta parameter

    # Generate gamma distributed velocities
    velocities = stats.gamma.rvs(shape, scale=scale, size=n_particles)

    # Truncate at maximum velocity (72.5 km/s)
    v_max = 72.5
    velocities = np.minimum(velocities, v_max)

    # ===== CREATE PARTICLE ARRAY =====
    particles = []
    for i in range(n_particles):
        particle = {
            "r": radii[i],      # radius in microns
            "v": velocities[i],  # velocity in km/s
            "active": True
        }
        particles.append(particle)

    return particles

# Generate particles
n_particles = 1000
particles = generate_particles(n_particles)

print(f"Generated {n_particles} particles")
print(f"Radius range: {min(p['r'] for p in particles):.2f} - {max(p['r'] for p in particles):.2f} microns")
print(f"Velocity range: {min(p['v'] for p in particles):.2f} - {max(p['v'] for p in particles):.2f} km/s")
print(f"Mean radius: {np.mean([p['r'] for p in particles]):.2f} microns")
print(f"Mean velocity: {np.mean([p['v'] for p in particles]):.2f} km/s")
