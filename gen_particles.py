"""
Generate cosmic dust particles using z-MIF distributions.
Based on research paper constraints for cosmic dust (0.1-100 microns).
"""

import numpy as np
from scipy import stats


def generate_particles(n_particles):
    """
    Generate cosmic dust particles following z-MIF distributions.

    Radius: Log-normal distribution (cosmic dust range: 0.1-100 μm)
    Velocity: Exponential distribution (most < 11 km/s)

    Parameters:
        n_particles (int): Number of particles to generate

    Returns:
        list of dicts: Each dict contains particle properties (r, v, active)
    """

    # ===== RADIUS DISTRIBUTION (Log-Normal) =====
    # Cosmic dust range: 0.1 to 100 microns
    R_MIN = 0.1    # microns
    R_MAX = 100.0  # microns

    # Log-normal parameters
    # Median around 2-3 microns, with spread to cover 0.1-100 range
    log10_radius_median = 0.4   # log10(2.5 μm) ≈ 0.4
    log10_radius_std = 0.9      # Standard deviation in log10 space

    # Convert to natural log parameters for scipy
    mu = np.log(10**log10_radius_median)
    sigma = log10_radius_std * np.log(10)

    # Generate log-normal distributed radii
    radii = stats.lognorm.rvs(s=sigma, scale=np.exp(mu), size=n_particles*2)

    # Truncate to cosmic dust range
    radii = radii[(radii >= R_MIN) & (radii <= R_MAX)]

    # Regenerate if needed to get enough particles
    while len(radii) < n_particles:
        additional = stats.lognorm.rvs(s=sigma, scale=np.exp(mu), size=n_particles)
        additional = additional[(additional >= R_MIN) & (additional <= R_MAX)]
        radii = np.concatenate([radii, additional])

    radii = radii[:n_particles]

    # ===== VELOCITY DISTRIBUTION (Exponential) =====
    # z-MIF shows exponential decay with most particles < 11 km/s
    v_scale = 8.0  # km/s (mean velocity, chosen so ~75% are below 11 km/s)

    # Generate exponentially distributed velocities
    velocities = stats.expon.rvs(scale=v_scale, size=n_particles)

    # Truncate at observed limits
    v_min = 0.1   # km/s (minimum)
    v_max = 72.5  # km/s (maximum observed)
    velocities = np.clip(velocities, v_min, v_max)

    # ===== CREATE PARTICLE ARRAY =====
    particles = []
    for i in range(n_particles):
        particle = {
            "r": radii[i],        # radius in microns
            "v": velocities[i],   # velocity in km/s
            "active": True
        }
        particles.append(particle)

    print(f"\nGenerated {n_particles} particles")
    print(f"\nRADIUS (Cosmic Dust: 0.1-100 μm):")
    print(f"  Range:  {np.min(radii):.2f} - {np.max(radii):.2f} μm")
    print(f"  Mean:   {np.mean(radii):.2f} μm")
    print(f"  Median: {np.median(radii):.2f} μm")

    print(f"\nVELOCITY:")
    print(f"  Range:  {np.min(velocities):.2f} - {np.max(velocities):.2f} km/s")
    print(f"  Mean:   {np.mean(velocities):.2f} km/s")
    print(f"  Median: {np.median(velocities):.2f} km/s")

    below_11kms = sum(1 for v in velocities if v < 11)
    pct_below_11kms = 100 * below_11kms / len(particles)
    print(f"  Below 11 km/s: {below_11kms} ({pct_below_11kms:.1f}%)")

    return particles


# Test/demo code
if __name__ == "__main__":
    n_particles = 1000
    particles = generate_particles(n_particles)

    radii = [p['r'] for p in particles]
    velocities = [p['v'] for p in particles]

    print(f"\nGenerated {n_particles} particles")
    print(f"\nRADIUS (Cosmic Dust: 0.1-100 μm):")
    print(f"  Range:  {np.min(radii):.2f} - {np.max(radii):.2f} μm")
    print(f"  Mean:   {np.mean(radii):.2f} μm")
    print(f"  Median: {np.median(radii):.2f} μm")

    print(f"\nVELOCITY:")
    print(f"  Range:  {np.min(velocities):.2f} - {np.max(velocities):.2f} km/s")
    print(f"  Mean:   {np.mean(velocities):.2f} km/s")
    print(f"  Median: {np.median(velocities):.2f} km/s")

    below_11kms = sum(1 for v in velocities if v < 11)
    pct_below_11kms = 100 * below_11kms / len(particles)
    print(f"  Below 11 km/s: {below_11kms} ({pct_below_11kms:.1f}%)")
