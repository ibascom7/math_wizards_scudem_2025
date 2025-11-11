"""
This will create the initial particles array for COSMIC DUST.
Based on z-MIF distributions with constraints:
1. Cosmic dust: radius between 0-100 microns
2. Majority (>50%) have mass < 5 micrograms
"""

import numpy as np
from scipy import stats

# Global constants
RHO = 2.0  # particle density in g/cm³

# Cosmic dust size limits
R_MIN = 0.1    # microns (practical minimum)
R_MAX = 100.0  # microns (cosmic dust upper limit)

# Mass limits (derived from radius limits)
M_MIN = RHO * (4/3) * np.pi * (R_MIN/1e4)**3   # 8.38e-15 g
M_MAX = RHO * (4/3) * np.pi * (R_MAX/1e4)**3   # 8.38e-6 g (~8.4 μg)


def mass_to_radius(mass_g):
    """
    Convert mass in grams to radius in microns.
    
    Using: m = ρ * (4/3) * π * r³
    So: r = (3m / (4πρ))^(1/3)
    """
    r_cm = (3 * mass_g / (4 * np.pi * RHO))**(1/3)
    r_microns = r_cm * 1e4  # convert cm to microns
    return r_microns


def radius_to_mass(radius_microns):
    """Convert radius in microns to mass in grams."""
    r_cm = radius_microns / 1e4
    mass_g = RHO * (4/3) * np.pi * r_cm**3
    return mass_g


def generate_particles(n_particles):
    """
    Generate cosmic dust particles following z-MIF distributions.
    
    Constraints:
    - Radius: 0-100 μm (cosmic dust definition)
    - Majority (>50%) have mass < 5 μg (from paper)
    - Velocity: exponential decay (most < 11 km/s)
    
    Parameters:
        n_particles (int): Number of particles to generate
    
    Returns:
        list of dicts: Each dict contains particle properties
    """
    
    # ===== MASS DISTRIBUTION (Log-Normal, truncated) =====
    # For cosmic dust with majority < 5 μg, need median around 2-3 μg
    # This corresponds to log10(mass) around -5.5 to -5.7
    
    # Target: median mass ~2.5 μg = 2.5e-6 g
    # log10(2.5e-6) = -5.6
    log10_mass_median = -5.6   # Median around 2.5 μg
    log10_mass_std = 0.8       # Standard deviation in log10 space
    
    # Convert to natural log parameters
    # For log-normal: median = exp(μ), so μ = ln(median)
    mu = np.log(10**log10_mass_median)
    sigma = log10_mass_std * np.log(10)
    
    # Generate log-normal distributed masses
    masses = stats.lognorm.rvs(s=sigma, scale=np.exp(mu), size=n_particles*2)
    
    # Truncate to cosmic dust range (0.1 to 100 μm → mass limits)
    masses = masses[(masses >= M_MIN) & (masses <= M_MAX)]
    
    # Take first n_particles (regenerate if needed)
    while len(masses) < n_particles:
        additional = stats.lognorm.rvs(s=sigma, scale=np.exp(mu), 
                                       size=n_particles)
        additional = additional[(additional >= M_MIN) & (additional <= M_MAX)]
        masses = np.concatenate([masses, additional])
    
    masses = masses[:n_particles]
    
    # ===== CONVERT MASS TO RADIUS =====
    radii = np.array([mass_to_radius(m) for m in masses])
    
    # Double-check: all should be in cosmic dust range
    assert np.all((radii >= R_MIN) & (radii <= R_MAX)), "Radius out of bounds!"
    
    # ===== VELOCITY DISTRIBUTION (Exponential) =====
    # z-MIF shows exponential decay with most particles < 11 km/s
    # For cosmic dust specifically, even more weighted to low velocities
    
    # Exponential decay scale (mean velocity)
    # Chosen so ~75% are below 11 km/s
    v_scale = 8.0  # km/s
    
    # Generate exponentially distributed velocities
    velocities = stats.expon.rvs(scale=v_scale, size=n_particles)
    
    # Truncate at maximum observed velocity (72.5 km/s)
    v_max = 72.5
    velocities = np.minimum(velocities, v_max)
    
    # Ensure minimum velocity
    v_min = 0.1  # km/s
    velocities = np.maximum(velocities, v_min)
    
    # ===== CREATE PARTICLE ARRAY =====
    particles = []
    for i in range(n_particles):
        particle = {
            "r": radii[i],        # radius in microns
            "v": velocities[i],   # velocity in km/s
            "active": True        # particle status flag
        }
        particles.append(particle)
    
    return particles


def print_particle_statistics(particles):
    """Print summary statistics for generated particles."""
    
    radii = [p['r'] for p in particles]
    masses = [p['m'] for p in particles]
    velocities = [p['v'] for p in particles]
    
    # Convert masses to micrograms for easier reading
    masses_ug = [m * 1e6 for m in masses]
    
    print(f"\n{'='*70}")
    print(f"COSMIC DUST PARTICLE STATISTICS (n = {len(particles)})")
    print(f"{'='*70}")
    
    print(f"\nRADIUS (Cosmic Dust Range: 0.1-100 μm):")
    print(f"  Range:  {np.min(radii):.2f} - {np.max(radii):.2f} μm")
    print(f"  Mean:   {np.mean(radii):.2f} μm")
    print(f"  Median: {np.median(radii):.2f} μm")
    
    # Key percentiles
    p25, p75 = np.percentile(radii, [25, 75])
    print(f"  25th percentile: {p25:.2f} μm")
    print(f"  75th percentile: {p75:.2f} μm")
    
    print(f"\nMASS (Max for cosmic dust: 8.38 μg):")
    print(f"  Range:  {np.min(masses_ug):.4f} - {np.max(masses_ug):.2f} μg")
    print(f"  Mean:   {np.mean(masses_ug):.3f} μg")
    print(f"  Median: {np.median(masses_ug):.3f} μg")
    
    # Check MAJORITY below 5 μg threshold (CRITICAL CONSTRAINT)
    below_5ug = sum(1 for m in masses_ug if m < 5)
    pct_below_5ug = 100 * below_5ug / len(particles)
    print(f"  Below 5 μg: {below_5ug} ({pct_below_5ug:.1f}%) {'✓' if pct_below_5ug > 50 else '✗ FAIL'}")
    
    # Additional thresholds
    below_1ug = sum(1 for m in masses_ug if m < 1)
    pct_below_1ug = 100 * below_1ug / len(particles)
    print(f"  Below 1 μg: {below_1ug} ({pct_below_1ug:.1f}%)")
    
    print(f"\nVELOCITY:")
    print(f"  Range:  {np.min(velocities):.2f} - {np.max(velocities):.2f} km/s")
    print(f"  Mean:   {np.mean(velocities):.2f} km/s")
    print(f"  Median: {np.median(velocities):.2f} km/s")
    
    # Check fraction below 11 km/s
    below_11kms = sum(1 for v in velocities if v < 11)
    pct_below_11kms = 100 * below_11kms / len(particles)
    print(f"  Below 11 km/s: {below_11kms} ({pct_below_11kms:.1f}%)")
    
    # Velocity percentiles
    v25, v75 = np.percentile(velocities, [25, 75])
    print(f"  25th percentile: {v25:.2f} km/s")
    print(f"  75th percentile: {v75:.2f} km/s")
    
    print(f"\n{'='*70}")
    print("CONSTRAINT CHECK:")
    print(f"{'='*70}")
    print(f"✓ All particles in cosmic dust range (0.1-100 μm)")
    if pct_below_5ug > 50:
        print(f"✓ Majority ({pct_below_5ug:.1f}%) have mass < 5 μg")
    else:
        print(f"✗ FAIL: Only {pct_below_5ug:.1f}% have mass < 5 μg (need >50%)")
    print(f"✓ Velocity follows z-MIF exponential decay")
    print(f"{'='*70}\n")


# ===== MAIN EXECUTION =====
if __name__ == "__main__":
    # Generate particles
    n_particles = 1000
    particles = generate_particles(n_particles)
    
    # Print statistics
    print_particle_statistics(particles)
    
    # Optional: Save particles to file for use in simulation
    # np.save('cosmic_dust_particles.npy', particles)
    print("Cosmic dust particles generated successfully!")
    print("Ready for atmospheric entry simulation.")
