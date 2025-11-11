"""
This is going to be a library of helper functions
"""

import numpy as np



"""
Calculate the mass of the particle
Assuming that particle density is 2 g/cm^3.
Assuming all particles are spheres.

Returns:
    mass in kilograms (kg) for SI unit consistency
"""
def calculate_mass(particle):
    # Density = Mass / Volume
    radius_cm = particle["r"] * 1e-4  # microns to cm
    average_density = 2.0  # g/cm³
    volume = (4.0/3.0) * np.pi * (radius_cm**3)  # cm³
    mass_grams = average_density * volume  # grams
    mass_kg = mass_grams / 1000.0  # Convert to kilograms for SI consistency
    return mass_kg



"""
Using the Ideal Gas Law
Calculate the air density at a given height. Temperature varies by height.
Temperatures and pressures of Upper Stratosphere, Lower Stratosphere, and Troposphere are used
"""
def calculate_air_density(height):
    if height > 25000:
        T = -131.21 + (0.00299 * height)
        p = 2.488 * (((T + 273.1) / 216.6) ** (-11.388))
        air_density = p / (0.2869 * (T + 273.1))
    elif height >= 11000 and height <= 25000:
        T = -56.46
        p = 22.65 * np.exp(1.73 - (0.000157 * height))
        air_density = p / (0.2869 * (T + 273.1))
    else:  # height < 11000
        T = 15.04 - (0.00649 * height)
        p = 101.29 * (((T + 273.1) / 288.08) ** 5.256)
        air_density = p / (0.2869 * (T + 273.1))
    return air_density





"""
Calculate the air viscocity at a given height
Viscocity is dynamic through atmosphere. Viscocity coefficient is constant with pressure but depends on the temperature of air which varies by 
layer of the atmosphere.
"""
def calculate_air_viscosity(height):
    T_0 = 518.7
    air_viscosity_0 = 3.62e-7
    if height > 25000:
        T = -131.21 + (0.00299 * height)
        T_abs = T + 459.67  # Convert Fahrenheit to Rankine (absolute scale)
        air_viscosity = air_viscosity_0 * ((T_abs / T_0) ** 1.5) * ((T_0 + 198.72) / (T_abs + 198.72))
    elif height >= 11000 and height <= 25000:
        T = -56.46
        T_abs = T + 459.67  # Convert Fahrenheit to Rankine
        air_viscosity = air_viscosity_0 * ((T_abs / T_0) ** 1.5) * ((T_0 + 198.72) / (T_abs + 198.72))
    else:  # height < 11000
        T = 15.04 - (0.00649 * height)
        T_abs = T + 459.67  # Convert Fahrenheit to Rankine
        air_viscosity = air_viscosity_0 * ((T_abs / T_0) ** 1.5) * ((T_0 + 198.72) / (T_abs + 198.72))

    return air_viscosity




"""
Calculate the Reynolds number of the particle

Reynolds number: Re = ρ * v * L / μ
where L is the characteristic length (diameter for spheres)

Dimensionless number indicating flow regime:
- Re < 1: Viscous/Stokes flow
- 1 < Re < 1000: Intermediate
- Re > 1000: Turbulent flow

Parameters:
    particle: dict with 'r' (radius in microns) and 'v' (velocity in km/s)
    air_density: kg/m³
    air_viscosity: Pa·s (kg/(m·s))

Returns:
    Reynolds number (dimensionless)
"""
def calculate_reynolds(particle, air_density, air_viscosity):
    # Convert to SI units
    v_ms = particle["v"] * 1000        # km/s to m/s
    diameter_m = particle["r"] * 2e-6  # radius in microns to diameter in meters

    Re = (air_density * v_ms * diameter_m) / air_viscosity
    return Re




"""
Calculate the Knudsen number of the particle

Knudsen number: Kn = λ / L
where λ is the mean free path of gas molecules
and L is the characteristic length (radius for spheres)

Dimensionless number indicating gas rarefaction:
- Kn < 0.01: Continuum flow (normal fluid dynamics)
- 0.01 < Kn < 0.1: Slip flow (small correction needed)
- 0.1 < Kn < 10: Transition regime
- Kn > 10: Free molecular flow

Parameters:
    particle: dict with 'r' (radius in microns)

Returns:
    Knudsen number (dimensionless)
"""
def calculate_knudsen(particle):
    # Mean free path of air molecules at standard conditions
    # λ ≈ 65 nm = 0.065 microns
    # (varies with altitude, but this is a reasonable approximation)
    lambda_microns = 0.065  # Mean free path in microns
    r_microns = particle["r"]

    Kn = lambda_microns / r_microns
    return Kn





"""
Calculate drag force on a particle

Uses different drag models depending on Reynolds number (Re):
- Re < 1: Stokes regime (linear drag)
  - With Knudsen correction for rarefied gas (high altitude)
- 1 < Re < 1000: Intermediate regime (empirical correlation)
- Re > 1000: High Reynolds regime (quadratic drag)

Parameters:
    particle: dict with 'r' (radius in microns) and 'v' (velocity in km/s)
    height: altitude in meters

Returns:
    drag force in Newtons
"""
def calculate_drag(particle, height):

    air_viscosity = calculate_air_viscosity(height)
    air_density = calculate_air_density(height)

    # Convert particle properties to SI units
    r_meters = particle["r"] * 1e-6  # microns to meters
    v_ms = particle["v"] * 1000       # km/s to m/s

    # Avoid issues with zero/tiny values
    if v_ms < 0.001 or r_meters < 1e-10:
        return 0

    # Determine regime using Reynolds number
    Re = calculate_reynolds(particle, air_density, air_viscosity)

    if Re < 1:
        # ===== STOKES REGIME (low Re) =====
        # Basic formula: F = 6πμrv

        Kn = calculate_knudsen(particle)

        if Kn < 0.1:
            # Continuum flow - standard Stokes drag
            drag = 6 * np.pi * air_viscosity * r_meters * v_ms
        else:
            # Rarefied flow - apply Cunningham slip correction
            # Common at high altitudes where air is thin
            alpha = 1.207
            beta = 0.376
            gamma = 0.332
            C_slip = 1 + (Kn * (alpha + (beta * np.exp(-gamma / Kn))))
            drag = (6 * np.pi * air_viscosity * r_meters * v_ms) / C_slip

    elif 1 <= Re < 1000:
        # ===== INTERMEDIATE REGIME =====
        # Empirical correlation for transition region
        C_d = (24.0 / Re) + (4.0 / np.sqrt(Re)) + 0.4
        A = np.pi * (r_meters ** 2)  # Cross-sectional area
        drag = 0.5 * C_d * air_density * A * (v_ms ** 2)

    else:
        # ===== HIGH REYNOLDS REGIME (Re >= 1000) =====
        # Turbulent flow - constant drag coefficient
        C_d = 0.44  # Typical for sphere in turbulent flow
        A = np.pi * (r_meters ** 2)
        drag = 0.5 * C_d * air_density * A * (v_ms ** 2)

    return drag



"""
Calculate the kinetic energy of the particle
"""
def calculate_ke(particle):
    # KE = (1/2)mv^2
    mass = calculate_mass(particle)
    kinetic_energy = 0.5 * mass * (particle["v"]**2)
    return kinetic_energy



"""
Calculate the ablation of the particle
We are defining this as the loss of radius due to heat from atmospheric friction.

This function returns dr/dh (change in radius per unit height)

Physical model:
- Kinetic energy from drag converts to heat
- Heat causes evaporation/melting of particle material
- Only a fraction (epsilon) of drag energy goes into heating the particle
- Rest is radiated away or heats surrounding air

Assumptions:
- Uniform particle density (2 g/cm³ - typical for cosmic dust)
- Ignore chemical composition (treat as generic silicate)
- Constant heat of vaporization

Parameters:
    particle: dict with 'r' (radius in microns) and 'v' (velocity in km/s)
    F_drag: drag force in Newtons

Returns:
    dr/dh: rate of radius change in microns/km (negative value = shrinking)
"""
def calculate_ablation(particle, F_drag):
    # Physical constants
    epsilon = 0.3  # Heat transfer efficiency (fraction of drag energy that heats particle)
    L_vap = 6e6  # J/kg - Heat of vaporization for silicates (~ 6000 kJ/kg)
    rho_particle = 2000  # kg/m³ - Particle density (2 g/cm³)

    # Convert units
    r_meters = particle['r'] * 1e-6  # microns to meters
    v_ms = particle['v'] * 1000  # km/s to m/s

    # Avoid division by zero
    if particle['v'] < 0.00001 or r_meters < 1e-10:  # Very slow or tiny particles
        return 0

    # Energy going into heating per unit height (J/m)
    # heating_rate = epsilon * F_drag * v / v_height
    # Since dh/dt = -v, we have: heating per height = epsilon * F_drag * v / v = epsilon * F_drag
    heating_per_height = epsilon * F_drag  # J/m (energy per meter of altitude change)

    # Mass loss per unit height (kg/m)
    dm_dh = (-heating_per_height) / L_vap

    # Convert mass loss to radius change
    # m = (4/3)πr³ρ → dm = 4πr²ρ dr
    # dr/dh = (dm/dh) / (4πr²ρ)
    if r_meters > 0:
        dr_dh = dm_dh / (4 * np.pi * (r_meters ** 2) * rho_particle)
    else:
        dr_dh = 0

    # Convert back to microns/km
    # dr_dh is in m/m (dimensionless)
    # To get microns/km: multiply by 1e6 (m→microns) and divide by 1000 (per-m → per-km)
    dr_dh_microns_per_km = dr_dh * 1e6 / 1000  # m/m → microns/km = dr_dh * 1000

    return dr_dh_microns_per_km  # Will be negative (radius decreases)




"""
Check if particle fragments due to ram pressure

Ram pressure is the dynamic pressure from moving through air:
    P_ram = (1/2) * ρ_air * v²

When this pressure exceeds the material strength, the particle shatters.

Physical model:
- Cosmic dust particles have tensile strength ~ 10⁵ to 10⁶ Pa (varies by composition)
- When ram pressure > strength, particle fragments
- Fragments are typically 2-5 smaller pieces
- Each fragment has roughly equal mass (on average)

Parameters:
    particle: dict with 'r' (radius in microns), 'v' (velocity in km/s)
    air_density: kg/m³

Returns:
    list of new fragment particles if fragmentation occurs, empty list otherwise
    Original particle should be marked inactive if it fragments
"""
def check_fragmentation(particle, air_density):
    # Material properties
    strength = 1e5  # Pa - Tensile strength of cosmic dust (lower bound)
                    # Range: 10⁵ to 10⁶ Pa depending on composition

    # Convert units
    v_ms = particle['v'] * 1000  # km/s to m/s

    # Calculate ram pressure (Pa)
    ram_pressure = 0.5 * air_density * (v_ms ** 2)

    # Check if fragmentation occurs
    if ram_pressure > strength:
        # Particle fragments!
        # Decide number of fragments (2-4 typically)
        n_fragments = np.random.randint(2, 5)

        # Conservation of mass: each fragment has mass/n_fragments
        # For spheres: if m_fragment = m_parent/n, then r_fragment = r_parent / n^(1/3)
        fragment_radius = particle['r'] / (n_fragments ** (1/3))

        # Fragments have similar velocity with some spread
        fragments = []
        for i in range(n_fragments):
            # Velocity varies slightly (±10%)
            velocity_variation = np.random.uniform(0.9, 1.1)

            fragment = {
                'r': fragment_radius,
                'v': particle['v'] * velocity_variation,
                'active': True
            }
            fragments.append(fragment)

        return fragments
    else:
        # No fragmentation
        return []