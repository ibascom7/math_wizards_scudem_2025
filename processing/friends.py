"""
This is going to be a library of helper functions
"""

import numpy as np

"""
Calculate the mass of the particle
Assuming that particle density is 2 g/cm^3.
Assuming all particles are spheres.
"""
def calculate_mass(particle):
    # Density = Mass / Volume
    radius_cm = particle["r"] * 1e-4
    average_density = 2.0
    volume = (4/3) * np.pi * radius_cm**3
    mass = average_density * volume
    return mass

"""
Calculate the kinetic energy of the particle
"""
def calculate_reynolds(particle, air_density, air_viscocity):
    Re = air_density*particle["v"]*(particle["r"]*2)/air_viscocity


def calculate_ke(particle):
    # KE = (1/2)mv^2
    mass = calculate_mass(particle)
    kinetic_energy = 0.5 * mass * (particle["v"]**2)
    return kinetic_energy


p = {"r": 5, "v": 2}


print(calculate_ke(p))