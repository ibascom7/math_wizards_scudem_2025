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
    average_density = 2.0
    volume = (4/3) * np.pi * particle["r"]**3
    mass = average_density * volume
    return mass

"""
Calculate the kinetic energy of the particle
"""
def calculate_ke(particle):
    # KE = (1/2)mv^2
    mass = calculate_mass(particle)
    kinetic_energy = 0.5 * mass * (particle["v"]**2)
    return kinetic_energy

p = {"r": 5, "v": 2}
print(calculate_ke(p))