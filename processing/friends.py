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
def calculate_ke(particle):
    # KE = (1/2)mv^2
    mass = calculate_mass(particle)
    kinetic_energy = 0.5 * mass * (particle["v"]**2)
    return kinetic_energy

def calculate_air_density(height):
    # Equation for air density given a height
    air_density = 100
    # Air density should increase as height decreases 
    return air_density

def drag_coefficent(particle, air_density):
    # Equation or method of finding drag coefficient based on air density
    D_c = .1
    return D_c

def drag_force(particle, drag_coef, air_density):
    # Cases for different sizes
    # Equations for drag of each case
    drag = 1000000
    return drag

p = {"r": 5, "v": 2}


print(calculate_ke(p))