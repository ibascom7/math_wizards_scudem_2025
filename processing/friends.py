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


def calculate_air_density(height):
    if height > 25000:
        T = -131.21 + 0.00299*height
        p = 2.488 * ((T + 273.1)/216.6)**(-11.388)
        air_density = p / (0.2869 * (T + 273.1))
    elif height > 11000 & height < 25000:
        T = -56.46
        p = 22.65 * np.exp(1.73 - 0.000157*height)
        air_density = p / (0.2869 * (T + 273.1))
    elif height < 11000:
        T = 15.04 - 0.00649*height
        p = 101.29 * ((T + 273.1)/288.08)**(5.256)
        air_density = p / (0.2869 * (T + 273.1))
    return air_density

def calculate_air_viscocity(height):
    T_0 = 518.7
    air_viscocity_0 = 3.62e-7
    if height > 25000:
        T = -131.21 + 0.00299*height
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))  
    elif height > 11000 & height < 25000:
        T = -56.46
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))

    elif height < 11000:
        T = 15.04 - 0.00649*height
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))
        
    return air_viscocity

def calculate_reynolds(particle, air_density, air_viscocity):
    Re = air_density*particle["v"]*(particle["r"]*2)/air_viscocity

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