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
Using the Ideal Gas Law
Calculate the air density at a given height. Temperature varies by height.
Temperatures and pressures of Upper Stratosphere, Lower Stratosphere, and Troposphere are used
"""
def calculate_air_density(height):
    if height > 25000:
        T = -131.21 + 0.00299*height
        p = 2.488 * ((T + 273.1)/216.6)**(-11.388)
        air_density = p / (0.2869 * (T + 273.1))
    elif height > 11000 and height < 25000:
        T = -56.46
        p = 22.65 * np.exp(1.73 - 0.000157*height)
        air_density = p / (0.2869 * (T + 273.1))
    elif height < 11000:
        T = 15.04 - 0.00649*height
        p = 101.29 * ((T + 273.1)/288.08)**(5.256)
        air_density = p / (0.2869 * (T + 273.1))
    return air_density
"""
Calculate the air viscocity at a given height
Viscocity is dynamic through atmosphere. Viscocity coefficient is constant with pressure but depends on the temperature of air which varies by 
layer of the atmosphere.
"""
def calculate_air_viscocity(height):
    T_0 = 518.7
    air_viscocity_0 = 3.62e-7
    if height > 25000:
        T = -131.21 + 0.00299*height
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))  
    elif height > 11000 and height < 25000:
        T = -56.46
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))

    elif height < 11000:
        T = 15.04 - 0.00649*height
        air_viscocity = air_viscocity_0*((T/T_0)**(1.5))*((T_0 + 198.72)/(T + 198.72))

    return air_viscocity


""" 
Calculate the Reynolds number of the particle
Spherical assumption comes back, we compute Reynolds number based on radius
"""
def calculate_reynolds(particle, air_density, air_viscocity):
    Re = air_density*particle["v"]*(particle["r"]*2)/air_viscocity
    return Re
"""
Calculate the Knudsen number of the particle
Spherical assumption comes back, we compute Reynolds number based on radius
"""
def calculate_knudsen(particle):
    lbda = 65e-9
    Kn = lbda / particle["r"]
    return Kn

"""
We assume the fluid is a continuous medium, allows us to use the classical fluid dynamics equations with Stokes' law

Calculate drag of a particle
"""
def calculate_drag(particle, air_viscocity, air_density):
    """
    Determine regime
    """
    Re = calculate_reynolds(particle, air_viscocity, air_density)
    if Re < 1:
        """Stokes Regime"""
        Kn = calculate_knudsen(particle)
        if Kn < 0.1:
            drag = 6*np.pi*air_viscocity*particle["r"]*particle["v"]
        elif Kn >= 0.1:
            alpha = 1.207
            beta = 0.376
            gamma  = 0.332
            C_slip = 1 + Kn*(alpha + beta*np.exp(-gamma/Kn))
            drag = C_slip*6*np.pi*air_viscocity*particle["r"]*particle["v"]
    elif Re > 1 & Re < 100:
        """Intermediate Regime"""
        C_d = 24/Re + 4/np.sqrt(Re) + 0.4
        drag = (0.5)*(C_d)*(air_density)*(np.pi*(particle["r"]**2)*(particle["v"]**2))


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