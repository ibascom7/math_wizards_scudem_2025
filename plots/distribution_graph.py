"""
This will create a graph of the distribution of given particles
by radius and velocity 
"""

import matplotlib as plt

def plot_distribution(particles):
    radii = [particle["r"] for particle in particles]
    velocities = [particle["v"] for particle in particles]
    
    fig, ax = plt.subplots()

    ax.scatter(radii, velocities)

    ax.set_xlabel("Radius (microns)")
    ax.set_ylabel("Velocity (km/s)")
    ax.set_title("Distribution of Particles by Radius and Velocity")

    plt.show()
