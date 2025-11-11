"""
This will create a graph of the distribution of given particles
by radius and velocity

Parameters:
    particles: list of dicts with 'r' (radius in microns) and 'v' (velocity in km/s)
    title: optional custom title for the plot
    save_path: optional path to save the plot (e.g., 'initial_dist.png')
               If not provided, will display the plot instead
"""

import matplotlib.pyplot as plt

def plot_distribution(particles, title="Distribution of Particles by Radius and Velocity", save_path=None):
    radii = [particle["r"] for particle in particles]
    velocities = [particle["v"] for particle in particles]

    fig, ax = plt.subplots()

    ax.scatter(radii, velocities, alpha=0.6, s=10)  # Added alpha for transparency, smaller dots

    ax.set_xlabel("Radius (microns)")
    ax.set_ylabel("Velocity (km/s)")
    ax.set_title(title)
    ax.set_xlim(0, 100) 
    ax.set_ylim(0, 72.5)

    # Add grid for better readability
    ax.grid(True, alpha=0.3)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
        plt.close()  # Close figure to free memory
    else:
        plt.show()
