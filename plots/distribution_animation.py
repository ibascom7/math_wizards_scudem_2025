"""
Take the history list of the sim run
and make a graph of the distribution at each height
and then make a final animation of all the graphs.
Attempting to convey fragmentation and ablation.
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


"""
Sugoi Sugoi
Assumes history is a dictionary with heights as keys
Args:
    history: dict with heights as keys, each containing 'radii', 'velocities', 'n_particles'
    save_path: optional path to save animation (e.g., 'animation.gif' or 'animation.mp4')
"""
def kawaii_anime(history, save_path=None):
    fig, ax = plt.subplots()
    
    scatter = ax.scatter([], [])

    ax.set_xlabel("Radius (microns)")
    ax.set_ylabel("Velocity (km/s)")

    # Gives the heights as a list from 100 to 0
    heights = sorted(history.keys(), reverse=True)

    # Calculate axis limits from all data
    all_radii = []
    all_velocities = []
    for h in history.values():
        all_radii.extend(h['radii'])
        all_velocities.extend(h['velocities'])

    ax.set_xlim(0, max(all_radii) * 1.1)  # Add 10% padding
    ax.set_ylim(0, max(all_velocities) * 1.1)

    def update(frame):
          # frame is passed as 0, 1, 2, ... by FuncAnimation()
          height = heights[frame]

          radii = history[height]["radii"]
          velocities = history[height]["velocities"]

          # Create properly shaped array for scatter plot
          if radii:
              points = np.column_stack([radii, velocities])
          else:
              points = np.empty((0, 2))  # Empty array with shape (0, 2)

          scatter.set_offsets(points)

          ax.set_title(f"Height: {height:.1f} km")

          return scatter,
    
    anim = animation.FuncAnimation(
         fig,
         update,
         frames = len(heights),
         interval = 100,
         blit = True
    )

    if save_path:
        print(f"Saving animation to {save_path}...")
        if save_path.endswith('.gif'):
            anim.save(save_path, writer='pillow', fps=10)
        elif save_path.endswith('.mp4'):
            anim.save(save_path, writer='ffmpeg', fps=10)
        else:
            print(f"Warning: Unknown file format. Supported: .gif, .mp4")
            anim.save(save_path, writer='pillow', fps=10)
        print("Animation saved!")
    else:
        plt.show()

    return anim