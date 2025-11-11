"""
Starting at the given height with given particles,
it takes the history dict of the closest simulation run
and shows it in reverse all the way to height 100
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def gen_reverse_animation(history, save_path=None):
    """
    Create animation showing particles going backward in time (up in altitude).
    Shows how particles "grow" and speed up as we reverse the ablation process.

    Args:
        original_particles: list of particle dicts at original_height (not used, kept for signature)
        original_height: height where particles were observed (km)
        history: aggregated history dict from shooter_mcgavin
                 Format: dict[height] -> {'radii': [...], 'velocities': [...], 'n_particles': ...}
        save_path: optional path to save animation (e.g., 'reverse.gif')

    Returns:
        animation object
    """
    if not history:
        print("No history to animate!")
        return None

    n_particles = max(h.get('n_particles', 0) for h in history.values())
    print(f"Animating {n_particles} particles going backward in time...")

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter([], [], s=2, alpha=0.6, c='blue')

    ax.set_xlabel("Radius (microns)", fontsize=12)
    ax.set_ylabel("Velocity (km/s)", fontsize=12)
    ax.set_title("Reverse Trajectory: Particles Growing & Accelerating", fontsize=14)

    # Set axis limits
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 72.5)

    # Heights sorted from LOWEST to HIGHEST (reverse direction)
    # Animation will show particles going from low altitude UP to 100 km
    heights = sorted(history.keys())  # Ascending order

    # Add text annotation for particle count
    particle_count_text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                                   verticalalignment='top', fontsize=10,
                                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    def update(frame):
        """Update function for animation."""
        height = heights[frame]

        radii = history[height]["radii"]
        velocities = history[height]["velocities"]

        # Create scatter plot data
        if radii:
            points = np.column_stack([radii, velocities])
        else:
            points = np.empty((0, 2))

        scatter.set_offsets(points)

        # Update title and particle count
        ax.set_title(f"Reverse Trajectory | Height: {height:.1f} km | "
                     f"Going {'DOWN' if frame == 0 else 'UP'} in altitude",
                     fontsize=14)
        particle_count_text.set_text(f'Particles: {len(radii)}')

        return scatter, particle_count_text

    # Create animation
    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(heights),
        interval=100,  # 100ms per frame
        blit=True,
        repeat=True
    )

    # Save or display
    if save_path:
        print(f"Saving reverse animation to {save_path}...")
        if save_path.endswith('.gif'):
            anim.save(save_path, writer='pillow', fps=10)
        elif save_path.endswith('.mp4'):
            anim.save(save_path, writer='ffmpeg', fps=10)
        else:
            print(f"Warning: Unknown file format. Supported: .gif, .mp4")
            anim.save(save_path, writer='pillow', fps=10)
        print("Reverse animation saved!")
    else:
        plt.show()

    return anim


def gen_side_by_side_animation(forward_history, reverse_history, save_path=None):
    """
    Create side-by-side animation showing forward and reverse processes.

    Args:
        forward_history: history from main.py (descending)
        reverse_history: aggregated history from shooter_mcgavin
        save_path: optional path to save

    Returns:
        animation object
    """

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    scatter1 = ax1.scatter([], [], s=2, alpha=0.6, c='red')
    scatter2 = ax2.scatter([], [], s=2, alpha=0.6, c='blue')

    # Setup left plot (forward)
    ax1.set_xlabel("Radius (microns)")
    ax1.set_ylabel("Velocity (km/s)")
    ax1.set_title("Forward: Particles Falling & Ablating")
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 72.5)

    # Setup right plot (reverse)
    ax2.set_xlabel("Radius (microns)")
    ax2.set_ylabel("Velocity (km/s)")
    ax2.set_title("Reverse: Particles Rising & Growing")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 72.5)

    # Get heights
    forward_heights = sorted(forward_history.keys(), reverse=True)  # Descending
    reverse_heights = sorted(reverse_history.keys())  # Ascending

    max_frames = max(len(forward_heights), len(reverse_heights))

    def update(frame):
        """Update both plots."""
        # Update forward plot (if frame available)
        if frame < len(forward_heights):
            h_forward = forward_heights[frame]
            radii_f = forward_history[h_forward]["radii"]
            vels_f = forward_history[h_forward]["velocities"]

            if radii_f:
                points_f = np.column_stack([radii_f, vels_f])
            else:
                points_f = np.empty((0, 2))

            scatter1.set_offsets(points_f)
            ax1.set_title(f"Forward | h={h_forward:.1f} km | n={len(radii_f)}")

        # Update reverse plot (if frame available)
        if frame < len(reverse_heights):
            h_reverse = reverse_heights[frame]
            radii_r = reverse_history[h_reverse]["radii"]
            vels_r = reverse_history[h_reverse]["velocities"]

            if radii_r:
                points_r = np.column_stack([radii_r, vels_r])
            else:
                points_r = np.empty((0, 2))

            scatter2.set_offsets(points_r)
            ax2.set_title(f"Reverse | h={h_reverse:.1f} km | n={len(radii_r)}")

        return scatter1, scatter2

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=max_frames,
        interval=100,
        blit=True,
        repeat=True
    )

    if save_path:
        print(f"Saving side-by-side animation to {save_path}...")
        if save_path.endswith('.gif'):
            anim.save(save_path, writer='pillow', fps=10)
        elif save_path.endswith('.mp4'):
            anim.save(save_path, writer='ffmpeg', fps=10)
        print("Animation saved!")
    else:
        plt.show()

    return anim
