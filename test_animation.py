"""
Quick test for the animation function
Uses generate_particles and simulates changes over 100 height steps
"""

from plots.distribution_animation import kawaii_anime
from processing.init_particles import generate_particles
import numpy as np

# Generate initial particles at 100 km
particles = generate_particles(1000)

# Add 'active' flag to each particle
for p in particles:
    p['active'] = True

# Create history by simulating particles falling through atmosphere
history = {}

for height in np.arange(100, 0, -1):  # 100, 99, 98, ..., 1
    # Simulate ablation: randomly deactivate some particles (more at lower heights)
    ablation_rate = 0.005 * (101 - height)  # Increases as height decreases
    for p in particles:
        if p['active'] and np.random.random() < ablation_rate:
            p['active'] = False

    # Simulate radius reduction due to ablation
    for p in particles:
        if p['active']:
            p['r'] = max(0.01, p['r'] - np.random.uniform(0, 0.02))

    # Simulate velocity reduction due to drag
    for p in particles:
        if p['active']:
            p['v'] = max(1, p['v'] - np.random.uniform(0, 0.1))

    # Simulate fragmentation: occasionally split particles
    if height < 80 and np.random.random() < 0.1:
        active_particles = [p for p in particles if p['active']]
        if active_particles:
            parent = np.random.choice(active_particles)
            # Create 2 smaller fragments
            for _ in range(2):
                fragment = {
                    'r': parent['r'] * 0.5,
                    'v': parent['v'] * np.random.uniform(0.8, 1.2),
                    'active': True
                }
                particles.append(fragment)

    # Record state at this height
    history[height] = {
        'radii': [p['r'] for p in particles if p['active']],
        'velocities': [p['v'] for p in particles if p['active']],
        'n_particles': len([p for p in particles if p['active']])
    }

print(f"Created history with {len(history)} height steps")
print(f"Heights: {sorted(history.keys(), reverse=True)[:5]}... (showing first 5)")
print(f"Initial particles: {history[100]['n_particles']}")
print(f"Final particles: {history[1]['n_particles']}")

# Run the animation
anim = kawaii_anime(history, "test_animation.gif")
