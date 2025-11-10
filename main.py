"""
FUN FUN FUN TIME 
FUN FUN RUN TIME
FUN RUN RUN TIME
RUN RUN RUN TIME
"""

import numpy as np
from scipy.constants import g

from processing.init_particles import generate_particles as gen_p
import processing.friends as friend
from plots.distribution_graph import plot_distribution as plot_d

def go_do_thing_time():
    acc = 0
    history = {}
    particles = gen_p(100000)
    plot_d(particles)
    step_size = -0.5
    heights = np.arange(100, step_size, step_size)
    for height in heights:
        for particle in particles:
            m = friend.calculate_mass(particle)
            drag = friend.calculate_drag(particle, height)

            dv_dh = ((drag / m) - g) / particle["v"]
            particle["v"] = particle["v"] + dv_dh * step_size
            # if particle["v"] < 0.0001:
            #     particle["v"] = 0.0001

            dr_dh = friend.calculate_ablation(particle, drag)
            particle["r"] = particle["r"] + dr_dh * step_size
            if particle["r"] < 0.001:
                particle["active"] = False

        history[height] = {
        'radii': [p['r'] for p in particles if p['active']],
        'velocities': [p['v'] for p in particles if p['active']],
        'n_particles': len([p for p in particles if p['active']])
        }
        acc += 1
        print(f"made it {acc}")
    final_particles = [p for p in particles if p["active"]]
    plot_d(final_particles, title = "Final Destination")
    


go_do_thing_time()