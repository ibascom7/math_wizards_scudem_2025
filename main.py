"""
FUN FUN FUN TIME 
FUN FUN RUN TIME
FUN RUN RUN TIME
RUN RUN RUN TIME
"""

from processing.init_particles import generate_particles as gen_p

from plots.distribution_graph import plot_distribution as pd

def go_do_thing_time():
    particles = gen_p
    pd(particles)

go_do_thing_time()