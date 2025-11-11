"""
Given a height with certain particles defined by radius and velocity,
use the shooting method to find the optimal initial conditions.
Would return the history of the optimal sim run
and call the reverse_animation function.
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.constants import g
from processing.init_particles import generate_particles as gen_p
import processing.friends as friend
from plots.reverse_animation import gen_reverse_animation as gen_r


def simulate_particle_forward(r_initial, v_initial, h_start=100, h_end=50, step_size=0.1):
    """
    Simulate a single particle falling from h_start to h_end.
    Returns history of the simulation.

    Parameters:
        r_initial: initial radius in microns
        v_initial: initial velocity in km/s
        h_start: starting height in km
        h_end: ending height in km
        step_size: height step size in km (positive value)

    Returns:
        dict with 'final_state': (r, v) or None if burned up
                    'history': dict[height] -> {'r': radius, 'v': velocity}
    """
    particle = {
        "r": r_initial,
        "v": v_initial,
        "active": True
    }

    history = {}
    heights = np.arange(h_start, h_end, -step_size)

    for height in heights:
        # Record current state
        history[height] = {
            'r': particle['r'],
            'v': particle['v']
        }

        if not particle["active"] or particle["r"] <= 0:
            return {'final_state': None, 'history': history}

        m = friend.calculate_mass(particle)
        height_meters = height * 1000
        drag = friend.calculate_drag(particle, height_meters)
        v_ms = particle["v"] * 1000

        # Stop if velocity too low
        if particle["v"] < 0.01:
            return {'final_state': None, 'history': history}

        # Update velocity
        if v_ms > 10:
            dv_dh = ((drag / m) - g) / v_ms
            dv_kms = dv_dh * (-step_size)  # negative step
            particle["v"] = particle["v"] + dv_kms

            if particle["v"] <= 0:
                return {'final_state': None, 'history': history}

        # Update radius (ablation)
        dr_dh = friend.calculate_ablation(particle, drag)
        new_radius = particle["r"] - (dr_dh * (-step_size))

        if new_radius < 0.001:
            return {'final_state': None, 'history': history}
        else:
            particle["r"] = new_radius

    # Final state
    final_state = (particle["r"], particle["v"])
    history[h_end] = {'r': particle['r'], 'v': particle['v']}

    return {'final_state': final_state, 'history': history}


def objective_function(initial_state, target_height, target_r, target_v, h_start=100):
    """
    Objective function for optimization.
    Returns error between simulated and target state.
    """
    r_0, v_0 = initial_state

    # Sanity checks
    if r_0 <= 0 or v_0 <= 0 or r_0 > 1000 or v_0 > 72.5:
        return [1e10, 1e10]  # Infeasible

    # Run forward simulation
    result = simulate_particle_forward(r_0, v_0, h_start=h_start, h_end=target_height)

    if result['final_state'] is None:
        # Particle burned up - penalize heavily
        return [1e10, 1e10]

    r_final, v_final = result['final_state']

    # Compute errors (normalized)
    error_r = (r_final - target_r) / max(target_r, 0.1)  # Relative error
    error_v = (v_final - target_v) / max(target_v, 0.1)  # Relative error

    return [error_r, error_v]


def find_initial_conditions_single(target_height, target_r, target_v, h_start=100, initial_guess=None):
    """
    Find initial radius and velocity at h_start for a single particle.

    Returns:
        dict with 'r_0', 'v_0', 'success', 'residual', 'history'
    """
    # Automatic initial guess: particle was larger and faster
    if initial_guess is None:
        r_guess = target_r * 1.3
        v_guess = target_v * 1.4
        initial_guess = [r_guess, v_guess]

    # Bounds: reasonable physical limits
    bounds = ([0.1, 0.1], [1000, 72.5])  # r: 0.1-1000 μm, v: 0.1-72.5 km/s

    # Solve using least squares
    result = least_squares(
        objective_function,
        x0=initial_guess,
        args=(target_height, target_r, target_v, h_start),
        bounds=bounds,
        method='trf',
        ftol=1e-6,
        xtol=1e-6,
        max_nfev=500
    )

    r_0, v_0 = result.x

    # Get full history by running forward simulation with optimal params
    sim_result = simulate_particle_forward(r_0, v_0, h_start=h_start, h_end=target_height)

    return {
        'r_0': r_0,
        'v_0': v_0,
        'success': result.success,
        'residual': np.linalg.norm(result.fun),
        'history': sim_result['history']
    }


def shooter_mcgavin(particles, height, h_start=100, verbose=True):
    """
    Find initial conditions at h_start for particles observed at given height.

    Parameters:
        particles: list of dicts with 'r' and 'v' at target height
        height: height where particles are observed (km)
        h_start: initial height to find conditions for (default 100 km)
        verbose: print progress

    Returns:
        dict with:
            'initial_particles': list of particles at h_start
            'history': aggregated history dict[height] -> {'radii': [...], 'velocities': [...], 'n_particles': ...}
            'n_success': number of successfully inverted particles
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"SHOOTING METHOD: Finding initial conditions at h={h_start} km")
        print(f"Given {len(particles)} particles at h={height} km")
        print(f"{'='*70}\n")

    initial_particles = []
    individual_histories = []
    n_success = 0

    for i, particle in enumerate(particles):
        target_r = particle['r']
        target_v = particle['v']

        result = find_initial_conditions_single(
            height, target_r, target_v, h_start=h_start
        )

        if result['success']:
            initial_particles.append({
                'r': result['r_0'],
                'v': result['v_0'],
                'active': True
            })
            individual_histories.append(result['history'])
            n_success += 1
        else:
            # Store failed result with original values
            initial_particles.append({
                'r': target_r,
                'v': target_v,
                'active': False
            })
            individual_histories.append({})

        if verbose and (i + 1) % 10 == 0:
            success_rate = 100 * n_success / (i + 1)
            print(f"Progress: {i + 1}/{len(particles)} particles | "
                  f"Success rate: {success_rate:.1f}%")

    # Aggregate individual histories into single history dict (like main.py format)
    aggregated_history = {}
    for particle_history in individual_histories:
        if not particle_history:  # Skip failed inversions
            continue

        for h, state in particle_history.items():
            if h not in aggregated_history:
                aggregated_history[h] = {'radii': [], 'velocities': []}

            aggregated_history[h]['radii'].append(state['r'])
            aggregated_history[h]['velocities'].append(state['v'])

    # Add particle counts to each height
    for h in aggregated_history:
        aggregated_history[h]['n_particles'] = len(aggregated_history[h]['radii'])

    if verbose:
        print(f"\n{'='*70}")
        print(f"RESULTS:")
        print(f"  Successfully inverted: {n_success}/{len(particles)} particles ({100*n_success/len(particles):.1f}%)")

        if n_success > 0:
            successful_initial = [p for p in initial_particles if p['active']]
            r_values = [p['r'] for p in successful_initial]
            v_values = [p['v'] for p in successful_initial]

            print(f"\nInitial distribution at h={h_start} km:")
            print(f"  Radius:   {np.min(r_values):.2f} - {np.max(r_values):.2f} μm "
                  f"(mean: {np.mean(r_values):.2f})")
            print(f"  Velocity: {np.min(v_values):.2f} - {np.max(v_values):.2f} km/s "
                  f"(mean: {np.mean(v_values):.2f})")
        print(f"{'='*70}\n")

    return {
        'initial_particles': initial_particles,
        'history': aggregated_history,
        'n_success': n_success
    }


# Example usage
if __name__ == "__main__":
    # Example particles at h=50 km
    partycles = gen_p(100)
    # Example height
    height = 50

    # Find initial conditions at h=100 km
    result = shooter_mcgavin(partycles, height, h_start=100)

    # The result contains:
    # - initial_particles: what the particles were at h=100 km
    # - history: single aggregated trajectory dict from h=100 to h=50
    #   Format: history[height] = {'radii': [...], 'velocities': [...], 'n_particles': N}
    # - n_success: how many inversions succeeded

    gen_r(result["history"], save_path="shooter.gif")

