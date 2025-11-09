# math_wizards_scudem_2025

We da best

## Setup

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

To activate the environment in the future:
```bash
source venv/bin/activate
```

To deactivate when done:
```bash
deactivate
```

## Problem Overview


Karman Line is where space starts at 100 km above sea level, so that is our initial height
We are also saying that at 100 km air density is practically 0


General Structure:

Given particles (a list of dictionaries) and height step size.
Create a list of heights. 
for height in heights:
    calculate_air_density(height)
    for particle in particles:
        calculate_ke(particle, air_density)
        calculate_drag(particle, air_density)
        sum_of_forces(ke, drag)
        calculate_ablation(particle, drag, sum_of_forces, air_density)
        check_fragmentation(particle, air_density)
    history[height] = {
        'radii': [p['r'] for p in particles if p['active']],
        'velocities': [p['v'] for p in particles if p['active']],
        'n_particles': len([p for p in particles if p['active']])
    }
Create animation of distribution of particles at each height.
Pretty pictures





What has to be done?
* Understand the Problem
* Make the model
* Write the simulation of the system
* Test simulation with existing data ([this](https://www.scientificamerican.com/article/antarctic-study-shows-how-much-space-dust-hits-earth-every-year/) for example)