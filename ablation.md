# Ablation Model for Cosmic Dust Particles

## Overview

Ablation is the process by which particles lose mass (and thus radius) as they pass through Earth's atmosphere due to frictional heating. This is a critical phenomenon for understanding which particles reach the surface versus which completely vaporize.

## Physical Mechanism

### Energy Transfer
As a particle moves through the atmosphere:
1. **Drag force** opposes the particle's motion
2. **Kinetic energy** is converted to heat through friction with air molecules
3. Only a **fraction (ε)** of this energy heats the particle (rest radiates away or heats air)
4. When particle temperature exceeds melting/vaporization point, material evaporates

### Key Physics Equation

The heating rate from drag:
```
Heating power = ε * F_drag * v
```

Where:
- `ε` ≈ 0.1 to 0.5 is the heat transfer efficiency
- `F_drag` is the drag force (N)
- `v` is particle velocity (m/s)

## Mathematical Model

### Mass Loss Rate

The rate at which mass evaporates depends on energy input:
```
dm/dt = -(ε * F_drag * v) / L_vaporization
```

Where `L_vaporization` is the heat of vaporization (J/kg) - the energy needed to evaporate 1 kg of material.

### Converting to Radius Change

Since mass and radius are related by:
```
m = (4/3)πr³ρ
```

Taking the derivative:
```
dm/dt = 4πr²ρ * dr/dt
```

Therefore:
```
dr/dt = (dm/dt) / (4πr²ρ)
     = -(ε * F_drag * v) / (4πr²ρ * L_vaporization)
```

### Height as Independent Variable

Since we use height `h` as our independent variable (not time `t`):
```
dr/dh = (dr/dt) / (dh/dt)
      = (dr/dt) / (-v)
      = (ε * F_drag) / (4πr²ρ * L_vaporization)
```

Note: The velocity cancels out when converting from time to height!

## Simplifying Assumptions

For the SCUDEM model, we make the following reasonable assumptions:

### 1. Uniform Particle Density
- **Assumption**: All particles have density ρ = 2 g/cm³ = 2000 kg/m³
- **Justification**: Cosmic dust is typically a mixture of silicates and metals with average density ~2 g/cm³
- **Impact**: Simplifies mass calculations without significantly affecting results

### 2. Ignore Chemical Composition
- **Assumption**: Treat all particles as generic "stony" material
- **Justification**: Most cosmic dust is similar composition (olivine, pyroxene, iron)
- **Impact**: Allows us to use a single value for heat of vaporization

### 3. Constant Heat of Vaporization
- **Assumption**: L_vap = 6000 kJ/kg for all particles
- **Justification**: This is typical for silicate materials
- **Impact**: Ablation rate scales uniformly with particle properties

### 4. Energy Efficiency Factor
- **Assumption**: ε = 0.3 (30% of drag energy heats particle)
- **Justification**: Based on meteorite entry studies; 70% radiates or heats air
- **Impact**: Can be tuned to match observational data

## Implementation Details

### Function Signature
```python
def calculate_ablation(particle, F_drag):
    """
    Returns: dr/dh in microns/km (negative = shrinking)
    """
```

### Physical Constants Used
```python
epsilon = 0.3           # Heat transfer efficiency (dimensionless)
L_vap = 6e6            # Heat of vaporization (J/kg) for silicates
rho_particle = 2000    # Particle density (kg/m³)
```

### Unit Conversions
The function handles:
- Input: radius in **microns**, velocity in **km/s**
- Internal calculations: **SI units** (meters, m/s, kg, J)
- Output: radius change in **microns/km**

### Edge Cases
The function handles:
- Very slow particles (v < 0.001 km/s): Returns 0 ablation
- Very small particles (r < 1e-10 m): Returns 0 ablation
- Prevents division by zero

## Physical Interpretation

### When does ablation matter?

**High ablation** occurs when:
- ✓ High velocity (v²) → more drag → more heating
- ✓ High air density → more drag
- ✓ Small radius → higher surface-to-volume ratio → easier to heat

**Low ablation** occurs when:
- ✓ Low velocity → less heating
- ✓ High altitude (low air density) → less drag
- ✓ Large radius → harder to heat entire particle

### Typical Behavior

For cosmic dust entering at ~20 km/s:
- **Large particles (>100 μm)**: Survive to lower altitudes, some reach surface
- **Medium particles (10-100 μm)**: Partially ablate, fragments may survive
- **Small particles (<10 μm)**: Completely ablate at high altitude (>60 km)

## Connection to Euler Method

In the main simulation loop, ablation is updated using Euler's method:

```python
r_(n+1) = r_n + (dr/dh) * Δh
```

Where:
- `dr/dh` comes from `calculate_ablation()`
- `Δh = -0.5 km` (stepping downward)
- `r` decreases over time (ablation shrinks particle)

## Validation and Tuning

The ablation model can be validated and tuned by:

1. **Comparing total mass deposition** to observed data (~5000 tons/year reaches Earth)
2. **Size distribution at ground** should match meteorite samples
3. **Altitude of ablation** should match meteor observations (typically 70-100 km)
4. **Adjusting ε** if needed to match empirical data

## References and Further Reading

- Heat of vaporization for silicates: 5-7 MJ/kg
- Meteoroid entry physics: Ceplecha et al. (1998)
- Cosmic dust flux: Plane (2012), Love & Brownlee (1993)
- Atmospheric entry modeling: Baldwin & Sheaffer (1971)

## Alternative Models

### Empirical Formula
A simpler empirical approach used in some models:
```
dr/dh = -σ * ρ_air * v² / ρ_particle
```
Where σ ≈ 10⁻¹² to 10⁻⁹ is an empirical ablation coefficient.

### Temperature-Threshold Model
More complex models only apply ablation when particle temperature exceeds melting point:
```
if T_particle > T_melting:
    apply ablation
else:
    no ablation
```

Our energy-based model implicitly assumes particles are hot enough to ablate when drag is significant.
