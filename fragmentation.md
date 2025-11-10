# Fragmentation Model for Cosmic Dust Particles

## Overview

Fragmentation is the process where particles break apart due to mechanical stress. For particles entering Earth's atmosphere at high velocity, **ram pressure** (dynamic pressure from moving through air) can exceed the particle's structural strength and cause it to shatter.

## Physical Mechanism

### Ram Pressure

As a particle moves through the atmosphere, it experiences dynamic pressure:

```
P_ram = (1/2) * ρ_air * v²
```

This is the "impact pressure" felt by the leading edge of the particle.

### Material Strength

Cosmic dust particles have **tensile strength** (resistance to being pulled apart):
- Loose aggregates: ~10⁴ to 10⁵ Pa
- Compact silicates: ~10⁵ to 10⁶ Pa
- Iron-nickel: ~10⁷ to 10⁸ Pa

### Fragmentation Criterion

```
IF P_ram > σ_strength THEN particle fragments
```

## When Does Fragmentation Occur?

### Typical Scenarios

For cosmic dust (strength ~ 10⁵ Pa):

**High altitude (h > 60 km):**
- Low air density (ρ ~ 10⁻⁴ kg/m³)
- Even at v = 20 km/s: P_ram ~ 2×10⁴ Pa
- **No fragmentation** (pressure too low)

**Medium altitude (h ~ 40 km):**
- Moderate air density (ρ ~ 10⁻² kg/m³)
- At v = 20 km/s: P_ram ~ 2×10⁶ Pa
- **Fragmentation likely!**

**Low altitude (h < 20 km):**
- High air density, but velocity reduced by drag
- If v slowed to < 5 km/s: P_ram ~ 10⁵ Pa
- **Borderline** (depends on remaining velocity)

### Critical Velocity

For a given strength and air density, there's a critical velocity:

```
v_critical = sqrt(2 * σ_strength / ρ_air)
```

Above this velocity, fragmentation occurs.

## Mathematical Model

### Conservation Laws

When a particle fragments into n pieces:

**Mass conservation:**
```
m_parent = Σ m_fragments
```

**Equal-mass assumption:**
For simplicity, assume equal-sized fragments:
```
m_fragment = m_parent / n
```

**Radius of fragments:**
Since m = (4/3)πr³ρ:
```
r_fragment = r_parent / n^(1/3)
```

### Fragment Velocities

Fragments inherit roughly the same velocity, with small variations:
```
v_fragment ≈ v_parent * (0.9 to 1.1)
```

Small spread accounts for:
- Different ejection angles
- Slight momentum transfer
- Turbulent effects

## Implementation Details

### Function Behavior

```python
fragments = check_fragmentation(particle, air_density)

if fragments:
    # Particle fragmented!
    particle['active'] = False  # Mark original as inactive
    for frag in fragments:
        particles.append(frag)  # Add new fragments
```

### Parameters to Tune

**Material strength** (σ_strength):
- Lower value (10⁴ Pa): Weak, porous aggregates → more fragmentation
- Higher value (10⁶ Pa): Strong, compact particles → less fragmentation

**Number of fragments**:
- Currently: 2-4 pieces (random)
- Could be made deterministic based on particle size or energy

### Edge Cases

The function handles:
- Particles below fragmentation threshold → returns empty list
- Random number of fragments (2-4) for variability
- Each fragment gets `'active': True` flag

## Integration with Main Simulation

### Example Loop Structure

```python
for height in heights:
    air_density = calculate_air_density(height)

    # Need to iterate carefully since we're adding particles
    i = 0
    while i < len(particles):
        particle = particles[i]

        if not particle['active']:
            i += 1
            continue

        # Check fragmentation FIRST (before other updates)
        fragments = check_fragmentation(particle, air_density)
        if fragments:
            particle['active'] = False  # Original breaks apart
            particles.extend(fragments)  # Add fragments to list
            i += 1
            continue

        # Normal updates for intact particles
        F_drag = calculate_drag(particle, height)

        # Update velocity (Euler method)
        dv_dh = ((F_drag/m - g) / particle['v'])
        particle['v'] += dv_dh * step_size

        # Update radius (ablation)
        dr_dh = calculate_ablation(particle, F_drag)
        particle['r'] += dr_dh * step_size

        # Check if particle too small
        if particle['r'] < 0.01:  # microns
            particle['active'] = False

        i += 1
```

## Physical Interpretation

### Why does fragmentation matter?

**Without fragmentation:**
- Large particles survive to low altitudes
- Fewer total particles
- Mass concentrated in original sizes

**With fragmentation:**
- Large particles break into smaller pieces at high altitude
- Many more small fragments
- Fragments ablate faster (higher surface/volume ratio)
- More realistic size distribution at ground level

### Cascade Effect

Fragmentation can create a cascade:
1. Large particle enters atmosphere
2. Ram pressure fragments it at ~60 km altitude
3. Fragments have smaller radius → less strength per area
4. Some fragments fragment again at lower altitudes
5. Eventually create many tiny particles that completely ablate

This matches observations of meteor "flares" (bright bursts) as particles break apart.

## Validation

The fragmentation model can be validated by:

1. **Meteor observations**: Bright flares indicate fragmentation events
2. **Size distribution**: Compare predicted vs. observed size distributions
3. **Altitude of fragmentation**: Should match observed meteor brightening altitudes (typically 60-80 km)
4. **Mass deposition**: Total mass reaching surface should match observations (~5000 tons/year)

## Complexity vs. Benefit Trade-off

### Pros of including fragmentation:
- ✅ More physically realistic
- ✅ Matches observed meteor behavior (flares, multiple trails)
- ✅ Creates interesting particle population dynamics
- ✅ Shows sophisticated modeling to SCUDEM judges

### Cons of including fragmentation:
- ❌ Adds complexity to simulation loop (need to handle growing particle list)
- ❌ More parameters to tune (strength, fragment number)
- ❌ Harder to explain in 10-minute video
- ❌ Takes time to implement and debug

## Simplified Alternative

If fragmentation is too complex, you can mention it as "future work":

> "In our current model, we assume particles remain intact. A more sophisticated
> model would include fragmentation when ram pressure exceeds material strength
> (~10⁵ Pa for cosmic dust). This would create a cascade of smaller fragments
> and better match observed meteor flare patterns."

This shows judges you understand the physics without implementing it.

## References

- Ram pressure and fragmentation: Baldwin & Sheaffer (1971)
- Meteor fragmentation observations: Ceplecha et al. (1998)
- Material strength of meteoroids: Popova et al. (2011)
- Tensile strength of cosmic dust: Güttler et al. (2010)
