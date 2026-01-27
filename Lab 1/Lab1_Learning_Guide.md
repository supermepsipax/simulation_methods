# Lab 1 Learning Guide: Numerical Methods for the Harmonic Oscillator

## Table of Contents
1. [The Big Picture](#the-big-picture)
2. [Understanding the Harmonic Oscillator](#understanding-the-harmonic-oscillator)
3. [Why Different Methods Matter](#why-different-methods-matter)
4. [Key Concept: Energy Preservation](#key-concept-energy-preservation)
5. [Understanding the Demo Code](#understanding-the-demo-code)
6. [Method 1: Symplectic Euler](#method-1-symplectic-euler)
7. [Method 2: Implicit Midpoint Method](#method-2-implicit-midpoint-method)
8. [Method 3: Verlet's Method (Optional)](#method-3-verlets-method-optional)
9. [MATLAB Implementation Patterns](#matlab-implementation-patterns)
10. [Interpreting Your Results](#interpreting-your-results)

---

## The Big Picture

### Why Are We Doing This?

When we have a differential equation describing physical motion (like a swinging pendulum or vibrating spring), we often can't solve it with pen and paper. Instead, we **approximate** the solution by stepping through time in small increments.

**The core question of this lab:** Not all stepping methods are created equal. Some methods work beautifully for certain problems and terribly for others. The harmonic oscillator is a *perfect test case* because:
- We know the exact answer (so we can check our approximations)
- It has a conserved quantity (energy) that good methods should preserve
- It reveals fundamental differences between numerical methods

### What You'll Discover

By the end of this lab, you should understand *why* standard Euler methods fail for oscillatory problems, and *what special properties* make symplectic methods succeed.

---

## Understanding the Harmonic Oscillator

### The Physical System

Imagine a mass on a spring, or a simple pendulum with small swings. The position oscillates back and forth forever (in the ideal, frictionless case).

**The ODE:**
```
d²x/dt² + ω²x = 0
```

Here `ω` (omega) controls how fast the oscillation happens. Larger ω = faster oscillation.

### Why Convert to a System?

Most numerical methods work with **first-order** ODEs, but we have a **second-order** ODE. The standard trick:

**Introduce a helper variable for velocity:**
- Let `u₁ = x` (position)
- Let `u₂ = dx/dt` (velocity)

Now we have two first-order equations instead of one second-order equation:
```
du₁/dt = u₂         (position changes according to velocity)
du₂/dt = -ω²u₁      (velocity changes according to the restoring force)
```

### The Matrix Form

MATLAB loves matrices. The system can be written as:
```
du/dt = A * u
```
where `A = [0, 1; -ω², 0]` and `u = [u₁; u₂]`

This is useful because matrix operations in MATLAB are clean and efficient.

---

## Why Different Methods Matter

### The Explicit Euler Problem

Look at the demo code's explicit Euler:
```matlab
y = y + h*A*y;
```

This says: "Use the current slope to predict the next position."

**What goes wrong?** Each step, the method overshoots slightly. For oscillators, these overshoots *compound* — the orbit spirals outward. The numerical "energy" grows without bound.

### The Implicit Euler Problem

The demo's implicit Euler:
```matlab
y = (eye(2) - h*A) \ y;
```

This says: "Find the point where, if we were there, the slope would bring us back to where we are."

**What goes wrong?** Each step, the method undershoots. The orbit spirals inward. The numerical "energy" decays to zero.

### The Symplectic Idea

*Symplectic* methods are designed to preserve the **geometric structure** of Hamiltonian systems (systems with conserved energy). They don't perfectly preserve energy, but they keep it bounded — it oscillates around the true value rather than drifting away.

**Key insight:** For long-time simulations of physical systems, a method that's "approximately right forever" beats a method that's "very accurate but drifts."

---

## Key Concept: Energy Preservation

### What is the Hamiltonian?

For the harmonic oscillator, total energy is:
```
H(t) = (velocity)² + ω²(position)²
     = u₂² + ω²u₁²
```

For the exact solution, `H(t) = ω²` (constant for all time).

### Why Energy Matters

- **Drifting energy** = your numerical solution is doing something physically impossible
- **Preserved energy** = your solution stays on the correct "energy surface" in phase space
- **Bounded but oscillating energy** = acceptable for symplectic methods

### Computing Energy in MATLAB

At each timestep, after updating `u₁` and `u₂`:
```matlab
E_current = omega^2 * u1^2 + u2^2;
```

Store these values to plot energy vs. time.

---

## Understanding the Demo Code

Before modifying the demo, understand what it does:

### Structure of the Demo

1. **Exact solution** (lines 3-24): Computes `cos(ωt)` and `-ω*sin(ωt)` analytically
2. **Explicit Euler** (lines 25-49): Forward stepping with `y = y + h*A*y`
3. **Implicit Euler** (lines 50-71): Backward stepping with matrix solve

### Key Variables to Understand

| Variable | Purpose |
|----------|---------|
| `h` | Time step size (0.1 in demo) |
| `N` | Number of steps (100 in demo) |
| `omega` | Oscillation frequency (2 in demo) |
| `Y` | Stores all solution values (grows each iteration) |
| `T` | Stores all time values |
| `E` | Stores energy at each step |

### The Growing Arrays Pattern

```matlab
Y = [y'];          % Initialize with first row
for k = 1:N
    % ... update y ...
    Y = [Y; y'];   % Append new row
end
```

This builds up your solution history for plotting.

---

## Method 1: Symplectic Euler

### The Key Difference

Look carefully at the formulas:
```
u₁,ₖ = u₁,ₖ₋₁ + h·u₂,ₖ₋₁      (uses OLD velocity)
u₂,ₖ = u₂,ₖ₋₁ - h·ω²·u₁,ₖ     (uses NEW position!)
```

Notice `u₁,ₖ` appears on the right side of the second equation. This is **semi-implicit**: position is updated explicitly, then velocity uses the *just-computed* new position.

### Why This Matters

This subtle change — using the updated position immediately — creates a method that respects the symplectic structure. It's a small modification to explicit Euler but with dramatically different long-term behavior.

### Implementation Approach

You need to update **component by component**, not using matrix multiplication:
```matlab
% Pseudocode structure:
u1_new = u1_old + h * u2_old;
u2_new = u2_old - h * omega^2 * u1_new;  % Note: u1_NEW here!
```

### Questions to Guide Your Implementation

- How do you initialize `u1` and `u2`?
- Where do you store results for plotting?
- When do you compute energy — before or after updating?

---

## Method 2: Implicit Midpoint Method

### The Concept

Instead of evaluating the derivative at the beginning (explicit) or end (implicit), evaluate it at the **midpoint**:

```
uₖ = uₖ₋₁ + h·A·((uₖ + uₖ₋₁)/2)
```

The unknown `uₖ` appears on both sides — this is an **implicit** method requiring you to solve for `uₖ`.

### Rearranging to Solve

The equation has `uₖ` on both sides. You need to isolate `uₖ`:

Starting from:
```
uₖ = uₖ₋₁ + (h/2)·A·uₖ + (h/2)·A·uₖ₋₁
```

Collect terms with `uₖ` on the left side...

**Hint:** You'll end up with something of the form:
```
(something) * uₖ = (something else) * uₖ₋₁
```

Which you solve as:
```matlab
u_new = left_matrix \ (right_matrix * u_old);
```

### The Identity Matrix

In MATLAB, `eye(2)` creates a 2×2 identity matrix:
```
[1  0]
[0  1]
```

You'll need this when rearranging the implicit equation.

### Questions to Guide Your Implementation

- What matrices appear on the left and right of your equation?
- Can you precompute these matrices outside the loop?
- How does this compare to the implicit Euler solve in the demo?

---

## Method 3: Verlet's Method (Optional)

### A Different Philosophy

Verlet works directly with the **second-order** equation instead of converting to a first-order system. It uses positions at three time levels: `k-1`, `k`, and `k+1`.

### Central Differences

The idea: approximate the second derivative using:
```
d²u₁/dt² ≈ (u₁,ₖ₊₁ - 2u₁,ₖ + u₁,ₖ₋₁) / h²
```

Substituting into the ODE and rearranging gives the update formula.

### The Bootstrap Problem

Verlet needs `u₁,₀` AND `u₁,₁` to start (two previous values). But we only have one initial condition. The assignment gives you a special formula to compute `u₁,₁`.

### Velocity Recovery

Verlet naturally gives you positions. Velocity is recovered using:
```
u₂,ₖ = (u₁,ₖ₊₁ - u₁,ₖ₋₁) / (2h)
```

### Questions to Guide Your Implementation

- How is the loop index different here? (Notice k goes to N-1, not N)
- What do you store: `u₁,ₖ₋₁`, `u₁,ₖ`, or both?
- When can you compute velocity values?

---

## MATLAB Implementation Patterns

### Initializing Your Solution Storage

```matlab
% Scalar storage (for component-wise methods like Symplectic Euler)
u1 = 1;  % Initial position
u2 = 0;  % Initial velocity
U1 = [u1];  % Will grow to store all positions
U2 = [u2];  % Will grow to store all velocities

% Vector storage (for matrix methods like IMM)
u = [1; 0];  % Column vector
U = [u'];    % First row of solution matrix
```

### The Time-Stepping Loop

```matlab
for k = 1:N
    % 1. Update your solution (method-specific)

    % 2. Store the new values
    U1 = [U1; u1];  % or U = [U; u'];

    % 3. Update time
    t = t + h;
    T = [T; t];

    % 4. Compute and store energy
    E = [E; omega^2*u1^2 + u2^2];
end
```

### Plotting Patterns

**Trajectories (position and velocity vs. time):**
```matlab
plot(T, U1, T, U2)
xlabel('t')
ylabel('u1, u2')
title('Your Title Here')
```

**Phase Portrait (velocity vs. position):**
```matlab
plot(U1, U2)
axis('equal')  % Important: makes circles look like circles!
xlabel('u1 (position)')
ylabel('u2 (velocity)')
```

**Energy vs. Time:**
```matlab
plot(T, E)
xlabel('t')
ylabel('Energy')
% Add a reference line for exact energy:
hold on
plot(T, omega^2 * ones(size(T)), 'r--')
hold off
```

### Solving Linear Systems

For implicit methods, you solve `M * x = b`:
```matlab
x = M \ b;  % Backslash operator
```

**Never** use `inv(M) * b` — backslash is more accurate and faster.

---

## Interpreting Your Results

### What to Look For: Trajectories

| Observation | Meaning |
|-------------|---------|
| Amplitude grows over time | Energy is being added (unstable) |
| Amplitude decays over time | Energy is being lost (dissipative) |
| Amplitude stays constant | Energy is preserved (good!) |
| Phase shifts over time | Frequency error (may be acceptable) |

### What to Look For: Phase Portrait

| Shape | Interpretation |
|-------|----------------|
| Outward spiral | Explicit Euler behavior — method is unstable |
| Inward spiral | Implicit Euler behavior — method is dissipative |
| Closed ellipse | Energy is preserved — symplectic behavior |
| Ellipse that doesn't quite close | Bounded energy oscillation (acceptable) |

The exact solution traces a perfect ellipse. Compare your numerical solutions to this.

### What to Look For: Energy Plot

| Pattern | Interpretation |
|---------|----------------|
| Exponential growth | Serious instability |
| Exponential decay | Artificial damping |
| Constant horizontal line | Perfect energy preservation |
| Small oscillations around constant | Symplectic behavior (bounded error) |

### Questions to Answer in Your Analysis

1. **Symplectic Euler:** Does the energy stay bounded? What shape is the phase portrait?

2. **Implicit Midpoint:** How does energy preservation compare to Symplectic Euler?

3. **Comparison:** Which method would you trust for a simulation running for 1000 time units instead of 10?

---

## Common Pitfalls and Debugging

### Symptom: Solution Explodes Immediately
- Check your initial conditions
- Check the sign of ω² term (should be negative in the u₂ equation)
- Check you're using `h` not `h²` where appropriate

### Symptom: Results Look Like Explicit/Implicit Euler
- For Symplectic Euler: Make sure you use `u1_new` (not `u1_old`) in the u₂ update
- For IMM: Make sure you rearranged the equation correctly

### Symptom: Energy Plot Shows Constant but Wrong Value
- Check your energy formula: `E = omega^2 * u1^2 + u2^2`
- Make sure you're storing energy *after* the update, not before

### Symptom: Phase Portrait is Stretched/Squished
- Use `axis('equal')` to get proper aspect ratio
- Check that you're plotting position on x-axis, velocity on y-axis

---

## Quick Reference: Method Formulas

### Symplectic Euler (Component Form)
```
u₁,ₖ = u₁,ₖ₋₁ + h·u₂,ₖ₋₁
u₂,ₖ = u₂,ₖ₋₁ - h·ω²·u₁,ₖ    ← uses NEW position
```

### Implicit Midpoint (Vector Form)
```
uₖ = uₖ₋₁ + h·A·(uₖ + uₖ₋₁)/2

Rearranged: (I - h/2·A)·uₖ = (I + h/2·A)·uₖ₋₁
```

### Verlet (Second-Order Form)
```
u₁,ₖ₊₁ = 2u₁,ₖ - u₁,ₖ₋₁ - h²ω²u₁,ₖ
u₂,ₖ = (u₁,ₖ₊₁ - u₁,ₖ₋₁)/(2h)

Bootstrap: u₁,₁ = u₁,₀ + h·u₂,₀ - (h²/2)·ω²·u₁,₀
```

---

## Further Thinking (Not Required)

If you want deeper understanding:

- **Why "symplectic"?** The word comes from Greek for "intertwined." Symplectic methods preserve the intertwined relationship between position and momentum in phase space.

- **Why does the midpoint work?** The midpoint method is "symmetric" in time — running it forward then backward returns exactly to the start. This time-symmetry relates to energy conservation.

- **Why does Verlet work?** Verlet is secretly a symplectic method too! Its clever use of central differences preserves the geometric structure.

---

Good luck with the lab! Remember: the goal is understanding *why* these methods behave differently, not just getting code that runs.
