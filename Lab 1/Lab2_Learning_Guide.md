# Lab 2 Learning Guide: Boundary Value Problem for Heat Conduction

## Table of Contents
1. [The Big Picture](#the-big-picture)
2. [Understanding the Physical Problem](#understanding-the-physical-problem)
3. [Key Concept: BVP vs IVP](#key-concept-bvp-vs-ivp)
4. [Finite Difference Discretization](#finite-difference-discretization)
5. [Building the Tridiagonal System](#building-the-tridiagonal-system)
6. [Handling Boundary Conditions](#handling-boundary-conditions)
7. [The Source Term Q(z)](#the-source-term-qz)
8. [MATLAB Implementation Patterns](#matlab-implementation-patterns)
9. [Interpreting Your Results](#interpreting-your-results)
10. [Common Pitfalls](#common-pitfalls)

---

## The Big Picture

### What's Different from Lab 1?

| Lab 1 (Harmonic Oscillator) | Lab 2 (Heat Conduction) |
|----------------------------|-------------------------|
| Initial Value Problem (IVP) | Boundary Value Problem (BVP) |
| Known: starting point | Known: both endpoints |
| March forward in time | Solve all points simultaneously |
| Time-stepping loop | Linear system solve |
| Variable: time t | Variable: position z |

### The Core Challenge

You know the temperature at both ends of the pipe (boundary conditions), but need to find the temperature everywhere in between. Unlike time-stepping where you compute one value after another, here **all interior temperatures depend on each other** — they must be solved simultaneously as a system of equations.

---

## Understanding the Physical Problem

### The Setup

Imagine a pipe of length L = 10:
- Fluid flows through at velocity v
- An electric coil heats the section between z = 1 and z = 3
- We want the steady-state temperature profile T(z)

```
z=0                    z=L
|----[===COIL===]------|
T₀=400    heating     T_out=300
          region
     a=1      b=3
```

### The ODE and Its Terms

```
-d/dz(κ·dT/dz) + vρC·dT/dz = Q(z)
```

| Term | Physical Meaning |
|------|------------------|
| `-d/dz(κ·dT/dz)` | **Diffusion**: Heat spreads from hot to cold regions |
| `vρC·dT/dz` | **Convection**: Flowing fluid carries heat downstream |
| `Q(z)` | **Source**: Electric coil adds heat |

### What Does v Control?

The parameter v (fluid velocity) determines which physical effect dominates:

- **v = 0**: Pure diffusion (no flow). Heat spreads symmetrically.
- **v > 0**: Convection present. Heat is carried in the flow direction.
- **Large v**: Convection dominates. Temperature profile becomes asymmetric.

This is why you're asked to compare different v values!

---

## Key Concept: BVP vs IVP

### Why Can't We Just "Step Forward"?

In Lab 1, you knew T(0) and marched forward in time. Here, you know T(0) AND T(L), but nothing in between.

If you tried stepping from z=0:
- To compute T at z=h, you'd need T at z=2h (the ODE has second derivatives)
- But you don't know T at z=2h yet!

### The Solution: Solve Everything at Once

Set up one equation for each interior grid point. All equations reference neighboring temperatures. Solve the coupled system simultaneously using linear algebra.

This leads to: **Au = b** where:
- **u** = vector of unknown temperatures at interior points
- **A** = matrix encoding how temperatures relate to neighbors
- **b** = vector with source terms and boundary condition contributions

---

## Finite Difference Discretization

### The Grid

Divide [0, L] into N+1 equal intervals:

```
z₀    z₁    z₂    ...   zₙ    zₙ₊₁
|-----|-----|-----|-----|-----|
0                             L
known  unknown interior     known
T₀     temperatures        T_out
```

**Grid spacing**: h = L / (N+1)

**Grid points**: zᵢ = i·h for i = 0, 1, 2, ..., N+1

### Approximating Derivatives

The ODE has first and second derivatives. Replace them with finite differences:

**First derivative** (how to choose?):
- Forward: `(Tᵢ₊₁ - Tᵢ) / h`
- Backward: `(Tᵢ - Tᵢ₋₁) / h`
- Central: `(Tᵢ₊₁ - Tᵢ₋₁) / (2h)`

**Second derivative** (standard central difference):
```
d²T/dz² ≈ (Tᵢ₊₁ - 2Tᵢ + Tᵢ₋₁) / h²
```

### Discretizing the Full ODE

Starting from:
```
-κ·d²T/dz² + vρC·dT/dz = Q(z)
```

At grid point zᵢ, substitute the finite difference approximations:
```
-κ·(Tᵢ₊₁ - 2Tᵢ + Tᵢ₋₁)/h² + vρC·(first derivative approx) = Q(zᵢ)
```

**Your task**: Rearrange this into the form:
```
(coeff for Tᵢ₋₁)·Tᵢ₋₁ + (coeff for Tᵢ)·Tᵢ + (coeff for Tᵢ₊₁)·Tᵢ₊₁ = Q(zᵢ)
```

These coefficients will form your tridiagonal matrix!

---

## Building the Tridiagonal System

### What is a Tridiagonal Matrix?

A matrix with nonzeros only on three diagonals:

```
[d₁  u₁  0   0   0 ]     [T₁]     [b₁]
[l₂  d₂  u₂  0   0 ]     [T₂]     [b₂]
[0   l₃  d₃  u₃  0 ]  *  [T₃]  =  [b₃]
[0   0   l₄  d₄  u₄]     [T₄]     [b₄]
[0   0   0   l₅  d₅]     [T₅]     [b₅]
```

- **d** = main diagonal (coefficient of Tᵢ)
- **u** = upper diagonal (coefficient of Tᵢ₊₁)
- **l** = lower diagonal (coefficient of Tᵢ₋₁)

### Why Tridiagonal?

Each equation only involves three consecutive temperatures: Tᵢ₋₁, Tᵢ, and Tᵢ₊₁. This comes directly from the finite difference stencil using nearest neighbors.

### Identifying Your Coefficients

From your discretized equation at point i, collect terms:

```
(something)·Tᵢ₋₁ + (something)·Tᵢ + (something)·Tᵢ₊₁ = Q(zᵢ)
     ↓                  ↓                ↓
    lᵢ                 dᵢ               uᵢ
```

**Hint**: The coefficients will involve combinations of κ, v, ρ, C, and h.

---

## Handling Boundary Conditions

### The Problem at Boundaries

At i=1 (first interior point), your equation references T₀.
At i=N (last interior point), your equation references Tₙ₊₁.

But T₀ and Tₙ₊₁ are **known** boundary values, not unknowns!

### The Solution: Move to Right-Hand Side

When a known value appears in your equation, move it to the right-hand side (the b vector).

**Example for i=1:**
```
l₁·T₀ + d₁·T₁ + u₁·T₂ = Q(z₁)
```

Since T₀ is known:
```
d₁·T₁ + u₁·T₂ = Q(z₁) - l₁·T₀
```

The term `-l₁·T₀` gets added to b₁.

**Similarly for i=N:** The term involving Tₙ₊₁ moves to bₙ.

### Structure of Your System

For N interior points:
- Matrix A is N×N (only unknowns)
- Vector u has N entries (T₁, T₂, ..., Tₙ)
- Vector b has N entries (source + boundary contributions at ends)

---

## The Source Term Q(z)

### The Piecewise Definition

```
Q(z) = 0                           if 0 ≤ z < a
Q(z) = Q₀·sin(π(z-a)/(b-a))        if a ≤ z ≤ b
Q(z) = 0                           if b < z ≤ L
```

With a=1, b=3: The coil is active only between z=1 and z=3, with a sinusoidal intensity profile (zero at edges, maximum in middle).

### Evaluating Q at Grid Points

For each grid point zᵢ, you need Q(zᵢ):

```matlab
% Pseudocode for computing Q at a single point z
if z < a
    Q = 0;
elseif z <= b
    Q = Q0 * sin(pi * (z - a) / (b - a));
else
    Q = 0;
end
```

### Building the Q Vector

You'll need Q evaluated at every interior grid point to form part of your b vector.

```matlab
% Pseudocode structure
Qvec = zeros(N, 1);
for i = 1:N
    z_i = i * h;
    Qvec(i) = ... % evaluate Q(z_i) using the piecewise formula
end
```

---

## MATLAB Implementation Patterns

### Setting Up Parameters

```matlab
% Given parameters
L = 10;
a = 1; b = 3;
Q0 = 50;
kappa = 0.5;
rho = 1;
C = 1;
T0 = 400;      % Left boundary
Tout = 300;    % Right boundary
v = 0;         % Start with v=0, then vary

% Discretization
N = 9;         % Number of interior points (try 9, 19, 39, 79)
h = L / (N + 1);
```

### Creating the Grid

```matlab
% All points including boundaries
z = 0:h:L;           % N+2 points total

% Interior points only (where we solve)
z_interior = h:h:L-h;  % or z(2:end-1)
```

### Building a Tridiagonal Matrix

**Option 1: Using `diag`**
```matlab
% Create vectors for each diagonal
main_diag = ... ;      % length N
upper_diag = ... ;     % length N-1
lower_diag = ... ;     % length N-1

A = diag(main_diag) + diag(upper_diag, 1) + diag(lower_diag, -1);
```

The second argument to `diag` specifies which diagonal:
- `diag(v, 0)` or `diag(v)`: main diagonal
- `diag(v, 1)`: one above main (upper)
- `diag(v, -1)`: one below main (lower)

**Option 2: Using sparse matrices (more efficient for large N)**
```matlab
A = sparse(N, N);
% Fill in entries...
```

### Building the Right-Hand Side Vector

```matlab
b = Qvec;  % Start with source term

% Add boundary contributions
b(1) = b(1) - (lower_coeff) * T0;
b(N) = b(N) - (upper_coeff) * Tout;
```

### Solving the System

```matlab
T_interior = A \ b;
```

### Assembling the Full Solution

```matlab
T_full = [T0; T_interior; Tout];  % Include boundary values
```

### Plotting

**Single curve:**
```matlab
plot(z, T_full)
xlabel('z')
ylabel('T(z)')
title('Temperature Distribution')
```

**Multiple N values on same plot:**
```matlab
hold on
% Plot for N=9, then N=19, etc.
% Use different colors or markers
hold off
legend('N=9', 'N=19', 'N=39', 'N=79')
```

**Subplots for different v values:**
```matlab
figure
subplot(2, 2, 1)
% Plot for v=0
title('v = 0')

subplot(2, 2, 2)
% Plot for v=0.1
title('v = 0.1')

subplot(2, 2, 3)
% Plot for v=0.5
title('v = 0.5')

subplot(2, 2, 4)
% Plot for v=1
title('v = 1')
```

---

## Interpreting Your Results

### Convergence Study (Varying N)

As N increases (more grid points, smaller h):
- The solution should become smoother
- Curves for different N should approach each other
- This demonstrates that your discretization is **converging** to the true solution

**What to look for:**
- Do the curves get closer together as N increases?
- Where are the biggest differences? (Often near the heat source or boundaries)

### Effect of Convection (Varying v)

**v = 0 (pure diffusion):**
- Heat spreads symmetrically from the source
- Temperature profile is relatively symmetric around the coil

**v > 0 (convection present):**
- Heat is "carried downstream" by the flow
- Temperature peak shifts in the flow direction
- Profile becomes asymmetric

**Large v (convection dominates):**
- Sharp temperature changes possible
- May see oscillations if discretization isn't appropriate (more on this below)

### Physical Sanity Checks

- Temperature should be between T₀=400 and T_out=300 at boundaries ✓
- Temperature should be elevated in/near the heating region
- No temperatures should be negative or unreasonably large

---

## Common Pitfalls

### Pitfall 1: Off-by-One Errors in Grid Indexing

MATLAB uses 1-based indexing, but your grid might be 0-based conceptually.

| Concept | MATLAB Index |
|---------|--------------|
| z₀ = 0 (boundary) | `z(1)` |
| z₁ (first interior) | `z(2)` |
| zₙ (last interior) | `z(N+1)` |
| zₙ₊₁ = L (boundary) | `z(N+2)` |

### Pitfall 2: Forgetting Boundary Contributions

If your solution looks wrong at the ends, check that you properly moved the boundary terms to the right-hand side.

### Pitfall 3: Sign Errors in Coefficients

Double-check signs when:
- Discretizing the negative sign in `-κ·d²T/dz²`
- Moving terms to the other side of the equation
- Computing the coefficients

### Pitfall 4: Wrong Diagonal Lengths

- Main diagonal: N elements
- Upper/lower diagonals: N-1 elements

If you get dimension mismatch errors, check your diagonal vector lengths.

### Pitfall 5: Inconsistent h

Make sure you use h consistently:
- `h = L / (N + 1)` if N is the number of **interior** points
- `h = L / N` if N is the number of **intervals**

The assignment says N+1 intervals with N interior points, so h = L/(N+1).

---

## Quick Reference: Problem Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| L | 10 | Pipe length |
| a | 1 | Coil start position |
| b | 3 | Coil end position |
| Q₀ | 50 | Heating intensity |
| κ | 0.5 | Heat conduction coefficient |
| ρ | 1 | Fluid density |
| C | 1 | Heat capacity |
| T₀ | 400 | Inlet temperature (z=0) |
| T_out | 300 | Outlet temperature (z=L) |
| v | 0, 0.1, 0.5, 1 | Fluid velocity (vary this) |
| N | 9, 19, 39, 79 | Interior grid points (vary this) |

---

## Summary: Implementation Checklist

1. [ ] Set up all parameters and compute h
2. [ ] Derive the discretized equation and identify coefficients
3. [ ] Build the tridiagonal matrix A
4. [ ] Compute Q(z) at all interior grid points
5. [ ] Build the right-hand side vector b (including boundary terms)
6. [ ] Solve A·T = b
7. [ ] Assemble full solution including boundary values
8. [ ] Plot and verify results make physical sense
9. [ ] Repeat for different N values (convergence study)
10. [ ] Repeat for different v values (convection study)

---

## Further Thinking (Not Required)

- **Why tridiagonal?** These systems can be solved in O(N) operations using the Thomas algorithm, making them very efficient even for large N.

- **Upwinding**: For large v, using a central difference for the first derivative can cause oscillations. "Upwind" schemes use backward differences for stability. This is a topic in convection-dominated problems.

- **Physical interpretation of κ and v**: The ratio vρC·h / κ (called the Péclet number when properly scaled) determines whether diffusion or convection dominates locally.

---

Good luck! The key insight is that BVPs require thinking about the **whole domain at once**, not stepping through sequentially.
