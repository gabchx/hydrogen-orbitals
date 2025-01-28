### **Code Explanation**

This code generates 3D visualizations of hydrogen-like atomic orbitals using Python. The orbitals are represented as probability density clouds, showing the spatial distribution of electrons in the atom. The visualization includes 3D scatter plots and 2D cross-sections of the electron cloud.

#### **1. Key Components**

- **Radial Wavefunction (`radial`):**

  - Computes the radial part of the hydrogen-like orbital using associated Laguerre polynomials.
  - Normalized using quantum numbers `n` (principal) and `l` (azimuthal).
  - Formula derived from the Schrödinger equation for hydrogen-like atoms.

- **Spherical Harmonics (`sph_harm`):**

  - Represents the angular part of the wavefunction using spherical harmonics.
  - `sph_harm(m, l, phi, theta)` matches the quantum mechanical convention for angular momentum.

- **Wavefunction in Cartesian Coordinates (`wave_func_cart`):**

  - Converts Cartesian coordinates `(x, y, z)` to spherical coordinates `(r, θ, φ)`.
  - Combines radial and angular components to compute the full wavefunction.

- **Orbital Class (`Orbital_3D`):**

  - Initializes a 3D grid of points.
  - Computes the wavefunction and probability density (`prob`).
  - Normalizes the probability to ensure total probability = 1.

- **Transition Class (`Transition_3D`):**
  - Manages animation events (transitions between orbitals or pauses).
  - Interpolates wavefunctions during transitions using a time-dependent linear combination:
    ```python
    psi = c1 * orb1.psi + c2 * orb2.psi  # c1 = cos(angle), c2 = sin(angle)
    ```
  - Renders frames by sampling 3D points according to the probability distribution and plotting them as a scatter plot.

#### **2. Animation Workflow**

1. **Frame Generation:**

   - For each frame, the code checks if it's a `transition` (interpolating between two orbitals) or a `wait` (static orbital).
   - During transitions, the probability density smoothly evolves from the initial to the target orbital.
   - Points are randomly sampled from the probability distribution to visualize the orbital's electron cloud.

2. **Visualization:**
   - **3D Scatter Plot:** Shows the electron density cloud.
   - **2D Cross-Sections:** Display slices through the `xy` and `xz` planes using contour plots.
   - Quantum numbers `(n, l, m)` are displayed dynamically.

## Sources :

- [Quantum Orbitals](https://en.wikipedia.org/wiki/Atomic_orbital)
- [Spherical Harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics)
- [Laguerre Polynomials](https://en.wikipedia.org/wiki/Laguerre_polynomials)
- [Hydrogen Atom Wavefunctions](https://en.wikipedia.org/wiki/Hydrogen_atom)
- [Quantum Mechanics with Python: Hydrogen Wavefunctions and Electron Density Plots](https://ssebastianmag.medium.com/computational-physics-with-python-hydrogen-wavefunctions-electron-density-plots-8fede44b7b12)
- [qmsolve](https://github.com/quantum-visualizations/qmsolve/tree/main)
