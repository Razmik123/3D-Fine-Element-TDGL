# Verification and validation plan

## 1. Policy

Verification asks whether the equations are solved correctly; validation asks whether the selected equations represent the experiment or accepted physical limit. Every benchmark must archive the configuration, mesh checksum, commit, MATLAB release, tolerances, nondimensional scales, and raw scalar metrics. A picture is supporting evidence, never the pass criterion.

Reference tolerances below are initial targets. They may be tightened as asymptotic behavior is established, but may not be relaxed merely to make a regression pass. All convergence fits use at least three asymptotic resolutions and report the fitted slope and confidence/fit residual.

## 2. Algebraic and element-level verification

### 2.1 Mesh topology and exact sequence

- Verify unique, globally oriented edges and faces under random element permutations.
- Verify signed incidence identities (CG=0) and (DC=0) exactly in integer arithmetic.
- Compare computed connected components and Betti numbers with cube, hollow cube, ring, torus, and disconnected meshes.
- Confirm material and boundary tags survive import/refinement.

**Pass:** exact incidence identities; invariant matrices/fields under cell reordering; known topology recovered.

### 2.2 Basis and quadrature

- Lagrange partition of unity, nodal interpolation, and polynomial reproduction.
- Nédélec edge moment interpolation, tangential continuity, mapped curl, and orientation sign.
- Raviart–Thomas/BDM normal continuity if enabled.
- Numerical quadrature against analytic monomial integrals and higher-order reference quadrature.

**Pass:** errors at roundoff for representable fields; expected quadrature order; no orientation-dependent result.

### 2.3 Jacobian and residual

Check every Newton block with complex-step differentiation where analytic, or centered finite differences over a step-size sweep. Test both complex and real-block representations.

**Pass:** directional derivative discrepancy reaches the expected roundoff/truncation minimum and is less than (10^{-7}) in normalized norm on representative states.

## 3. Manufactured-solution verification

### 3.1 Smooth coupled 3-D cube

Choose smooth nonzero complex $\psi(\mathbf x,t)$, vector potential $\mathbf A(\mathbf x,t)$, scalar potential $\varphi(\mathbf x,t)$, and gauge multiplier on the unit cube. Derive forcing, boundary data, initial data, and interface-free coefficients symbolically. Use a field with all components active and nonzero curl/divergence before gauge projection.

Measure:

- $\|\psi-\psi_h\|_{L^2}$ and $H^1$ seminorm;
- $\|\mathbf A-\mathbf A_h\|_{L^2}$ and $H(\mathrm{curl})$ norm;
- $\|\varphi-\varphi_h\|_{H^1}$;
- current, gauge, and nonlinear residuals.

For P1 Lagrange and first-order Nédélec on quasi-uniform meshes, require the theoretical first-order energy-norm behavior and the supported higher (L^2) behavior when regularity/duality permits. Rates, not a single tolerance, are decisive.

### 3.2 Discontinuous-material interface

Manufacture a solution in two subdomains with jumps in $\mu$, $\sigma$, and GL coefficients, satisfying the intended tangential/normal interface conditions. Use both fitted planar and curved interfaces.

**Pass:** optimal expected convergence with no degradation in interface jump residual; correct weighted flux continuity.

### 3.3 Nonconvex domain

Use an L-shaped prism or reentrant polyhedron with the known reduced Maxwell regularity. Compare the edge-element method with a nodal-vector reference only as a negative control.

**Pass:** edge/mixed formulation converges to a stable physical solution; refinement near the reentrant edge reduces error; no claim of a smooth-domain rate where regularity forbids it.

## 4. Temporal convergence and nonlinear stability

Use a spatial mesh fine enough that time error dominates and a smooth forced transient.

- Backward Euler: expected global order 1.
- BDF2: expected global order 2 after consistent startup.
- Crank–Nicolson/linearized CN: expected global order 2 in its documented stability regime.

Repeat with $\Delta t,\Delta t/2,\Delta t/4,\Delta t/8$, measuring field and observable error. Test phase-slip and vortex-entry events separately because nonsmooth event timing can reduce observed order.

**Pass:** fitted orders within 0.15 of theory for smooth tests; nonlinear tolerances at least one order below discretization error; rejected steps do not alter the converged trajectory beyond the estimator tolerance.

## 5. Gauge invariance and constraint tests

Select a smooth $\chi(\mathbf x,t)$ compatible with the boundary formulation and transform $\psi,\mathbf A,\varphi$. Run from both representations.

Compare $|\psi|$, $\mathbf B$, $\mathbf E$, currents, energy, voltage, force-like observables, and vortex lines.

Also verify:

- Coulomb-gauge residual decreases at the discretization rate;
- scalar and harmonic nullspaces are reported and constrained, not regularized accidentally;
- the temporal-gauge validation branch agrees in physical observables where both gauges apply.

**Pass:** physical differences scale with discretization/solver tolerance; no O(1) dependence on $\chi$ or mesh cut.

## 6. Current conservation

For closed and terminal-driven configurations, measure elementwise divergence residual, flux jumps, net boundary flux, and terminal balance:

$$
R_I=\frac{|I_++I_-|}{\max(|I_+|,|I_-|,I_{\rm scale})}.
$$

Integrate charge/current continuity over arbitrary unions of elements. Compare current computed from constitutive fields and from terminal multipliers.

**Pass:** global balance below (10^{-10}) for direct small solves or commensurate with iterative tolerance; local weak residual converges at the expected rate; no unexplained current through insulating boundaries.

## 7. Energy and work balance

### Undriven relaxation

With time-independent zero sources and compatible boundaries, the discrete free energy should be nonincreasing for the fully implicit baseline up to nonlinear/linear tolerance.

**Pass:** no step increases normalized energy by more than (10^{-10}) on small direct-solve tests; cumulative dissipation matches energy loss within temporal/discretization error.

### Driven system

With changing applied field or transport current, verify

$$
\Delta\mathcal G + \int \mathcal D\,dt = W_{\rm source}+\text{discretization error}.
$$

**Pass:** balance defect converges to zero with time and space refinement. Energy monotonicity is not incorrectly demanded when sources do work.

## 8. Physical benchmarks

### 8.1 Meissner screening in a slab

For a thick planar slab in the London/small-field regime, compare the interior profile with

$$
B(x)=B_a\frac{\cosh(x/\lambda)}{\cosh(d/(2\lambda))}
$$

for a slab centered at zero with thickness $d$. Use a surrounding vacuum and apply $B_a$ remotely.

**Metrics:** profile $L^2/L^\infty$ error, fitted penetration depth, surface current, magnetization, vacuum-padding dependence.

**Pass:** fitted $\lambda$ and magnetization converge to the London result; doubling vacuum padding changes target observables by less than the stated tolerance (initially 0.5%).

### 8.2 Superconducting sphere in vacuum

Place a sphere of radius $R$ in a much larger vacuum domain. In the London regime compare the magnetic moment/susceptibility with the analytical spherical screening expression, commonly written in normalized form proportional to

$$
-\frac32\left[1-3\frac{\lambda}{R}\coth\left(\frac{R}{\lambda}\right)
+3\left(\frac{\lambda}{R}\right)^2\right],
$$

with the exact SI normalization documented in the test implementation. Test several $\lambda/R$.

**Pass:** moment, axial field, and surface-current distribution converge with mesh and vacuum radius; extrapolation to infinite padding agrees with the analytical normalization within 1% in the resolved London regime.

### 8.3 Single straight vortex

Use a large cylinder or periodic transverse cell with one flux quantum. Compare radial $|\psi(r)|$ with a high-accuracy 1-D GL boundary-value solution (and Clem approximation only as a secondary reference), and far-core field with the London $K_0(r/\lambda)$ behavior.

**Pass:** phase winding $2\pi$; integrated flux approaches $\Phi_0$; core and field profiles converge; vortex energy per length approaches the reference with domain-size extrapolation.

### 8.4 Multiply connected ring/torus and fluxoid quantization

Initialize several winding sectors in a superconducting ring and torus. Integrate fluxoid on multiple homologous contours and use at least two independently generated cut systems.

**Pass:** fluxoid differs from $n\Phi_0$ only by discretization/quadrature error; homologous contours agree; results are independent of cut choice; $n$ remains fixed absent a resolved phase slip and changes by an integer when a phase slip occurs.

### 8.5 S–N proximity benchmark

For a planar S–N bilayer near $T_c$, compare the normal-side decay of $\psi$ with the linearized GL exponential and enforce the chosen transparent or finite-barrier interface conditions. Separately test the no-proximity model with $\psi$ absent in N but normal current present.

**Pass:** fitted decay length and interface amplitude/flux match the analytic coefficients; current is continuous; the two physical models remain distinguishable in input and output metadata.

### 8.6 Insulating cavity/pinning inclusion

Use a spherical cavity and compare vortex attachment/detachment trends and free-energy ordering with Doria et al. Repeat under mesh rotation.

**Pass:** qualitative state ordering reproduced; quantitative observables converge; cavity boundary carries no pair current and no spurious topology-dependent force.

### 8.7 Transport and voltage

Use a uniform normal conductor first, then an S–N or superconducting wire below critical current. Compare Ohm's law in the normal limit, zero/low voltage in a stationary superconducting state, and published current-driven vortex-motion trends.

**Pass:** imposed integral current recovered, terminal currents balance, voltage is gauge invariant, and steady line-integral voltage is path independent to tolerance.

## 9. Exterior-domain verification

For a fixed specimen and source, solve with geometrically increasing vacuum domains and at least two outer-boundary shapes. Compare fields near the specimen, moment, energy, and vortex state. Fit the truncation trend where possible.

When FEM–BEM is available, compare it to the vacuum-box extrapolation on identical interior meshes.

**Pass:** reported observables include an exterior error estimate; boundary shape does not change the extrapolated answer; FEM–BEM and box extrapolation agree within combined error bars.

For a later full-wave branch, use separate plane-wave scattering/cavity tests and verify PML or DtN truncation. These tests do not validate the MQS exterior, and vice versa.

## 10. Mesh adaptation verification

- Seed a moving curved vortex line and verify refinement follows core and screening layer.
- Compare adaptive versus uniform error at equal degrees of freedom.
- Refine/coarsen repeatedly and check gauge, current, energy, and fluxoid before/after transfer.
- Check corner/interface indicators against reference over-refined solutions.

**Pass:** adaptive meshes reduce selected observable error materially at fixed cost; transfer errors are below the requested time-step error; no change in winding sector without a physical phase slip.

## 11. Published numerical cross-checks

Reproduce small, well-specified cases from the collected literature:

- Hong–Ma–Xu reentrant-domain edge-element tests;
- Gao–Sun mixed-form convergence cases;
- Gunter–Kaper–Leaf implicit relaxation comparison;
- Winiecki–Adams current-carrying/self-field case at feasible scale;
- Sadovskyy et al. bulk periodic/pinning case with prescribed-field mode;
- Lara et al. qualitative 3-D vortex geometry, explicitly noting their different boundary physics;
- Alstrøm et al. 2-D complex-geometry state as a dimensional-reduction check.

Agreement is expected only after matching nondimensionalization, gauge-invariant observables, boundary conditions, and the authors' omission/inclusion of self-field. A mismatch caused by deliberately improved exterior physics should be explained, not tuned away.

## 12. Regression matrix

Every merge to the production branch should run:

| Tier | Content | Target runtime |
|---|---|---:|
| Smoke | assembly, one time step, sparse solve, tiny tetra mesh | < 2 min |
| Unit | topology, elements, Jacobian, IO, invariance | < 10 min |
| Numerical | manufactured spatial/temporal convergence | < 45 min |
| Physical | slab, sphere, vortex, S–N, ring | scheduled/nightly |
| Scale | million-DOF memory, iterative convergence, restart | release/HPC |

A benchmark failure must report the first violated physical metric and retain enough diagnostic state to reproduce it without committing large output files.
