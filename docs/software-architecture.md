# Software architecture

## 1. Design principle

The solver is divided at mathematical interfaces rather than file-size boundaries:

```text
external geometry
      ↓
canonical mesh and topology
      ↓
materials + experiment specification
      ↓
compiled finite-element problem
      ↓
nonlinear/time solver
      ↓
observers and derived fields
```

The solver layer never interprets raw CAD names. The geometry layer never selects physical boundary conditions. Observers do not change the evolving state.

## 2. Canonical data objects

### Mesh

Created by `tdgl.mesh.fromArrays`, `tdgl.mesh.fromPdeMesh`, or a verification geometry generator. It contains:

- `nodes`: $N_v\times3$ coordinates;
- `cells`: positively oriented $N_t\times4$ tetrahedra;
- `regionIds`: material-volume identifiers;
- `topology.edges`, `topology.faces`, and cell-to-entity maps;
- orientation signs for local-to-global edge and face bases;
- face adjacency, boundary faces, and material interfaces;
- affine tetrahedron data and physical barycentric gradients;
- semantic boundary tags.

The canonical global edge orientation is smaller node ID to larger node ID. A global face is oriented by ascending node IDs. These rules make topology reproducible across runs and element reorderings.

### Material

`tdgl.materials.create` produces an explicit dimensionless record containing material kind, conductivity, permeability, GL activity, and the $a,b,K,u$ coefficients. Reference defaults exist for tests. Calibrated simulations must override them and record the SI scale conversion.

### Experiment

An experiment is an external struct containing:

- remote electromagnetic outer-boundary tags;
- boundary-condition records;
- terminal records;
- explicit interface models;
- time schedule and drive functions;
- initial conditions and requested observers.

`tdgl.problem.compile` resolves all named surfaces into face, edge, and node IDs. It rejects incomplete material assignments, overlapping terminals, unbalanced constant terminal currents, missing voltage references, incomplete exterior boundaries, and unstated GL-active interface laws.

### State

The final coupled state will contain:

```text
psi                 complex nodal order parameter
A                   real Nédélec edge potential
phi                 real nodal electric scalar potential
gaugeMultiplier     real nodal mixed-gauge variable
terminalMultipliers integral terminal constraints
cycleMultipliers    harmonic/fluxoid constraints
```

## 3. Geometry and topology layer

The mesh compiler does not merely hold tetrahedra. It constructs the discrete de Rham chain:

$$
\mathbb R^{N_v}\xrightarrow{G}
\mathbb R^{N_e}\xrightarrow{C}
\mathbb R^{N_f}\xrightarrow{D}
\mathbb R^{N_t}.
$$

The integer identities

$$
CG=0,\qquad DC=0
$$

are checked whenever a mesh is compiled. A failure is a fatal topology/orientation error.

General geometry enters through node/tetrahedron arrays supplied by a mesher. PDE Toolbox meshes can be converted with `tdgl.mesh.fromPdeMesh`. Boundary semantics are deliberately assigned by an external tagger because coordinate guessing is not safe for arbitrary devices.

## 4. Finite-element layer

The order parameter uses P1 Lagrange functions initially. The vector potential uses lowest-order first-family Nédélec functions

$$
\mathbf N_{ij}=\lambda_i\nabla\lambda_j-
\lambda_j\nabla\lambda_i.
$$

The global Nédélec coefficient is the oriented line integral

$$
A_e=\int_e\mathbf A\cdot d\mathbf l.
$$

This is why edge orientation is part of the mesh contract rather than an assembly detail.

`tdgl.assembly.covariantP1` assembles

$$
\int K\,
\overline{(\nabla-i\mathbf A)v}
\cdot(\nabla-i\mathbf A)u\,dV,
$$

using the Nédélec representation of $\mathbf A$. The nonlinear order-parameter routine assembles $b|\psi|^2\psi$ and its full real $2N_v\times2N_v$ Newton Jacobian.

## 5. Current solver kernels

### Magnetostatic kernel

`tdgl.solvers.solveMagnetostatic` solves

$$
\nabla\times(\mu^{-1}\nabla\times\mathbf A)=\mathbf J
$$

with tangential outer-boundary data and a mixed Coulomb constraint. The saddle system uses the edge mass matrix and discrete gradient to represent $(\mathbf A,\nabla q)=0$.

This kernel currently requires tangential data on the complete exterior boundary. Later variants will add scattered-field, source-coil, transient-conductivity, and FEM–BEM closures.

### Order-parameter kernel

`tdgl.solvers.stepOrderParameter` solves a fully implicit backward-Euler step for prescribed $\mathbf A$ and $\phi$:

$$
\frac{u}{\Delta t}\left(\psi^{n+1}-
e^{-i\phi^{n+1}\Delta t}\psi^n\right)
- (\nabla-i\mathbf A)^2\psi^{n+1}
-a\psi^{n+1}+b|\psi^{n+1}|^2\psi^{n+1}=0.
$$

The temporal link makes the nodal time-gauge transformation exact for the discrete step. A damped Newton method solves the real/imaginary block system. `tdgl.time.integrateOrderParameter` repeats the step and records observers.

This is a prescribed-field subsystem, not yet the final self-consistent TDGL–Maxwell solve.

## 6. Boundaries and terminals

Boundary records are physical declarations such as:

```matlab
bc = tdgl.problem.boundaryCondition( ...
    "outer field","A","tangential-dirichlet",outerTags,potential);
```

Terminals are similarly declarative:

```matlab
source = tdgl.problem.terminal("source","left","current",I);
drain  = tdgl.problem.terminal("drain","right","current",-I);
```

Current excitation is an integral surface constraint. A uniform pointwise current density will only be used when the experiment explicitly requests it.

An S–N or S–vacuum interface touching a GL-active material requires an explicit `tdgl.problem.interfaceModel`; the compiler does not silently choose transparent versus de Gennes coupling.

## 7. Voltage and time-dependent measurements

The scalar potential difference alone is not gauge invariant when $\mathbf A$ changes in time. For an oriented measurement path $C$ from $a$ to $b$, `tdgl.observe.pathVoltage` evaluates

$$
V_{ab}=\phi(a)-\phi(b)-
\frac{1}{\Delta t}\left[
\int_C\mathbf A^{n+1}\cdot d\mathbf l-
\int_C\mathbf A^n\cdot d\mathbf l
\right].
$$

The order-parameter electrochemical quantity is evaluated as

$$
\mu_{\mathrm{gi}}=\phi+\partial_t\arg\psi.
$$

Surface averages and terminal measurements are observers. In a time-varying magnetic field, a reported voltage must identify the terminal surfaces and lead/path convention.

## 8. Next coupled milestone

The next implementation checkpoint adds a monolithic residual/Jacobian with blocks for

$$
(\operatorname{Re}\psi,\operatorname{Im}\psi,
\mathbf A,\phi,p,\lambda_{\mathrm{terminal}}).
$$

It will add superconducting-current assembly, normal conductivity, the scalar-potential current-continuity equation, Maxwell time mass, terminal multipliers, and energy/work diagnostics. Only after its manufactured and conservation tests pass will proximity, large-domain topology, and adaptivity be added.
