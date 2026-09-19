# 3-D finite-element TDGL–Maxwell

This repository is developing a modular MATLAB solver for three-dimensional time-dependent Ginzburg–Landau physics coupled to electromagnetic fields in superconductors, normal metals, insulators, and surrounding vacuum.

The production direction is:

- complex nodal finite elements for the order parameter $\psi$;
- Nédélec edge elements for the vector potential $\mathbf A$;
- nodal scalar potential $\phi$ and a mixed gauge constraint;
- an explicitly solved exterior vacuum rather than applied magnetic data on the superconducting surface;
- explicit interface, terminal, topology, and measurement definitions.

## Current implementation status

The repository now contains a verified finite-element foundation and a
self-consistent staggered TDGL–MQS reference solver. It is not yet the final
monolithic, topology-aware production solver.

Implemented:

- canonical tetrahedral mesh compiler with positive orientation;
- globally oriented vertices, edges, faces, and cells;
- exact incidence operators satisfying $CG=0$ and $DC=0$;
- P1 Lagrange and lowest-order first-family Nédélec bases;
- sparse nodal mass/stiffness and edge mass/curl-curl assembly;
- mixed Coulomb-constrained magnetostatic vacuum solve;
- declarative materials, interfaces, boundary conditions, and terminals;
- gauge-aware path voltage and electrochemical-potential observers;
- fully implicit nonlinear TDGL order-parameter stepping for prescribed $\mathbf A$ and $\phi$;
- transient MQS edge/nodal solve with superconducting screening, Ohmic current,
  current continuity, and a mixed Coulomb constraint;
- self-consistent staggered TDGL–MQS time stepping with convergence diagnostics;
- equipotential current, voltage, ground, and floating-terminal constraints,
  including prescribed total rather than pointwise current;
- restricted GL-active cell sets for the no-proximity material model;
- finite-cylinder/vacuum geometry and a remote homogeneous-field vortex-entry
  experiment;
- append-only compressed HDF5 trajectories, atomic restart checkpoints, mesh
  fingerprints, and provenance manifests;
- explicit CPU, GPU, and automatic sparse-solve execution policies;
- fixed-step time-series integration and observer recording;
- gauge-invariant path voltage, terminal electrochemical voltage, and edge
  electric-field extraction;
- automated MATLAB unit tests.

Not yet implemented:

- monolithic Newton coupling and scalable block preconditioning;
- assembled finite-barrier/de Gennes interface terms beyond the natural
  zero-flux no-proximity boundary;
- scalable cohomology basis for large multiply connected meshes;
- adaptive 3-D remeshing and FEM–BEM exterior coupling.

This distinction is intentional: tests must validate each foundation before it is used in the coupled nonlinear solver.

## Quick start

From MATLAB in the repository root:

```matlab
startup
information = tdgl.version
```

Run all tests:

```matlab
addpath tests
results = run_all_tests;
table(results)
```

Run the first examples:

```matlab
run examples/vacuum_uniform_field.m
run examples/uniform_order_parameter_relaxation.m
run examples/coupled_uniform_relaxation.m
```

Run the cylinder experiment:

```matlab
run experiments/run_cylinder_vortex_entry.m
```

## Source layout

| Package | Responsibility |
|---|---|
| `tdgl.geometry` | Verification geometry generators |
| `tdgl.mesh` | Mesh import, compilation, validation, and physical tags |
| `tdgl.topology` | Incidence operators, paths, and topology diagnostics |
| `tdgl.elements` | Reference/physical finite-element basis and quadrature |
| `tdgl.assembly` | Sparse finite-element operators and nonlinear residuals |
| `tdgl.materials` | Explicit material records |
| `tdgl.problem` | Experiment, interface, boundary, and terminal compilation |
| `tdgl.state` | Validated coupled field-state construction |
| `tdgl.boundary` | Boundary selection and finite-element interpolation |
| `tdgl.physics` | Applied/source field definitions |
| `tdgl.solvers` | Linear, mixed, and nonlinear field solves |
| `tdgl.time` | Time integration and observer scheduling |
| `tdgl.observe` | Gauge-aware voltage and other measurements |
| `tdgl.post` | Derived fields, energy, and later vortex diagnostics |
| `tdgl.io` | Streaming trajectories, restart checkpoints, provenance |
| `tdgl.compute` | CPU/GPU selection and linear solves |
| `tdgl.experiments` | Reproducible run and offline-analysis workflows |

Detailed design and data flow are documented in
[`docs/software-architecture.md`](docs/software-architecture.md). The governing
equations are in [`docs/mathematical-model.md`](docs/mathematical-model.md), and
the implemented coupled discretization is described in
[`docs/coupled-solver.md`](docs/coupled-solver.md).
Experiment storage is documented in
[`docs/experiments-and-data.md`](docs/experiments-and-data.md), GPU/CPU behavior
in [`docs/compute-backends.md`](docs/compute-backends.md), and the cylinder case
in [`docs/cylinder-vortex-benchmark.md`](docs/cylinder-vortex-benchmark.md).

## Scientific invariants

The implementation must preserve or explicitly test:

- gauge-independent physical observables;
- exact discrete curl of a gradient and divergence of a curl;
- current conservation;
- correct material-interface traces;
- flux and fluxoid constraints in multiply connected domains;
- energy decay for undriven implicit relaxation;
- exterior-domain convergence for finite vacuum truncations.

No module may silently infer a proximity law, terminal excitation, magnetic outer boundary, or winding sector.
