# 3-D finite-element TDGL–Maxwell

This repository is developing a modular MATLAB solver for three-dimensional time-dependent Ginzburg–Landau physics coupled to electromagnetic fields in superconductors, normal metals, insulators, and surrounding vacuum.

The production direction is:

- complex nodal finite elements for the order parameter $\psi$;
- Nédélec edge elements for the vector potential $\mathbf A$;
- nodal scalar potential $\phi$ and a mixed gauge constraint;
- an explicitly solved exterior vacuum rather than applied magnetic data on the superconducting surface;
- explicit interface, terminal, topology, and measurement definitions.

## Current implementation status

The repository currently contains the verified finite-element foundation, not yet the complete monolithic TDGL–MQS production solver.

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
- fixed-step time-series integration and observer recording;
- automated MATLAB unit tests.

Not yet implemented:

- monolithic feedback from superconducting/normal current into transient Maxwell;
- the production scalar-potential current-continuity block;
- normal-metal proximity assembly;
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
| `tdgl.boundary` | Boundary selection and finite-element interpolation |
| `tdgl.physics` | Applied/source field definitions |
| `tdgl.solvers` | Linear, mixed, and nonlinear field solves |
| `tdgl.time` | Time integration and observer scheduling |
| `tdgl.observe` | Gauge-aware voltage and other measurements |
| `tdgl.post` | Derived fields, energy, and later vortex diagnostics |

Detailed design and data flow are documented in [`docs/software-architecture.md`](docs/software-architecture.md). The governing equations are in [`docs/mathematical-model.md`](docs/mathematical-model.md).

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
