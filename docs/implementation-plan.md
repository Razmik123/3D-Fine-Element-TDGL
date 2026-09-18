# MATLAB implementation plan

## 1. Executive recommendation

Use MATLAB PDE Toolbox for geometry import, tetrahedral mesh prototyping, scalar reference solves, and visualization, but **not** as the production discretization engine. The production solver requires custom tetrahedral exact-sequence assembly:

- complex Lagrange elements for the order parameter;
- Nédélec edge elements for the vector potential;
- nodal scalar-potential and gauge spaces;
- optional Raviart–Thomas/BDM face elements for conservative flux/current variables;
- explicit topology, cycle, and harmonic-field handling.

The first end-to-end physics target is magnetoquasistatic TDGL coupled to a surrounding truncated vacuum domain. FEM–BEM exterior coupling and full electromagnetic-wave coupling are later, separate milestones.

No production solver is implemented in this research checkpoint.

## 2. PDE Toolbox suitability assessment

The tested local version is MATLAB R2022b with Partial Differential Equation Toolbox 3.9. The [documented general PDE form](https://www.mathworks.com/help/pde/ug/equations-you-can-solve.html) is a scalar/system divergence-form workflow with nodal unknowns. The [documented mesh data](https://www.mathworks.com/help/pde/ug/mesh-data.html) consists of linear or quadratic triangles/tetrahedra and does not support mixed cell types. The [geometry and mesh workflow](https://www.mathworks.com/help/pde/geometry-and-mesh.html) is nevertheless useful for import, region labeling, 3-D tetrahedral mesh prototyping, interpolation, visualization, and scalar elliptic/parabolic reference solves.

It is insufficient as the core production formulation for these reasons:

1. **Maxwell conformity.** A vector represented by three nodal components is not an (H(\mathrm{curl})\)-conforming Nédélec field. It can generate spurious modes and wrong limits on reentrant/nonconvex domains where \(\mathbf A\notin H^1\).
2. **Mixed de Rham spaces.** The required nodal–edge–face–cell exact sequence and mixed block operators are not exposed as native PDE Toolbox trial spaces.
3. **Multiply connected topology.** Gauge fixing, harmonic fields, cohomology bases, fluxoid sectors, and cuts require explicit mesh-topology operators.
4. **3-D adaptivity.** The documented legacy [`adaptmesh`](https://www.mathworks.com/help/pde/ug/adaptmesh.html) workflow is two-dimensional; production 3-D vortex-core refinement and coarsening require an external/custom mesh loop.
5. **Exterior Maxwell domain.** PDE Toolbox does not by itself provide the required MQS FEM–BEM coupling or a topology-aware unbounded-domain closure.
6. **Block solvers.** Robust Newton–Krylov and field-split preconditioners require access to all coupled blocks and nullspaces.

For convex, simply connected toy problems, a PDE Toolbox nodal prototype can be a useful independent comparison. It must not be promoted to the main solver based only on matching such smooth tests.

## 3. Proposed repository structure

```text
3D_Finite_elemnt_TDGL/
├── README.md
├── papers/                       # legal source PDFs, manifest, BibTeX
├── docs/                         # model, literature, implementation, verification
├── config/                       # versioned YAML/JSON/MATLAB run specifications
├── geometry/                     # parameterized CAD/CSG descriptions
├── meshes/                       # small canonical meshes and metadata only
├── src/
│   └── +tdgl/
│       ├── +mesh/                # import, orientation, tags, refinement transfer
│       ├── +topology/            # incidence, cycles, cohomology, cuts
│       ├── +elements/            # Lagrange, Nédélec, RT/BDM bases/quadrature
│       ├── +assembly/            # local kernels and sparse global assembly
│       ├── +physics/             # material laws, TDGL, Maxwell, interfaces
│       ├── +solvers/             # Newton/Krylov, nullspaces, preconditioners
│       ├── +time/                # BE, BDF2, CN, adaptivity and events
│       ├── +post/                # fields, energy, vortices, fluxoid, voltage
│       └── +io/                  # configuration, checkpoints, provenance
├── tests/
│   ├── unit/                     # element, topology, quadrature tests
│   ├── manufactured/             # forced exact solutions
│   ├── benchmarks/               # physical and published benchmarks
│   └── matlab_environment_smoke.m
├── examples/                     # small, reproducible documented runs
├── third_party/                  # license-recorded vendored code, if approved
└── results/                      # generated output; ignored except tiny gold data
```

Large meshes, transient fields, MATLAB preferences, caches, and duplicate paper files must remain outside Git or under ignored generated directories. Every result file should record commit hash, configuration, material scales, mesh checksum, MATLAB release, solver tolerances, and random seed.

## 4. Architectural decisions

### Mesh and topology

- Accept tetrahedral meshes from PDE Toolbox initially; evaluate Gmsh for robust physical groups, boundary layers, and local 3-D adaptation.
- Store vertex, global edge, face, and cell orientation once. Element mappings must preserve orientation signs.
- Build integer incidence matrices (G,C,D) and unit-test (CG=0), (DC=0).
- Label every cell and boundary face by material/interface/terminal role; never infer physics from coordinate tolerances during assembly.
- Compute connected components, first Betti number, cycle basis, and harmonic basis. Treat changes in topology during mesh conversion as fatal.

### Assembly

- Derive weak forms from [`mathematical-model.md`](mathematical-model.md) and keep dimensional scaling in one immutable object.
- Use vectorized element batches to build sparse triplets; avoid repeated sparse insertion.
- Provide independent high-order quadrature checks for nonlinear terms.
- Split applied and scattered potentials where this improves outer boundary data and conditioning.
- Make complex algebra explicit. If MATLAB solver/preconditioner limitations require it, expose a real (2\times2) block representation with equivalence tests.

### Nonlinear and linear solvers

- Baseline: monolithic damped Newton for backward Euler, with residual-based line search and a Picard fallback for poor initial guesses.
- Exploit block structure: order parameter, edge potential, electric potential, gauge multiplier, terminal/cycle multipliers.
- Prototype with MATLAB sparse direct factorization on small meshes. Scale with GMRES/MINRES as appropriate and physics-based block preconditioners.
- Project or constrain all known nullspaces. Never rely on a tiny diagonal perturbation as the final gauge treatment.
- Report residuals by physical block and constraint, not only a global norm.

### Time integration

- First: backward Euler, fully implicit, adaptive step rejection based on nonlinear convergence and embedded/step-doubling error estimates.
- Second: BDF2 with backward-Euler startup.
- Third: linearized Crank–Nicolson or IMEX only after an energy/work-balance test shows the intended stability.
- Limit time-step growth near phase slips, vortex entry, rapidly changing drive, and mesh adaptation.

### Mesh resolution and adaptation

Initial engineering targets, to be replaced by convergence evidence, are

- (h\le\xi/4) in vortex-core regions, giving roughly 8 cells across a (2\xi) core diameter;
- (h\le\lambda/4) in penetration-depth boundary layers;
- additional refinement at S–N interfaces, sharp corners, holes, terminals, and material jumps;
- graded coarsening through remote vacuum, constrained by magnetic far-field accuracy.

Use residual indicators with separate TDGL, curl-curl, current-continuity, interface, and goal-oriented components. Mesh transfer must preserve complex phase as well as possible, reproject the gauge, and preserve loop winding/fluxoid. A refined mesh is accepted only if incidence identities and material tags survive.

## 5. Phased delivery plan and gates

### Phase 0 — Research and reproducibility foundation

**Deliverables:** legal paper library, bibliography, mathematical specification, environment record, verification plan, repository policy.

**Gate:** all documents reviewed; no applied field placed directly on the finite superconductor boundary in the target model; MQS and full-wave formulations clearly distinguished.

### Phase 1 — Tetrahedral exact-sequence foundation

Implement mesh import, global orientations, quadrature, P1/P2 Lagrange and first-order Nédélec basis functions, incidence matrices, interpolation, and basic visualization.

**Gate:** polynomial reproduction; commuting-diagram tests; (CG=0), (DC=0) exactly; element orientation invariance; nonconvex and toroidal meshes retain topology.

### Phase 2 — Standalone exterior Maxwell/MQS solver

Implement curl-curl, conductivity mass, scalar electric potential, Coulomb gauge, source currents, remote boundary conditions, material transmission, and vacuum padding.

**Gate:** manufactured (H(\mathrm{curl})) convergence, current conservation, gauge invariance, sphere/slab screening limits, vacuum-size convergence, and no spurious null modes.

### Phase 3 — Standalone TDGL on prescribed compatible fields

Implement complex \(\psi\), covariant gradient, de Gennes conditions, nonlinear reaction, backward Euler, Newton/Picard, energy and vortex diagnostics.

**Gate:** manufactured convergence, zero-field relaxation, single-vortex flux/winding, published bounded-domain comparisons, and energy decay without drive.

### Phase 4 — Monolithic coupled TDGL–MQS in vacuum

Couple all Jacobian blocks, retain physical \(\varphi\), add applied-field/source splitting, restart files, and field-split preconditioning.

**Gate:** self-consistent Meissner screening, sphere in vacuum, vortex entry, energy/work balance, gauge-transform invariance, and mesh/time/vacuum convergence.

### Phase 5 — Materials, interfaces, and transport

Add normal-metal conductivity, no-proximity and GL-proximity alternatives, insulating cavities, terminals, integral-current constraints, and voltage extraction.

**Gate:** planar S–N decay/interface benchmark, global/local current balance, transport voltage path consistency, and insulating-inclusion comparison.

### Phase 6 — Multiply connected superconductors

Add cycle/cohomology basis generation, harmonic constraints, winding-sector initialization, fluxoid tracking, and phase-slip detection.

**Gate:** ring/torus fluxoid quantization, independence from mesh cuts, preservation of winding absent phase slip, and correct sector change during a resolved phase slip.

### Phase 7 — 3-D adaptivity and scalability

Add residual/goal indicators, refine/coarsen loop, conservative state transfer, parallel element assembly, iterative solvers, profiling, and GPU kernels only where measured beneficial.

**Gate:** error versus degrees-of-freedom improvement, invariant preservation after remesh, memory/performance budgets, and reproducible restart.

### Phase 8 — Unbounded and full-wave extensions

First implement MQS FEM–BEM exterior coupling and compare with expanding vacuum boxes. Only then, as a separate physics module, add displacement current with PML/DtN/radiation conditions.

**Gate:** exterior-method cross-validation; for full wave, cavity/scattering benchmarks and a documented regime in which MQS and full wave agree.

### Phase 9 — Validation release

Freeze reference configurations and small gold datasets, write user documentation, automate the full verification matrix, and archive reproducibility metadata.

**Gate:** every acceptance criterion in [`verification-plan.md`](verification-plan.md) passes on a clean checkout.

## 6. Proposed checkpoint sequence

The following commits keep reviewable scientific boundaries:

1. `docs: add curated TDGL research library`
2. `docs: specify 3D TDGL Maxwell formulation`
3. `docs: add MATLAB implementation and verification plans`
4. `feat(mesh): add oriented tetrahedral topology`
5. `feat(fem): add Lagrange and Nedelec element kernels`
6. `feat(maxwell): add mixed MQS vacuum solver`
7. `feat(tdgl): add implicit order-parameter evolution`
8. `feat(coupling): couple TDGL to exterior MQS fields`
9. `feat(materials): add interfaces proximity and transport`
10. `feat(topology): add fluxoid and harmonic-cycle constraints`
11. `feat(adapt): add 3D residual mesh adaptation`
12. `test(validation): freeze physical benchmark suite`

At every checkpoint: run unit and convergence tests, record MATLAB version, review generated-file status, and push only after the tree is clean. Never commit local preference files, crash dumps, generated caches, full transient output, or duplicate PDFs.

## 7. Principal risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| Nodal vector potential on nonsmooth geometry | spurious or wrong solution | Nédélec (H(\mathrm{curl})) production space |
| Gauge/harmonic nullspaces | singular or topology-dependent solve | mixed constraint plus explicit cohomology basis |
| Applied field on sample surface | suppressed screening/demagnetization | solve surrounding vacuum; apply data remotely or use coils |
| Under-resolved \(\xi\) or \(\lambda\) | wrong vortex force and entry field | local refinement plus observable convergence |
| Uncontrolled vacuum truncation | biased magnetization/self-field | padding study, then FEM–BEM |
| Semi-implicit instability | unphysical energy increase | fully implicit baseline and energy/work gate |
| GL proximity used outside validity | misleading S–N prediction | label regime; compare with analytic GL limit; defer microscopic model |
| MATLAB memory/block fill | scale ceiling | matrix-free kernels, field splits, profiling, external mesh/BEM libraries only with license review |

## 8. Immediate next implementation task

After scientific review of these documents, Phase 1 should begin with a tiny hand-checkable two-tetrahedron mesh. The first code checkpoint should contain only topology, orientations, basis evaluation, quadrature, and exact-sequence unit tests—not coupled physics. This isolates the most consequential 3-D finite-element conventions before nonlinear TDGL complexity is introduced.
