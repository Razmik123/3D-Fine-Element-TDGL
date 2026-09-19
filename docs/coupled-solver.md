# Coupled TDGL–MQS reference solver

## Purpose and present scope

The first self-consistent implementation is an auditable staggered reference
solver. It advances the complex order parameter, vector potential, physical
electric scalar potential, and Coulomb-gauge multiplier. The electromagnetic
field is solved on the complete mesh, including surrounding vacuum. Applied
magnetic data therefore belong on the remote outer boundary, not the
superconductor surface.

This implementation is magnetoquasistatic. It omits displacement current and
does not model electromagnetic-wave propagation.

## One implicit time step

For fixed electromagnetic potentials, `stepOrderParameter` solves

$$
\frac{M_u}{\Delta t}
\left(\psi^{n+1}-e^{-i\phi^{n+1}\Delta t}\psi^n\right)
+K_A\psi^{n+1}-M_a\psi^{n+1}
+f_b(\psi^{n+1})=0.
$$

The exponential is the temporal gauge link. `K_A` is the covariant P1
stiffness matrix assembled from the Nédélec representation of the vector
potential. The nonlinear term and its full real/imaginary Jacobian use
fourth-order tetrahedral quadrature.

For fixed order parameter, `stepElectromagnetic` solves

$$
\left(\frac{M_\sigma}{\Delta t}+K_{\mu}+M_{K|\psi|^2}\right)A^{n+1}
+M_\sigma G\phi^{n+1}+M_1G p
=\frac{M_\sigma}{\Delta t}A^n+f_{\rm phase}+f_{\rm ext},
$$

where

$$
K_\mu=\kappa^2\,K_{\mathrm{curl},\mu^{-1}},
\qquad
f_{\rm phase}=\int K\,\operatorname{Im}(\overline\psi\nabla\psi)
\cdot N_e\,dV.
$$

The scalar-potential block is the discrete divergence of total current,

$$
G_c^T\left[
f_{\rm phase}-M_{K|\psi|^2}A^{n+1}
-M_\sigma\left(
\frac{A^{n+1}-A^n}{\Delta t}+G\phi^{n+1}
\right)+f_{\rm ext}
\right]=0.
$$

Only nodes belonging to conducting connected components carry scalar-potential
unknowns. One reference value per component is inserted when the caller has
not supplied one. The mixed Coulomb constraint is

$$
G_g^T M_1 A^{n+1}=0.
$$

`stepCoupledStaggered` alternates these two solves to a fixed-point tolerance.
Optional under-relaxation affects intermediate iterates only; the returned
state is an unrelaxed, boundary-consistent subproblem solution.

### Electrical terminals

`tdgl.boundary.scalarPotentialSpace` performs a nodal coordinate reduction.
Every current or floating electrode has one scalar-potential degree of freedom,
so it is equipotential without prescribing its unknown voltage. For a current
terminal, the aggregated continuity row enforces

$$
\int_{\Gamma_I}\mathbf J\cdot\mathbf n\,dS=I.
$$

The pointwise current profile is therefore an output of the solve. Voltage and
ground terminals are strong equipotential data. One current terminal may serve
as the additive potential reference; its current follows from component-wise
balance. Prescribed currents are checked separately on every disconnected
conducting component.

### No-proximity domains

An `ActiveCellMask`, or cellwise `model.glActive`, restricts every TDGL mass,
reaction, covariant-gradient, nonlinear, and superconducting-current term to
the selected fitted submesh. Nodes belonging exclusively to inactive cells are
held at zero. This implements the no-proximity baseline with its natural
zero-normal-covariant-flux interface. A finite de Gennes extrapolation length
still requires the planned interface surface term.

## State and model contract

Create a state with:

```matlab
state = tdgl.state.create(mesh,psi,A,phi);
```

The coupled material/source model contains:

```matlab
model = struct( ...
    'u',u,'a',a,'b',b,'K',K, ...
    'conductivity',sigma,'muInv',muInverse, ...
    'kappa',kappa,'sourceCurrent',sourceCurrent);
```

Coefficients may be scalar or cellwise arrays. `sourceCurrent` may be a
constant three-vector or a function evaluated at quadrature points.
`tdgl.problem.materialFields` expands compiled region materials into cellwise
arrays while keeping the GL-active mask explicit.

## Time series and measurements

`tdgl.time.integrateCoupled` accepts static structs or time-dependent model and
outer-boundary providers. Its observers receive the previous and current
state, which supports:

- edge electric field from $-\partial_t A-G\phi$;
- gauge-invariant path voltage from the electric-field line integral;
- superconducting terminal voltage from surface-averaged
  $\phi+\partial_t\arg\psi$;
- magnetic flux density and weak superconducting, normal, and total currents.

## Deliberate limitations of this checkpoint

The following are not silently approximated:

- de Gennes and finite-barrier interface terms in the assembled residual;
- harmonic/cohomology constraints for multiply connected production meshes;
- monolithic Newton coupling and scalable block preconditioning;
- adaptive remeshing and FEM–BEM exterior closure.

Until those modules pass their gates, the staggered solver is a reference
implementation for simply connected, fitted meshes and explicit
coefficient-based GL regions.
