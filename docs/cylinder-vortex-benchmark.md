# Finite superconducting cylinder in a remote homogeneous field

## Geometry and forcing

The benchmark uses a finite superconducting cylinder centered inside a larger
coaxial vacuum cylinder. Vacuum padding exists radially and above and below the
sample. The applied field is homogeneous and axial at the remote exterior:

$$
\mathbf B_a(t)=B_a(t)\,\mathbf e_z,
\qquad
\mathbf A_a=\tfrac12\mathbf B_a\times\mathbf r.
$$

Tangential $A_a$ is imposed only on the outer vacuum boundary. No magnetic
condition is placed on the superconductor surface. The computed field near the
sample is therefore free to become nonuniform through screening,
demagnetization, and vortex penetration.

## Protocol

1. Generate a fitted two-region tetrahedral mesh with region 1 as
   superconductor and region 2 as vacuum.
2. Initialize the superconducting nodes near $|\psi|=1$ with deterministic
   small noise and set $A=0$.
3. Ramp the remote homogeneous axial field to the configured maximum.
4. At every implicit step, solve the GL field and the complete sample-plus-vacuum
   MQS field self-consistently.
5. Save restart checkpoints independently of the lower-frequency trajectory.
6. Reconstruct $B$, $|\psi|$, and phase-winding vortex segments from any saved
   snapshot.

The default is an intentionally substantial production configuration. Mesh,
time-step, vacuum-padding, and random-seed studies are required before treating
an observed entry field as a physical result.

## Regression tests

`tests/unit/testCylinderBenchmark.m` checks:

- the presence of distinct sample and vacuum regions and a remote outer surface;
- exact reproduction of a uniform field in the same curved vacuum mesh when
  superconducting response is disabled;
- gauge-invariant detection of a known straight phase-winding vortex.

The full entry trajectory is a scheduled physical benchmark rather than a short
unit test. Its acceptance metrics are first-entry time/field, vortex count and
paths, minimum $|\psi|$, sample magnetization, current conservation, coupling
residual, and convergence under mesh, time-step, and vacuum-padding refinement.
