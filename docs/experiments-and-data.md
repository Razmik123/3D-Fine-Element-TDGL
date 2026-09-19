# Experiments, checkpoints, and saved data

## Run lifecycle

An experiment is defined by a versioned MATLAB configuration function. The
cylinder example is [`config/cylinder_vortex_entry_config.m`](../config/cylinder_vortex_entry_config.m).
It fixes geometry, material parameters, field ramp, time step, nonlinear
tolerances, execution policy, storage cadence, and random seed.

Start it from the repository root with:

```matlab
startup
configuration = cylinder_vortex_entry_config;
output = tdgl.experiments.runCylinderVortexEntry(configuration);
```

The convenience script `experiments/run_cylinder_vortex_entry.m` also creates a
final PNG. Production runs should normally be launched with `matlab -batch` so
the configuration and terminal log can be archived by a scheduler.

## Run directory format

Every run receives a unique directory under `results/<experiment>/<timestamp>/`:

```text
manifest.json          human-readable status, configuration, code and hardware
configuration.mat      exact MATLAB configuration, including function handles
mesh.mat               canonical oriented tetrahedral mesh and physical regions
trajectory.h5          append-only primary-field snapshots
checkpoint.mat         atomically replaced latest restart state
final-state.png        optional derived visualization
```

`trajectory.h5` stores time, global step number, real and imaginary parts of
$\psi$, Nédélec edge degrees of freedom for $A$, nodal $\phi$, and compact
coupling diagnostics. Datasets are chunked, compressed, and unlimited in time.
Writing `time` last is the snapshot commit marker: an interrupted partial write
is not counted as a complete snapshot.

Only primary state is mandatory. Magnetic field, electric field, currents,
energy, voltage, magnetization, and vortex segments are recomputed from the
stored mesh and consecutive primary states. This prevents large duplicated
arrays while retaining the information needed for new analyses.

## Snapshot cadence versus checkpoint cadence

These are independent:

- `StoreEvery` controls durable trajectory snapshots used for analysis.
- `CheckpointEvery` controls how many integration steps may be lost after an
  interruption. The checkpoint is written through a temporary file and then
  atomically replaces the previous checkpoint.

For a large run, a useful starting point is a checkpoint every 5–20 steps and a
trajectory snapshot every 20–100 steps. Fast vortex-entry or phase-slip windows
need a denser snapshot cadence.

## Restart

Set the saved directory and extend `StopTime`:

```matlab
configuration = cylinder_vortex_entry_config;
configuration.storage.ResumeDirectory = "results/cylinder-vortex-entry/...";
configuration.time.StopTime = 50;
output = tdgl.experiments.runCylinderVortexEntry(configuration);
```

The stored mesh is reused exactly; it is not regenerated. A SHA-256 mesh
fingerprint protects against combining incompatible degrees of freedom.

## Analysis without recalculation

Load one snapshot directly, without reading the full trajectory:

```matlab
snapshot = tdgl.io.readSnapshot(runDirectory,25);
analysis = tdgl.experiments.analyzeCylinderRun(runDirectory,25);
tdgl.visualization.plotCylinderSnapshot(analysis);
```

The analysis reconstructs cellwise $B$, $|\psi|$, pierced faces, and connected
within-cell vortex segments. More observers can be added later without rerunning
the PDE as long as they depend on stored primary fields and archived inputs.

## Provenance and retention

The manifest records the Git commit, solver version, MATLAB release, CPU/GPU
capabilities, mesh fingerprint, configuration, status, and snapshot count.
Generated run directories remain ignored by Git. Archive important completed
runs in institutional storage together with their terminal log; do not commit
large trajectories to the source repository.
