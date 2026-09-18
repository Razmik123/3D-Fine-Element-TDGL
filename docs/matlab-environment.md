# MATLAB environment verification

## Verification date and scope

Environment checked on 2026-09-19 in batch mode from the repository root. The checked-in test is [`tests/matlab_environment_smoke.m`](../tests/matlab_environment_smoke.m). It is an environment/finite-element smoke test, not production TDGL code.

MATLAB was given a repository-local temporary preference directory so the test would not depend on or modify the user's normal MATLAB preferences. That directory and other generated MATLAB state are excluded by `.gitignore`.

## Installed MATLAB and relevant licenses

| Component | Detected version/status |
|---|---|
| MATLAB | 9.13.0.2049777 (R2022b), 64-bit Windows (`PCWIN64`) |
| Partial Differential Equation Toolbox | 3.9; license test passed |
| Parallel Computing Toolbox | 7.7; license test passed |
| Optimization Toolbox | 9.4; license test passed |

These are the licenses directly relevant to the proposed workflow. A license test reports present/available at test time; it does not guarantee that a floating license will always be free on another machine or date.

## CPU and GPU

| Resource | Batch result |
|---|---|
| Host logical processors exposed by environment | 12 |
| MATLAB `feature('numcores')` | 6 |
| GPUs visible to MATLAB | 1 |
| GPU | NVIDIA GeForce RTX 3060 Laptop GPU |
| CUDA compute capability | 8.6 |

The operating-system query for the precise CPU model was denied by host permissions, so no model name is claimed. The MATLAB-reported core availability is sufficient for planning. GPU visibility does not imply that sparse complex edge-element assembly or factorization will benefit; GPU kernels should be adopted only after profiling representative block operations.

## Batch command

The successful run used the installed `matlab` executable in noninteractive batch mode, with `MATLAB_PREFDIR` set to the ignored project-local `.matlab-pref` directory, and evaluated:

```matlab
run('tests/matlab_environment_smoke.m')
```

On first launch MATLAB warned that it could not create the default personal work folder `C:\Users\Razmik\Documents\MATLAB`. The project-local preference setup avoided dependence on that folder, and all test assertions passed.

## Sparse linear-system test

The script constructs a reproducible sparse symmetric matrix of size 1200 with a positive diagonal shift, solves with MATLAB's sparse backslash, and checks the normalized residual.

| Metric | Result |
|---|---:|
| Matrix dimension | 1200 |
| Nonzeros | 5,480 |
| Relative residual $\|Ax-b\|/\|b\|$ | $1.734\times10^{-16}$ |
| Status | PASS |

This verifies the basic sparse direct-solve path. It does not test the indefinite mixed saddle-point structure or large-scale iterative preconditioners required by production TDGL–Maxwell.

## 3-D tetrahedral finite-element smoke test

The script uses PDE Toolbox to:

1. create a unit cube;
2. generate a linear tetrahedral mesh with `Hmax = 0.35`;
3. solve $-\nabla^2u=1$ with homogeneous Dirichlet data on all faces;
4. assert finite values and a positive interior maximum.

| Metric | Result |
|---|---:|
| Mesh nodes | 90 |
| Tetrahedral elements | 243 |
| Minimum solution | 0 |
| Maximum solution | 0.0506315 |
| Status | PASS |

This verifies licensed 3-D geometry, tetrahedral meshing, scalar PDE assembly, and solution. It does **not** establish that PDE Toolbox has Nédélec edge elements or the mixed spaces needed for the production solver; the assessment in [`implementation-plan.md`](implementation-plan.md) concludes that custom exact-sequence assembly is required.

## Reproduction and expected output

The final output record was:

```text
MATLAB_VERSION=9.13.0.2049777 (R2022b)
RELEASE=2022b
COMPUTER=PCWIN64
HOST_PROCESSORS=12
FEATURE_NUMCORES=6
PDE_LICENSE=1
PCT_LICENSE=1
OPT_LICENSE=1
GPU_COUNT=1
GPU_NAME=NVIDIA GeForce RTX 3060 Laptop GPU
GPU_COMPUTE_CAPABILITY=8.6
SPARSE_N=1200
SPARSE_NNZ=5480
SPARSE_RELRES=1.734e-16
TET_NODES=90
TET_ELEMENTS=243
TET_UMIN=0
TET_UMAX=0.0506315
ALL_TESTS=PASS
```

## Environment recommendation

Keep R2022b as the reproducibility baseline until the element/topology unit tests exist. Before production implementation, also test the chosen external tetrahedral mesher, confirm Parallel Computing Toolbox worker startup, measure sparse complex/real-block memory, and establish whether the target workstation or cluster has a supported iterative/preconditioning path. Do not make the production architecture depend on GPU execution.
