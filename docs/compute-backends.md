# CPU and GPU execution

## Policy

Execution is explicit:

```matlab
cpu  = tdgl.compute.execution("cpu");
gpu  = tdgl.compute.execution("gpu",'AllowFallback',false);
auto = tdgl.compute.execution("auto", ...
    'MinimumUnknowns',15000,'MemoryFraction',0.3);
```

The policy is passed through the coupled solver to every large sparse linear
solve. Explicit GPU mode transfers the sparse matrix and right-hand side to the
selected CUDA device, solves there, and gathers the solution. Automatic mode
uses the GPU only when:

- a supported double-precision device exists;
- the system exceeds the configured crossover size;
- a conservative sparse-factor memory estimate fits the allowed fraction of
  currently available device memory.

If an automatic GPU solve fails because of factor fill or an unsupported matrix
case, the default policy falls back to CPU and records that fact in diagnostics.
Set `AllowFallback=false` when benchmarking so a requested GPU run cannot
silently become a CPU run.

## What is accelerated now

The sparse magnetostatic, electromagnetic saddle, and nonlinear Newton linear
systems support CPU or GPU execution. Mesh construction and element assembly
remain on CPU in this reference implementation. This is deliberate: transferring
thousands of tiny tetrahedral matrices individually is slower than CPU assembly.
A future high-throughput assembly kernel should batch elements and benchmark the
complete step, not only a device kernel.

`tdgl.compute.environment` reports the runtime capabilities. On the development
machine the verified device is an NVIDIA GeForce RTX 3060 Laptop GPU with CUDA
compute capability 8.6, double precision, and approximately 6 GB physical
memory. The CPU remains preferable for small verification cases because GPU
transfer and factorization setup dominate.

## Reproducibility

CPU and GPU answers are tested by residual, not bitwise equality. Sparse
factorization order and floating-point reduction order can differ. Scientific
regressions compare gauge-invariant observables using tolerances below the
discretization error.
