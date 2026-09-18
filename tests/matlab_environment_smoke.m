% Reproducible environment smoke test for planning the 3D TDGL solver.
% This is not production solver code.

fprintf('MATLAB_VERSION=%s\n', version);
fprintf('RELEASE=%s\n', version('-release'));
fprintf('COMPUTER=%s\n', computer);
fprintf('HOST_PROCESSORS=%s\n', getenv('NUMBER_OF_PROCESSORS'));
try
    fprintf('FEATURE_NUMCORES=%d\n', feature('numcores'));
catch exception
    fprintf('FEATURE_NUMCORES_ERROR=%s\n', exception.message);
end

fprintf('PDE_LICENSE=%d\n', license('test', 'PDE_Toolbox'));
fprintf('PCT_LICENSE=%d\n', license('test', 'Distrib_Computing_Toolbox'));
fprintf('OPT_LICENSE=%d\n', license('test', 'Optimization_Toolbox'));

try
    count = gpuDeviceCount;
    fprintf('GPU_COUNT=%d\n', count);
    if count > 0
        gpu = gpuDevice;
        fprintf('GPU_NAME=%s\n', gpu.Name);
        fprintf('GPU_COMPUTE_CAPABILITY=%s\n', gpu.ComputeCapability);
    end
catch exception
    fprintf('GPU_ERROR=%s\n', exception.message);
end

rng(7);
n = 1200;
A = sprandsym(n, 0.004, 0.2, 1) + 10*speye(n);
b = randn(n, 1);
x = A\b;
relativeResidual = norm(A*x - b)/norm(b);
fprintf('SPARSE_N=%d\n', n);
fprintf('SPARSE_NNZ=%d\n', nnz(A));
fprintf('SPARSE_RELRES=%.3e\n', relativeResidual);
assert(relativeResidual < 1e-10, 'Sparse linear solve residual is too large.');

model = createpde(1);
model.Geometry = multicuboid(1, 1, 1);
mesh = generateMesh(model, 'Hmax', 0.35, 'GeometricOrder', 'linear');
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 1);
applyBoundaryCondition(model, 'dirichlet', ...
    'Face', 1:model.Geometry.NumFaces, 'u', 0);
solution = solvepde(model);
u = solution.NodalSolution;
fprintf('TET_NODES=%d\n', size(mesh.Nodes, 2));
fprintf('TET_ELEMENTS=%d\n', size(mesh.Elements, 2));
fprintf('TET_UMIN=%.6g\n', min(u));
fprintf('TET_UMAX=%.6g\n', max(u));
assert(all(isfinite(u)), '3D tetrahedral solution contains nonfinite values.');
assert(max(u) > 0, '3D tetrahedral Poisson solution is not positive inside.');

fprintf('ALL_TESTS=PASS\n');
