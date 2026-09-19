%VACUUM_UNIFORM_FIELD Reproduce a uniform magnetic field in a vacuum box.
repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repositoryRoot);
startup;
mesh = tdgl.geometry.boxMesh([3 3 3],[-1 1;-1 1;-1 1]);
targetB = [0 0 0.2];
potential = tdgl.physics.uniformFieldPotential(targetB);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,potential);
solution = tdgl.solvers.solveMagnetostatic(mesh,struct( ...
    'muInv',1,'sourceCurrent',[0 0 0],'boundary',boundary));

fprintf('Maximum cellwise B error: %.3e\n', ...
    max(vecnorm(solution.cellMagneticFluxDensity-targetB,2,2)));
fprintf('Free residual: %.3e\n',solution.diagnostics.freeResidual);
