function tests = testMagnetostaticSolver
%TESTMAGNETOSTATICSOLVER Verify a uniform field in a vacuum box.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testUniformMagneticFieldReproduction(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2],[-1 1;-1 1;-1 1]);
targetB = [0.17 -0.08 0.23];
potential = tdgl.physics.uniformFieldPotential(targetB,[0 0 0]);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,potential);
model = struct('muInv',1,'sourceCurrent',[0 0 0],'boundary',boundary);
solution = tdgl.solvers.solveMagnetostatic(mesh,model);

error = solution.cellMagneticFluxDensity-targetB;
verifyLessThan(testCase,max(vecnorm(error,2,2)),1e-10);
verifyLessThan(testCase,solution.diagnostics.freeResidual,1e-10);
verifyLessThan(testCase,solution.diagnostics.gaugeResidual,1e-10);
end

function testNamedBoundarySelection(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
faces = tdgl.boundary.selectFaces(mesh,["xmin","xmax"]);
verifyEqual(testCase,numel(faces),4);
end
