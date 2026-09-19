function tests = testCylinderBenchmark
%TESTCYLINDERBENCHMARK Verify remote uniform-field cylinder infrastructure.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testCylinderHasSampleVacuumAndRemoteBoundary(testCase)
mesh = tdgl.geometry.cylinderInVacuum(0.6,1.0,1.5,2.4,'Hmax',0.7);
verifyTrue(testCase,any(mesh.regionIds == 1));
verifyTrue(testCase,any(mesh.regionIds == 2));
verifyNotEmpty(testCase,mesh.interfaceFaceIds);
outer = tdgl.boundary.selectFaces(mesh, ...
    ["outer-side","outer-top","outer-bottom"]);
verifyEqual(testCase,sort(outer), ...
    sort(double(mesh.topology.boundaryFaceIds)));
end

function testRemoteUniformFieldIsExactlyRepresentedInVacuum(testCase)
mesh = tdgl.geometry.cylinderInVacuum(0.5,0.8,1.3,2.0,'Hmax',0.75);
targetB = [0 0 0.17];
potential = tdgl.physics.uniformFieldPotential(targetB);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,potential);
solution = tdgl.solvers.solveMagnetostatic(mesh,struct( ...
    'muInv',1,'sourceCurrent',[0 0 0],'boundary',boundary), ...
    'Execution',tdgl.compute.execution("cpu"));
verifyLessThan(testCase,max(vecnorm( ...
    solution.cellMagneticFluxDensity-targetB,2,2)),2e-10);
end

function testVortexFaceDetectionIsGaugeInvariant(testCase)
mesh = tdgl.geometry.boxMesh([4 4 2],[-1 1;-1 1;-0.5 0.5]);
x = mesh.nodes(:,1)-0.11;
y = mesh.nodes(:,2)-0.17;
psi = (x+1i*y)./max(hypot(x,y),0.2);
A = zeros(double(mesh.topology.nEdges),1);
base = tdgl.post.vortexFaces(mesh,psi,A,'AmplitudeTolerance',1e-12);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
chi = 0.13*x-0.07*y+0.04*mesh.nodes(:,3);
transformed = tdgl.post.vortexFaces( ...
    mesh,exp(1i*chi).*psi,A+G*chi,'AmplitudeTolerance',1e-12);
verifyEqual(testCase,transformed.allWinding,base.allWinding);
verifyGreaterThan(testCase,nnz(base.allWinding),0);
segments = tdgl.post.vortexSegments(mesh,base);
verifyGreaterThan(testCase,size(segments.startPoint,1),0);
end
