function tests = testElectromagneticSolver
%TESTELECTROMAGNETICSOLVER Verify transient MQS and current continuity.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testSupercurrentAssemblyForUniformState(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
psi = exp(0.37i)*ones(size(mesh.nodes,1),1);
[phaseLoad,densityMass] = tdgl.assembly.supercurrent(mesh,psi,1);
[edgeMass,~] = tdgl.assembly.edgeMassCurlCurl(mesh,1,0);

verifyLessThan(testCase,norm(phaseLoad),2e-13);
verifyLessThan(testCase,norm(densityMass-edgeMass,'fro'),2e-13);
verifyLessThan(testCase,norm(densityMass-densityMass.','fro'),2e-13);
end

function testVacuumUniformFieldReproduction(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2],[-1 1;-1 1;-1 1]);
targetB = [0.19 -0.07 0.11];
potential = tdgl.physics.uniformFieldPotential(targetB,[0 0 0]);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,potential);
model = referenceModel(mesh);
model.conductivity = 0;
model.K = 0;

step = tdgl.solvers.stepElectromagnetic( ...
    mesh,zeros(double(mesh.topology.nEdges),1), ...
    zeros(size(mesh.nodes,1),1),model, ...
    'TimeStep',0.1,'Boundary',boundary);

error = step.cellMagneticFluxDensity-targetB;
verifyLessThan(testCase,max(vecnorm(error,2,2)),2e-10);
verifyLessThan(testCase,step.diagnostics.linearResidual,2e-10);
verifyLessThan(testCase,step.diagnostics.gaugeResidual,2e-10);
verifyEmpty(testCase,step.scalarPotential(step.scalarPotential ~= 0));
end

function testZeroConductingEquilibrium(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
boundary = zeroBoundary(mesh);
model = referenceModel(mesh);
previousA = zeros(double(mesh.topology.nEdges),1);
psi = ones(size(mesh.nodes,1),1);

step = tdgl.solvers.stepElectromagnetic( ...
    mesh,previousA,psi,model,'TimeStep',0.05,'Boundary',boundary);

verifyLessThan(testCase,norm(step.edgePotential),2e-12);
verifyLessThan(testCase,norm(step.scalarPotential),2e-12);
verifyLessThan(testCase,step.diagnostics.linearResidual,2e-12);
verifyLessThan(testCase,step.diagnostics.currentContinuityResidual,2e-12);
end

function testPhaseGradientRespectsWeakCurrentContinuity(testCase)
mesh = tdgl.geometry.boxMesh([3 2 2]);
boundary = zeroBoundary(mesh);
model = referenceModel(mesh);
psi = exp(0.31i*mesh.nodes(:,1));

step = tdgl.solvers.stepElectromagnetic( ...
    mesh,zeros(double(mesh.topology.nEdges),1),psi,model, ...
    'TimeStep',0.04,'Boundary',boundary);

verifyLessThan(testCase,step.diagnostics.linearResidual,2e-10);
verifyLessThan(testCase,step.diagnostics.currentContinuityResidual,2e-10);
verifyLessThan(testCase,step.diagnostics.gaugeResidual,2e-10);
end

function testRejectsPhiConstraintOutsideConductor(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
model = referenceModel(mesh);
model.conductivity = 0;
call = @() tdgl.solvers.stepElectromagnetic( ...
    mesh,zeros(double(mesh.topology.nEdges),1), ...
    zeros(size(mesh.nodes,1),1),model, ...
    'TimeStep',0.1,'Boundary',zeroBoundary(mesh), ...
    'FixedPhiNodeIds',1,'FixedPhiValues',0);
verifyError(testCase,call,'tdgl:solvers:PhiWithoutConductor');
end

function model = referenceModel(mesh)
model = struct();
model.conductivity = 1;
model.muInv = 1;
model.kappa = 1;
model.K = 1;
model.sourceCurrent = [0 0 0];
assert(double(mesh.topology.nEdges) > 0);
end

function boundary = zeroBoundary(mesh)
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,[0 0 0]);
end
