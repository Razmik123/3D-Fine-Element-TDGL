function tests = testCoupledSolver
%TESTCOUPLEDSOLVER Verify staggered self-consistent time stepping.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testUniformEquilibriumIsStationary(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
state = zeroFieldState(mesh,ones(size(mesh.nodes,1),1));
step = tdgl.solvers.stepCoupledStaggered( ...
    mesh,state,referenceModel(),'TimeStep',0.08, ...
    'Boundary',zeroBoundary(mesh));

verifyEqual(testCase,step.state.orderParameter,state.orderParameter, ...
    'AbsTol',3e-12);
verifyLessThan(testCase,norm(step.state.edgePotential),3e-12);
verifyLessThan(testCase,norm(step.state.scalarPotential),3e-12);
verifyTrue(testCase,step.diagnostics.converged);
verifyLessThan(testCase,step.diagnostics.electromagnetic.currentContinuityResidual, ...
    3e-12);
end

function testHomogeneousRelaxationRemainsSelfConsistent(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
state = zeroFieldState(mesh,0.3*ones(size(mesh.nodes,1),1));
step = tdgl.solvers.stepCoupledStaggered( ...
    mesh,state,referenceModel(),'TimeStep',0.07, ...
    'Boundary',zeroBoundary(mesh));

verifyGreaterThan(testCase,mean(abs(step.state.orderParameter)),0.3);
verifyLessThan(testCase,norm(step.state.edgePotential),3e-12);
verifyLessThan(testCase,norm(step.state.scalarPotential),3e-12);
verifyLessThan(testCase,step.diagnostics.history(end,4),2e-8);
end

function testAppliedFieldIsSolvedInVacuumDomain(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2],[-1 1;-1 1;-1 1]);
targetB = [0.13 -0.04 0.09];
potential = tdgl.physics.uniformFieldPotential(targetB,[0 0 0]);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,potential);
model = referenceModel();
model.conductivity = 0;
model.K = 0;
model.a = -1;
state = zeroFieldState(mesh,zeros(size(mesh.nodes,1),1));

step = tdgl.solvers.stepCoupledStaggered( ...
    mesh,state,model,'TimeStep',0.1,'Boundary',boundary);

error = step.cellMagneticFluxDensity-targetB;
verifyLessThan(testCase,max(vecnorm(error,2,2)),2e-10);
verifyLessThan(testCase,step.diagnostics.electromagnetic.linearResidual,2e-10);
verifyGreaterThanOrEqual(testCase,step.diagnostics.iterations,2);
end

function testCoupledIntegratorStoresFieldsAndObservations(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
state = zeroFieldState(mesh,ones(size(mesh.nodes,1),1));
observer = struct('Name',"meanMagnitude",'Evaluate', ...
    @(~,~,current,~,~,~,~) mean(abs(current.orderParameter)));
result = tdgl.time.integrateCoupled( ...
    mesh,state,referenceModel(),zeroBoundary(mesh), ...
    'StopTime',0.15,'TimeStep',0.05,'StoreEvery',2, ...
    'Observers',observer);

verifyEqual(testCase,result.time,[0 0.1 0.15],'AbsTol',2e-14);
verifySize(testCase,result.orderParameter,[size(mesh.nodes,1),3]);
verifySize(testCase,result.edgePotential,[double(mesh.topology.nEdges),3]);
verifyEqual(testCase,result.observations.meanMagnitude,ones(1,3), ...
    'AbsTol',3e-12);
verifyEqual(testCase,numel(result.stepDiagnostics),3);
end

function testElectricFieldForPureScalarGradient(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
phi = mesh.nodes(:,1)-2*mesh.nodes(:,3);
electricField = tdgl.post.edgeElectricField(mesh, ...
    zeros(double(mesh.topology.nEdges),1), ...
    zeros(double(mesh.topology.nEdges),1),phi,0.2);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
verifyEqual(testCase,electricField,-G*phi,'AbsTol',2e-14);
end

function testTerminalElectrochemicalVoltage(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
previous = ones(size(mesh.nodes,1),1);
phaseRate = 0.2+0.4*mesh.nodes(:,1);
current = exp(1i*0.1*phaseRate).*previous;
phi = 0.7-0.3*mesh.nodes(:,1);
voltage = tdgl.observe.terminalElectrochemicalVoltage( ...
    mesh,tdgl.boundary.selectFaces(mesh,"xmax"), ...
    tdgl.boundary.selectFaces(mesh,"xmin"), ...
    previous,current,phi,0.1);
verifyEqual(testCase,voltage,0.1,'AbsTol',2e-14);
end

function state = zeroFieldState(mesh,psi)
state = tdgl.state.create(mesh,psi, ...
    zeros(double(mesh.topology.nEdges),1),zeros(size(mesh.nodes,1),1));
end

function model = referenceModel
model = struct('u',1,'a',1,'b',1,'K',1, ...
    'conductivity',1,'muInv',1,'kappa',1,'sourceCurrent',[0 0 0]);
end

function boundary = zeroBoundary(mesh)
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,[0 0 0]);
end
