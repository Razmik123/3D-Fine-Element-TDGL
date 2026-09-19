function tests = testOrderParameterSolver
%TESTORDERPARAMETERSOLVER Verify implicit nonlinear TDGL substeps.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testEquilibriumStateIsStationary(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
model = referenceModel(mesh);
previous = ones(size(mesh.nodes,1),1);
step = tdgl.solvers.stepOrderParameter(mesh,previous,model,'TimeStep',0.1);
verifyEqual(testCase,step.orderParameter,previous,'AbsTol',2e-12);
verifyTrue(testCase,step.converged);
end

function testHomogeneousStateRelaxesAndEnergyFalls(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
model = referenceModel(mesh);
psi = 0.35*ones(size(mesh.nodes,1),1);
initialEnergy = tdgl.post.orderParameterEnergy(mesh,psi,model);

for index = 1:5
    step = tdgl.solvers.stepOrderParameter(mesh,psi,model,'TimeStep',0.08);
    psi = step.orderParameter;
end
finalEnergy = tdgl.post.orderParameterEnergy(mesh,psi,model);

verifyGreaterThan(testCase,mean(abs(psi)),0.35);
verifyLessThan(testCase,mean(abs(psi)),1);
verifyLessThan(testCase,finalEnergy.total,initialEnergy.total);
end

function testUniformScalarPotentialUsesTemporalLink(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
model = referenceModel(mesh);
model.scalarPotential = 0.7;
dt = 0.13;
previous = ones(size(mesh.nodes,1),1);
step = tdgl.solvers.stepOrderParameter(mesh,previous,model,'TimeStep',dt);
expected = exp(-1i*model.scalarPotential*dt)*previous;
verifyEqual(testCase,step.orderParameter,expected,'AbsTol',2e-12);
end

function testCovariantMatrixIsHermitian(testCase)
mesh = tdgl.geometry.boxMesh([2 1 1]);
rng(23);
edgePotential = randn(double(mesh.topology.nEdges),1);
matrix = tdgl.assembly.covariantP1(mesh,edgePotential,1);
verifyLessThan(testCase,norm(matrix-matrix','fro'),2e-12);
verifyGreaterThanOrEqual(testCase,min(eig(full(matrix))),-2e-12);
end

function testTimeIntegratorAndObserver(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
model = referenceModel(mesh);
initial = 0.4*ones(size(mesh.nodes,1),1);
observer = struct('Name',"meanMagnitude",'Evaluate', ...
    @(~,~,current,~,~,~) mean(abs(current)));
result = tdgl.time.integrateOrderParameter(mesh,initial,model, ...
    'StopTime',0.2,'TimeStep',0.05,'StoreEvery',2,'Observers',observer);
verifyEqual(testCase,result.time,[0 0.1 0.2],'AbsTol',1e-14);
verifyEqual(testCase,numel(result.observations.meanMagnitude),4);
verifyGreaterThan(testCase,result.observations.meanMagnitude(end),0.4);
verifyEqual(testCase,numel(result.stepDiagnostics),4);
end

function testNoProximityRestrictsOrderParameterToActiveCells(testCase)
mesh = tdgl.geometry.boxMesh([3 1 1]);
cellCenters = zeros(size(mesh.cells,1),1);
for cellId = 1:size(mesh.cells,1)
    cellCenters(cellId) = mean(mesh.nodes(double(mesh.cells(cellId,:)),1));
end
activeCells = cellCenters < 0.5;
activeNodes = unique(double(mesh.cells(activeCells,:)));
inactiveOnlyNodes = setdiff((1:size(mesh.nodes,1)).',activeNodes);
model = referenceModel(mesh);
previous = ones(size(mesh.nodes,1),1);

step = tdgl.solvers.stepOrderParameter( ...
    mesh,previous,model,'TimeStep',0.05,'ActiveCellMask',activeCells);

verifyEqual(testCase,step.orderParameter(inactiveOnlyNodes), ...
    zeros(numel(inactiveOnlyNodes),1),'AbsTol',2e-14);
verifyGreaterThan(testCase,min(abs(step.orderParameter(activeNodes))),0.9);
end

function testEmptyGLDomainReturnsZeroState(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
model = referenceModel(mesh);
step = tdgl.solvers.stepOrderParameter( ...
    mesh,ones(size(mesh.nodes,1),1),model, ...
    'TimeStep',0.1,'ActiveCellMask',false(size(mesh.cells,1),1));
verifyEqual(testCase,step.orderParameter, ...
    zeros(size(mesh.nodes,1),1),'AbsTol',0);
verifyEqual(testCase,step.iterations,0);
end

function model = referenceModel(mesh)
model = struct();
model.u = 1;
model.a = 1;
model.b = 1;
model.K = 1;
model.edgePotential = zeros(double(mesh.topology.nEdges),1);
model.scalarPotential = 0;
end
