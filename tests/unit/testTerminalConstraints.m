function tests = testTerminalConstraints
%TESTTERMINALCONSTRAINTS Verify equipotential integral-current terminals.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testCurrentTerminalsMeetIntegralConstraints(testCase)
mesh = tdgl.geometry.boxMesh([3 2 2]);
current = 0.04;
terminals = currentTerminalPair(mesh,current);
model = conductionModel();
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,[0 0 0]);

step = tdgl.solvers.stepElectromagnetic( ...
    mesh,zeros(double(mesh.topology.nEdges),1), ...
    zeros(size(mesh.nodes,1),1),model, ...
    'TimeStep',0.1,'Boundary',boundary,'Terminals',terminals);

verifyEqual(testCase,step.terminal.current,[current;-current], ...
    'AbsTol',2e-11);
verifyTrue(testCase,step.terminal.reference(1));
verifyEqual(testCase,step.scalarPotential(terminals(1).nodeIds), ...
    zeros(numel(terminals(1).nodeIds),1),'AbsTol',2e-12);
verifyLessThan(testCase,range(step.scalarPotential(terminals(2).nodeIds)),2e-12);
verifyLessThan(testCase,step.diagnostics.currentContinuityResidual,2e-10);
end

function testVoltageTerminalsAreStrongEquipotentials(testCase)
mesh = tdgl.geometry.boxMesh([3 2 2]);
terminals = [withNodes(tdgl.problem.terminal( ...
    "ground","xmin","ground",0),mesh,"xmin"); ...
    withNodes(tdgl.problem.terminal( ...
    "drive","xmax","voltage",0.3),mesh,"xmax")];
model = conductionModel();
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,[0 0 0]);

step = tdgl.solvers.stepElectromagnetic( ...
    mesh,zeros(double(mesh.topology.nEdges),1), ...
    zeros(size(mesh.nodes,1),1),model, ...
    'TimeStep',0.1,'Boundary',boundary,'Terminals',terminals);

verifyEqual(testCase,step.scalarPotential(terminals(1).nodeIds), ...
    zeros(numel(terminals(1).nodeIds),1),'AbsTol',2e-13);
verifyEqual(testCase,step.scalarPotential(terminals(2).nodeIds), ...
    0.3*ones(numel(terminals(2).nodeIds),1),'AbsTol',2e-13);
verifyLessThan(testCase,abs(sum(step.terminal.current)),2e-11);
end

function testUnbalancedIsolatedCurrentTerminalsAreRejected(testCase)
mesh = tdgl.geometry.boxMesh([2 1 1]);
terminals = currentTerminalPair(mesh,0.1);
terminals(2).excitation = -0.07;
mask = true(size(mesh.cells,1),1);
call = @() tdgl.boundary.scalarPotentialSpace(mesh,mask,terminals);
verifyError(testCase,call,'tdgl:boundary:UnbalancedComponentCurrent');
end

function terminals = currentTerminalPair(mesh,current)
terminals = [withNodes(tdgl.problem.terminal( ...
    "source","xmin","current",current),mesh,"xmin"); ...
    withNodes(tdgl.problem.terminal( ...
    "drain","xmax","current",-current),mesh,"xmax")];
end

function terminal = withNodes(terminal,mesh,tag)
faceIds = tdgl.boundary.selectFaces(mesh,tag);
terminal.faceIds = faceIds;
terminal.nodeIds = unique(double(mesh.topology.faces(faceIds,:)));
end

function model = conductionModel
model = struct('conductivity',1,'muInv',1,'kappa',1, ...
    'K',0,'sourceCurrent',[0 0 0]);
end
