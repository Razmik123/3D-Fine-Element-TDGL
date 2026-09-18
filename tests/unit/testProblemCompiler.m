function tests = testProblemCompiler
%TESTPROBLEMCOMPILER Verify explicit materials, interfaces, and terminals.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testBalancedCurrentTerminalsCompile(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
material = tdgl.materials.create("superconductor",'Name',"reference-S");
assignments = struct('regionId',1,'material',material);
terminals = [ ...
    tdgl.problem.terminal("source","xmin","current",2.5); ...
    tdgl.problem.terminal("drain","xmax","current",-2.5)];
experiment = baseExperiment();
experiment.terminals = terminals;

problem = tdgl.problem.compile(mesh,assignments,experiment);
verifyEqual(testCase,numel(problem.terminals),2);
verifyEqual(testCase,[problem.terminals.area],[1 1],'AbsTol',1e-14);
verifyEqual(testCase,sort(problem.outerFaceIds), ...
    sort(double(mesh.topology.boundaryFaceIds)));
end

function testUnbalancedCurrentIsRejected(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
assignments = struct('regionId',1, ...
    'material',tdgl.materials.create("superconductor"));
experiment = baseExperiment();
experiment.terminals = [ ...
    tdgl.problem.terminal("source","xmin","current",2); ...
    tdgl.problem.terminal("drain","xmax","current",-1)];
verifyError(testCase,@() tdgl.problem.compile(mesh,assignments,experiment), ...
    'tdgl:problem:UnbalancedTerminalCurrent');
end

function testGLInterfaceRequiresExplicitLaw(testCase)
base = tdgl.geometry.boxMesh([2 1 1]);
centers = zeros(size(base.cells,1),1);
for cellId = 1:size(base.cells,1)
    centers(cellId) = mean(base.nodes(double(base.cells(cellId,:)),1));
end
regionIds = ones(size(centers));
regionIds(centers >= 0.5) = 2;
mesh = tdgl.mesh.fromArrays(base.nodes,base.cells,'RegionIds',regionIds);
mesh = tdgl.mesh.tagBoxBoundary(mesh,[0 1;0 1;0 1]);
assignments(1) = struct('regionId',1, ...
    'material',tdgl.materials.create("superconductor"));
assignments(2) = struct('regionId',2, ...
    'material',tdgl.materials.create("normal-metal"));
experiment = baseExperiment();

verifyError(testCase,@() tdgl.problem.compile(mesh,assignments,experiment), ...
    'tdgl:problem:MissingInterfaceModel');

experiment.interfaceModels = tdgl.problem.interfaceModel([1 2],"de-gennes", ...
    'GammaB',0.4);
problem = tdgl.problem.compile(mesh,assignments,experiment);
verifyEqual(testCase,problem.interfaces.orderParameterType,"de-gennes");
verifyEqual(testCase,problem.interfaces.gammaB,0.4);
end

function experiment = baseExperiment
experiment = struct();
experiment.name = "unit-test";
experiment.outerBoundaryTags = ...
    ["xmin","xmax","ymin","ymax","zmin","zmax"];
experiment.boundaryConditions = struct([]);
experiment.terminals = struct([]);
experiment.interfaceModels = struct([]);
experiment.time = struct('start',0,'stop',1,'initialStep',0.1);
experiment.initialCondition = struct();
end
