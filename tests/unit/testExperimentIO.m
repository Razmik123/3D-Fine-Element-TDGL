function tests = testExperimentIO
%TESTEXPERIMENTIO Verify append, random-access read, and restart files.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
temporaryRoot = fullfile(repositoryRoot,'tmp','io-tests');
if ~isfolder(temporaryRoot), mkdir(temporaryRoot); end
testCase.TestData.temporaryRoot = temporaryRoot;
end

function testAppendReadAndResume(testCase)
runDirectory = string(tempname(testCase.TestData.temporaryRoot));
cleanup = onCleanup(@() removeRun(runDirectory));
mesh = tdgl.geometry.boxMesh([1 1 1]);
configuration = struct('name',"io-test",'drive',@(time,state) time+ ...
    mean(abs(state.orderParameter)));
writer = tdgl.io.ExperimentWriter(runDirectory,mesh,configuration, ...
    'CheckpointEvery',1);
state0 = tdgl.state.create(mesh,ones(size(mesh.nodes,1),1), ...
    zeros(double(mesh.topology.nEdges),1),zeros(size(mesh.nodes,1),1));
state1 = state0;
state1.orderParameter = 0.8*state0.orderParameter+0.1i;
state1.edgePotential(:) = 0.03;
diagnostics = struct('iterations',3,'history',[0.2 0.2;1e-9 1e-9]);
writer.append(0,0,state0,struct());
writer.append(0.1,1,state1,diagnostics);
writer.writeCheckpoint(0.1,1,state1,diagnostics);
writer.complete();

verifyTrue(testCase,isfile(fullfile(runDirectory,'trajectory.h5')));
verifyTrue(testCase,isfile(fullfile(runDirectory,'checkpoint.mat')));
snapshot = tdgl.io.readSnapshot(runDirectory,2);
verifyEqual(testCase,snapshot.time,0.1,'AbsTol',0);
verifyEqual(testCase,snapshot.state.orderParameter,state1.orderParameter, ...
    'AbsTol',0);
verifyEqual(testCase,snapshot.state.edgePotential,state1.edgePotential, ...
    'AbsTol',0);
checkpoint = tdgl.io.loadCheckpoint(runDirectory);
verifyEqual(testCase,checkpoint.stepIndex,1);
verifyEqual(testCase,checkpoint.state.orderParameter,state1.orderParameter, ...
    'AbsTol',0);

resumed = tdgl.io.ExperimentWriter(runDirectory,mesh,configuration, ...
    'CheckpointEvery',1,'Resume',true);
verifyEqual(testCase,resumed.SnapshotCount,2);
state2 = state1;
state2.scalarPotential(:) = 0.2;
resumed.append(0.2,2,state2,diagnostics);
resumed.complete();
snapshot = tdgl.io.readSnapshot(runDirectory,3);
verifyEqual(testCase,snapshot.state.scalarPotential,state2.scalarPotential, ...
    'AbsTol',0);
clear cleanup;
removeRun(runDirectory);
end

function testMeshFingerprintChangesWithGeometry(testCase)
meshA = tdgl.geometry.boxMesh([1 1 1]);
meshB = meshA;
meshB.nodes(1,1) = meshB.nodes(1,1)+0.01;
verifyNotEqual(testCase,tdgl.io.meshFingerprint(meshA), ...
    tdgl.io.meshFingerprint(meshB));
end

function testStreamingIntegratorSeparatesSnapshotsAndCheckpoints(testCase)
runDirectory = string(tempname(testCase.TestData.temporaryRoot));
cleanup = onCleanup(@() removeRun(runDirectory));
mesh = tdgl.geometry.boxMesh([1 1 1]);
writer = tdgl.io.ExperimentWriter(runDirectory,mesh,struct('name',"stream"), ...
    'CheckpointEvery',1);
state = tdgl.state.create(mesh,ones(size(mesh.nodes,1),1), ...
    zeros(double(mesh.topology.nEdges),1),zeros(size(mesh.nodes,1),1));
model = struct('u',1,'a',1,'b',1,'K',1,'conductivity',1, ...
    'muInv',1,'kappa',1,'sourceCurrent',[0 0 0]);
boundary = tdgl.boundary.tangentialDirichlet( ...
    mesh,mesh.topology.boundaryFaceIds,[0 0 0]);
result = tdgl.time.integrateCoupled(mesh,state,model,boundary, ...
    'StopTime',0.2,'TimeStep',0.05,'StoreEvery',2, ...
    'StoreHistory',false,'Writer',writer);

verifyEmpty(testCase,result.time);
times = h5read(fullfile(runDirectory,'trajectory.h5'),'/time');
verifyEqual(testCase,times,[0 0.1 0.2],'AbsTol',2e-14);
checkpoint = tdgl.io.loadCheckpoint(runDirectory);
verifyEqual(testCase,checkpoint.stepIndex,4);
verifyEqual(testCase,checkpoint.time,0.2,'AbsTol',2e-14);
clear cleanup;
removeRun(runDirectory);
end

function removeRun(runDirectory)
if isfolder(runDirectory), rmdir(runDirectory,'s'); end
end
