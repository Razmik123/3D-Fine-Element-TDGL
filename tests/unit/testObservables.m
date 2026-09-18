function tests = testObservables
%TESTOBSERVABLES Verify surface and gauge-invariant voltage measurements.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testSurfaceAverageOfLinearField(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
field = 3 + mesh.nodes(:,2) + 2*mesh.nodes(:,3);
faces = tdgl.boundary.selectFaces(mesh,"xmin");
average = tdgl.observe.surfaceAverageScalar(mesh,field,faces);
verifyEqual(testCase,average,4.5,'AbsTol',1e-13);
end

function testPathVoltageIsGaugeInvariant(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
startNode = 1;
endNode = size(mesh.nodes,1);
path = tdgl.topology.shortestEdgePath(mesh,startNode,endNode);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
rng(11);
Aprevious = randn(double(mesh.topology.nEdges),1);
Acurrent = randn(double(mesh.topology.nEdges),1);
phi = randn(size(mesh.nodes,1),1);
chiPrevious = 0.1*randn(size(mesh.nodes,1),1);
chiCurrent = 0.1*randn(size(mesh.nodes,1),1);
dt = 0.07;

reference = tdgl.observe.pathVoltage( ...
    mesh,phi,Aprevious,Acurrent,dt,path);
transformed = tdgl.observe.pathVoltage(mesh, ...
    phi-(chiCurrent-chiPrevious)/dt, ...
    Aprevious+G*chiPrevious,Acurrent+G*chiCurrent,dt,path);
verifyEqual(testCase,transformed,reference,'AbsTol',5e-13);
end

function testElectrochemicalPotentialIsGaugeInvariant(testCase)
rng(17);
n = 20;
psiPrevious = 0.5+rand(n,1) .* exp(1i*0.2*randn(n,1));
psiCurrent = 0.5+rand(n,1) .* exp(1i*0.2*randn(n,1));
phi = randn(n,1);
chiPrevious = 0.05*randn(n,1);
chiCurrent = 0.05*randn(n,1);
dt = 0.2;

reference = tdgl.observe.electrochemicalPotential( ...
    psiPrevious,psiCurrent,phi,dt);
transformed = tdgl.observe.electrochemicalPotential( ...
    psiPrevious.*exp(1i*chiPrevious), ...
    psiCurrent.*exp(1i*chiCurrent), ...
    phi-(chiCurrent-chiPrevious)/dt,dt);
verifyEqual(testCase,transformed,reference,'AbsTol',2e-14);
end
