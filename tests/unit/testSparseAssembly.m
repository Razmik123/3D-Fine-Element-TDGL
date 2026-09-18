function tests = testSparseAssembly
%TESTSPARSEASSEMBLY Verify symmetry, conservation, and scalar reproduction.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testP1Matrices(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
[mass,stiffness] = tdgl.assembly.p1MassStiffness(mesh);
onesVector = ones(size(mesh.nodes,1),1);
verifyEqual(testCase,full(sum(mass,'all')),1,'AbsTol',1e-13);
verifyLessThan(testCase,norm(stiffness*onesVector),1e-12);
verifyLessThan(testCase,norm(mass-mass.','fro'),1e-14);
verifyLessThan(testCase,norm(stiffness-stiffness.','fro'),1e-14);
verifyGreaterThan(testCase,min(diag(mass)),0);
end

function testEdgeMatrices(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
[mass,curlCurl] = tdgl.assembly.edgeMassCurlCurl(mesh);
verifyLessThan(testCase,norm(mass-mass.','fro'),1e-13);
verifyLessThan(testCase,norm(curlCurl-curlCurl.','fro'),1e-12);
verifyGreaterThan(testCase,min(eig(full(mass))),0);
end
