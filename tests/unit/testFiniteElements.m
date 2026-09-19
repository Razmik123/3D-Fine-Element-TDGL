function tests = testFiniteElements
%TESTFINITEELEMENTS Verify P1 and first-family Nedelec element kernels.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testQuadratureAndP1Basis(testCase)
[barycentric,weights] = tdgl.elements.tetraQuadrature(2);
[values,gradients] = tdgl.elements.lagrangeP1(barycentric);
verifyEqual(testCase,sum(weights),1/6,'AbsTol',1e-15);
verifyEqual(testCase,sum(values,2),ones(4,1),'AbsTol',1e-15);
verifyEqual(testCase,sum(gradients,1),zeros(1,3),'AbsTol',1e-15);

integralLambdaSquared = sum(weights.*values(:,1).^2);
verifyEqual(testCase,integralLambdaSquared,1/60,'AbsTol',1e-15);
end

function testNedelecEdgeMoments(testCase)
gradients = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
vertices = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
localEdges = [1 2;1 3;1 4;2 3;2 4;3 4];
moments = zeros(6);

for integrationEdge = 1:6
    i = localEdges(integrationEdge,1);
    j = localEdges(integrationEdge,2);
    barycentric = zeros(1,4);
    barycentric([i j]) = 0.5;
    tangent = vertices(j,:)-vertices(i,:);
    [basis,~] = tdgl.elements.nedelec1(barycentric,gradients);
    vectors = squeeze(basis(1,:,:)).';
    moments(integrationEdge,:) = vectors*tangent.';
end

verifyEqual(testCase,moments,eye(6),'AbsTol',1e-14);
end

function testPhysicalCurlIdentity(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
[~,curlCurl] = tdgl.assembly.edgeMassCurlCurl(mesh);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
verifyLessThan(testCase,norm(curlCurl*G,'fro'),1e-11);
verifyLessThan(testCase,norm(curlCurl-curlCurl.','fro'),1e-13);
end

function testFourthOrderQuadrature(testCase)
[barycentric,weights] = tdgl.elements.tetraQuadrature(4);
verifyEqual(testCase,sum(weights),1/6,'AbsTol',2e-15);
verifyEqual(testCase,sum(weights.*barycentric(:,1).^4),1/210, ...
    'AbsTol',2e-15);
verifyEqual(testCase,sum(weights.*barycentric(:,1).^2.* ...
    barycentric(:,2).^2),1/1260,'AbsTol',2e-15);
end
