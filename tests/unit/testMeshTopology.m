function tests = testMeshTopology
%TESTMESHTOPOLOGY Verify orientations and exact-sequence incidence maps.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testUnitBoxCountsAndTags(testCase)
mesh = tdgl.geometry.boxMesh([1 1 1]);
verifyEqual(testCase,size(mesh.nodes,1),8);
verifyEqual(testCase,size(mesh.cells,1),6);
verifyEqual(testCase,double(mesh.topology.nEdges),19);
verifyEqual(testCase,double(mesh.topology.nFaces),18);
verifyEqual(testCase,numel(mesh.topology.boundaryFaceIds),12);
verifyEqual(testCase,sort(unique(mesh.faceTags(mesh.faceTags ~= ""))), ...
    sort(["xmin";"xmax";"ymin";"ymax";"zmin";"zmax"]));
verifyGreaterThan(testCase,min(mesh.geometry.volume),0);
verifyEqual(testCase,sum(mesh.geometry.volume),1,'AbsTol',1e-14);
end

function testIncidenceComplexIsExact(testCase)
mesh = tdgl.geometry.boxMesh([2 2 2]);
[G,C,D] = tdgl.topology.incidenceMatrices(mesh);
verifyEqual(testCase,nnz(C*G),0);
verifyEqual(testCase,nnz(D*C),0);
verifyEqual(testCase,tdgl.topology.bettiNumbers(mesh),[1 0 0 0]);
end

function testInvertedCellIsReoriented(testCase)
nodes = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
mesh = tdgl.mesh.fromArrays(nodes,[1 3 2 4]);
verifyGreaterThan(testCase,mesh.geometry.detJ,0);
verifyEqual(testCase,mesh.geometry.volume,1/6,'AbsTol',1e-14);
end

function testDegenerateCellIsRejected(testCase)
nodes = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
verifyError(testCase,@() tdgl.mesh.fromArrays(nodes,[1 2 3 4]), ...
    'tdgl:mesh:DegenerateTetrahedron');
end

function testPredicateBoundaryTagger(testCase)
mesh = tdgl.mesh.fromArrays( ...
    [0 0 0;1 0 0;0 1 0;0 0 1],[1 2 3 4]);
mesh = tdgl.mesh.tagBoundary(mesh,"bottom", ...
    @(centroids,normals,vertices) abs(centroids(:,3)) < 1e-14);
verifyEqual(testCase,sum(mesh.faceTags == "bottom"),1);
end
