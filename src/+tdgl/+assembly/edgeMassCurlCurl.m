function [massMatrix, curlCurlMatrix] = edgeMassCurlCurl(mesh, massCoefficient, curlCoefficient)
%EDGEMASSCURLCURL Assemble lowest-order Nedelec mass and curl-curl matrices.

arguments
    mesh struct
    massCoefficient = 1
    curlCoefficient = 1
end

nCells = size(mesh.cells,1);
nEdges = double(mesh.topology.nEdges);
massCoefficient = tdgl.assembly.cellCoefficient( ...
    massCoefficient, nCells, 'massCoefficient');
curlCoefficient = tdgl.assembly.cellCoefficient( ...
    curlCoefficient, nCells, 'curlCoefficient');
[barycentric, referenceWeights] = tdgl.elements.tetraQuadrature(2);

nEntries = 36*nCells;
rows = zeros(nEntries,1);
columns = zeros(nEntries,1);
massValues = zeros(nEntries,1);
curlValues = zeros(nEntries,1);
cursor = 0;

for cellId = 1:nCells
    edgeIds = double(mesh.topology.cellEdges(cellId,:));
    signs = double(mesh.topology.cellEdgeSigns(cellId,:));
    gradients = mesh.geometry.gradLambda(:,:,cellId);
    [basis, curls] = tdgl.elements.nedelec1(barycentric, gradients);
    basis = basis .* reshape(signs,1,1,6);
    curls = curls .* signs.';

    determinant = mesh.geometry.detJ(cellId);
    localMass = zeros(6);
    for q = 1:numel(referenceWeights)
        vectors = squeeze(basis(q,:,:)).';
        localMass = localMass + referenceWeights(q)*determinant * ...
            (vectors*vectors.');
    end
    localMass = massCoefficient(cellId) * localMass;
    localCurl = curlCoefficient(cellId) * mesh.geometry.volume(cellId) * ...
        (curls*curls.');

    indices = cursor + (1:36);
    [localRows, localColumns] = ndgrid(edgeIds,edgeIds);
    rows(indices) = localRows(:);
    columns(indices) = localColumns(:);
    massValues(indices) = localMass(:);
    curlValues(indices) = localCurl(:);
    cursor = cursor + 36;
end

massMatrix = sparse(rows, columns, massValues, nEdges, nEdges);
curlCurlMatrix = sparse(rows, columns, curlValues, nEdges, nEdges);
end
