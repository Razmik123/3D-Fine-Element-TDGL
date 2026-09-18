function magneticFluxDensity = cellMagneticFluxDensity(mesh, edgePotential)
%CELLMAGNETICFLUXDENSITY Evaluate curl(A) in each tetrahedron.
%   Lowest-order Nedelec fields have a cellwise constant curl.

nEdges = double(mesh.topology.nEdges);
edgePotential = edgePotential(:);
if numel(edgePotential) ~= nEdges
    error('tdgl:post:EdgeFieldSizeMismatch', ...
        'edgePotential must contain one degree of freedom per edge.');
end

nCells = size(mesh.cells,1);
magneticFluxDensity = zeros(nCells,3);
centroid = [1 1 1 1]/4;
for cellId = 1:nCells
    gradients = mesh.geometry.gradLambda(:,:,cellId);
    [~, curls] = tdgl.elements.nedelec1(centroid, gradients);
    edgeIds = double(mesh.topology.cellEdges(cellId,:));
    signs = double(mesh.topology.cellEdgeSigns(cellId,:));
    localDegrees = signs(:) .* edgePotential(edgeIds);
    magneticFluxDensity(cellId,:) = localDegrees.' * curls;
end
end
