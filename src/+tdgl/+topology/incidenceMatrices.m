function [G, C, D] = incidenceMatrices(mesh)
%INCIDENCEMATRICES Return discrete gradient, curl, and divergence matrices.
%   G maps vertex cochains to edge cochains, C maps edges to faces, and D
%   maps faces to cells. With the canonical orientations, C*G and D*C are
%   exactly zero in integer arithmetic.

topology = mesh.topology;
nNodes = double(topology.nNodes);
nEdges = double(topology.nEdges);
nFaces = double(topology.nFaces);
nCells = double(topology.nCells);

edgeIds = (1:nEdges).';
edges = double(topology.edges);
G = sparse([edgeIds; edgeIds], [edges(:,1); edges(:,2)], ...
    [-ones(nEdges,1); ones(nEdges,1)], nEdges, nNodes);

faceIds = repmat((1:nFaces).',3,1);
faceEdges = double(topology.faceEdges(:));
faceSigns = double(topology.faceEdgeSigns(:));
C = sparse(faceIds, faceEdges, faceSigns, nFaces, nEdges);

cellIds = repmat((1:nCells).',4,1);
cellFaces = double(topology.cellFaces(:));
cellSigns = double(topology.cellFaceSigns(:));
D = sparse(cellIds, cellFaces, cellSigns, nCells, nFaces);
end
