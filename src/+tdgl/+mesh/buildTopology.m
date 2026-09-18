function topology = buildTopology(cells, nNodes)
%BUILDTOPOLOGY Build globally oriented tetrahedral edges and faces.
%   Global edges point from the smaller to the larger node ID. Global faces
%   use ascending node order. CellEdgeSigns and CellFaceSigns map local
%   oriented basis functions/chains to these global orientations.

cells = uint32(cells);
nCells = size(cells,1);

localEdges = uint32([1 2; 1 3; 1 4; 2 3; 2 4; 3 4]);
allEdges = zeros(6*nCells,2,'uint32');
edgeDirection = zeros(6*nCells,1,'int8');
for localId = 1:6
    rows = (localId-1)*nCells + (1:nCells);
    pair = cells(:,localEdges(localId,:));
    allEdges(rows,:) = sort(pair,2);
    edgeDirection(rows) = int8(2*(pair(:,1) < pair(:,2))-1);
end
[edges, ~, edgeIds] = unique(allEdges, 'rows', 'sorted');
cellEdges = uint32(reshape(edgeIds, nCells, 6));
cellEdgeSigns = reshape(edgeDirection, nCells, 6);

localFaces = uint32([2 3 4; 1 3 4; 1 2 4; 1 2 3]);
boundaryCoefficient = int8([1 -1 1 -1]);
allFaces = zeros(4*nCells,3,'uint32');
orientedFaceSigns = zeros(4*nCells,1,'int8');
for localId = 1:4
    rows = (localId-1)*nCells + (1:nCells);
    triple = cells(:,localFaces(localId,:));
    allFaces(rows,:) = sort(triple,2);
    parity = permutationParity3(triple);
    orientedFaceSigns(rows) = boundaryCoefficient(localId) .* parity;
end
[faces, ~, faceIds] = unique(allFaces, 'rows', 'sorted');
cellFaces = uint32(reshape(faceIds, nCells, 4));
cellFaceSigns = reshape(orientedFaceSigns, nCells, 4);

nFaces = size(faces,1);
faceCells = zeros(nFaces,2,'uint32');
faceLocalIds = zeros(nFaces,2,'uint8');
faceCellCounts = zeros(nFaces,1,'uint8');
for localId = 1:4
    for cellId = 1:nCells
        faceId = double(cellFaces(cellId,localId));
        slot = double(faceCellCounts(faceId)) + 1;
        if slot > 2
            error('tdgl:mesh:NonManifoldFace', ...
                'Face %d belongs to more than two tetrahedra.', faceId);
        end
        faceCells(faceId,slot) = uint32(cellId);
        faceLocalIds(faceId,slot) = uint8(localId);
        faceCellCounts(faceId) = uint8(slot);
    end
end

boundaryFaceIds = uint32(find(faceCellCounts == 1));
interiorFaceIds = uint32(find(faceCellCounts == 2));

faceEdgePairs = [faces(:,[1 2]); faces(:,[1 3]); faces(:,[2 3])];
[present, ids] = ismember(faceEdgePairs, edges, 'rows');
if ~all(present)
    error('tdgl:mesh:InternalTopologyError', ...
        'A face edge was not found in the global edge table.');
end
faceEdges = uint32(reshape(ids, nFaces, 3));
faceEdgeSigns = repmat(int8([1 -1 1]), nFaces, 1);

topology = struct();
topology.nNodes = uint32(nNodes);
topology.nEdges = uint32(size(edges,1));
topology.nFaces = uint32(nFaces);
topology.nCells = uint32(nCells);
topology.edges = edges;
topology.faces = faces;
topology.cellEdges = cellEdges;
topology.cellEdgeSigns = cellEdgeSigns;
topology.cellFaces = cellFaces;
topology.cellFaceSigns = cellFaceSigns;
topology.faceEdges = faceEdges;
topology.faceEdgeSigns = faceEdgeSigns;
topology.faceCells = faceCells;
topology.faceLocalIds = faceLocalIds;
topology.boundaryFaceIds = boundaryFaceIds;
topology.interiorFaceIds = interiorFaceIds;
topology.localEdges = localEdges;
topology.localFaces = localFaces;
end

function parity = permutationParity3(values)
inversions = (values(:,1) > values(:,2)) + ...
    (values(:,1) > values(:,3)) + (values(:,2) > values(:,3));
parity = int8(1 - 2*mod(inversions,2));
end
