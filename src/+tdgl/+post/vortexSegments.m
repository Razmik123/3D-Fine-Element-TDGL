function segments = vortexSegments(mesh,vortexFaceData)
%VORTEXSEGMENTS Connect pairs of pierced faces inside tetrahedra.
%   Complex cells with more than two pierced faces are reported separately
%   instead of being connected by an arbitrary topology-changing choice.

faceWinding = vortexFaceData.allWinding;
faceCenters = (mesh.nodes(double(mesh.topology.faces(:,1)),:)+ ...
    mesh.nodes(double(mesh.topology.faces(:,2)),:)+ ...
    mesh.nodes(double(mesh.topology.faces(:,3)),:))/3;
nCells = size(mesh.cells,1);
startPoint = zeros(nCells,3);
endPoint = zeros(nCells,3);
strength = zeros(nCells,1);
sourceCell = zeros(nCells,1);
complexCells = zeros(0,1);
count = 0;
for cellId = 1:nCells
    faceIds = double(mesh.topology.cellFaces(cellId,:));
    pierced = faceIds(faceWinding(faceIds) ~= 0);
    if numel(pierced) == 2
        count = count+1;
        startPoint(count,:) = faceCenters(pierced(1),:);
        endPoint(count,:) = faceCenters(pierced(2),:);
        strength(count) = max(abs(faceWinding(pierced)));
        sourceCell(count) = cellId;
    elseif numel(pierced) > 2
        complexCells(end+1,1) = cellId; %#ok<AGROW>
    end
end
segments = struct( ...
    'startPoint',startPoint(1:count,:), ...
    'endPoint',endPoint(1:count,:), ...
    'strength',strength(1:count), ...
    'cellIds',sourceCell(1:count), ...
    'complexCellIds',complexCells);
end
