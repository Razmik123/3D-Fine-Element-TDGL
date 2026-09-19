function components = conductingComponents(mesh,conductingCellMask)
%CONDUCTINGCOMPONENTS Return nodal components of the conducting submesh.

conductingCellMask = logical(conductingCellMask(:));
if numel(conductingCellMask) ~= size(mesh.cells,1)
    error('tdgl:topology:CellMaskSizeMismatch', ...
        'conductingCellMask must contain one value per tetrahedron.');
end
if ~any(conductingCellMask)
    components = cell(0,1);
    return;
end

cellEdges = double(mesh.topology.cellEdges(conductingCellMask,:));
edgeIds = unique(cellEdges(:));
edges = double(mesh.topology.edges(edgeIds,:));
nodeIds = unique(edges(:));
network = graph(edges(:,1),edges(:,2),[],size(mesh.nodes,1));
labels = conncomp(network,'OutputForm','vector');
activeLabels = unique(labels(nodeIds));
components = cell(numel(activeLabels),1);
for index = 1:numel(activeLabels)
    components{index} = nodeIds(labels(nodeIds) == activeLabels(index));
end
end
