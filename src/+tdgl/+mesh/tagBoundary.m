function mesh = tagBoundary(mesh, tag, selector, options)
%TAGBOUNDARY Assign a physical name to selected exterior faces.
%   SELECTOR receives (centroids,normals,vertexCoordinates) and returns one
%   logical value per exterior face. vertexCoordinates is F-by-3-by-3.

arguments
    mesh struct
    tag (1,1) string
    selector (1,1) function_handle
    options.AllowOverwrite (1,1) logical = false
end
if strlength(tag) == 0
    error('tdgl:mesh:EmptyBoundaryTag','Boundary tag cannot be empty.');
end

faceIds = double(mesh.topology.boundaryFaceIds);
faceNodes = double(mesh.topology.faces(faceIds,:));
vertexCoordinates = zeros(numel(faceIds),3,3);
for vertex = 1:3
    vertexCoordinates(:,:,vertex) = mesh.nodes(faceNodes(:,vertex),:);
end
centroids = mean(vertexCoordinates,3);
normals = mesh.boundaryNormals(faceIds,:);
selected = selector(centroids,normals,vertexCoordinates);
selected = selected(:);
if numel(selected) ~= numel(faceIds) || ~islogical(selected)
    error('tdgl:mesh:InvalidBoundarySelector', ...
        'Boundary selector must return one logical value per boundary face.');
end
selectedFaces = faceIds(selected);
if isempty(selectedFaces)
    error('tdgl:mesh:EmptyBoundarySelection', ...
        'Boundary selector for tag "%s" selected no faces.',tag);
end
if ~options.AllowOverwrite && any(mesh.faceTags(selectedFaces) ~= "")
    error('tdgl:mesh:BoundaryTagConflict', ...
        'Boundary selector overlaps an existing nonempty tag.');
end
mesh.faceTags(selectedFaces) = tag;
end
