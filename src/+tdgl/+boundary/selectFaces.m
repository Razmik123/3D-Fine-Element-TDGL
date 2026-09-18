function faceIds = selectFaces(mesh, tags)
%SELECTFACES Resolve one or more named boundary tags to global face IDs.

tags = string(tags(:));
if isempty(tags) || any(strlength(tags) == 0)
    error('tdgl:boundary:InvalidTagSelection', ...
        'At least one nonempty boundary tag is required.');
end

isBoundary = false(double(mesh.topology.nFaces),1);
isBoundary(double(mesh.topology.boundaryFaceIds)) = true;
selected = ismember(mesh.faceTags, tags) & isBoundary;
faceIds = find(selected);

missing = tags(~ismember(tags, unique(mesh.faceTags(selected))));
if ~isempty(missing)
    error('tdgl:boundary:UnknownTag', ...
        'Unknown or empty boundary tag selection: %s.', strjoin(missing, ', '));
end
end
