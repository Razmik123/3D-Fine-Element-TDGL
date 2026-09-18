function average = surfaceAverageScalar(mesh, nodalField, faceIds)
%SURFACEAVERAGESCALAR Area-average a P1 nodal field over triangle faces.

nodalField = nodalField(:);
faceIds = double(faceIds(:));
if numel(nodalField) ~= size(mesh.nodes,1)
    error('tdgl:observe:NodalFieldSizeMismatch', ...
        'nodalField must contain one value per mesh node.');
end
if isempty(faceIds)
    error('tdgl:observe:EmptySurface', ...
        'At least one face is required for a surface average.');
end

faceNodes = double(mesh.topology.faces(faceIds,:));
faceMean = mean(nodalField(faceNodes),2);
areas = mesh.faceAreas(faceIds);
average = sum(areas.*faceMean)/sum(areas);
end
