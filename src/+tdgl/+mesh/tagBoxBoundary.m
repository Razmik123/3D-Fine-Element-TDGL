function mesh = tagBoxBoundary(mesh, bounds, tolerance)
%TAGBOXBOUNDARY Label the six exterior faces of an axis-aligned box.

arguments
    mesh struct
    bounds (3,2) double
    tolerance (1,1) double {mustBePositive} = 1e-10
end

scale = max(bounds(:,2)-bounds(:,1));
absoluteTolerance = tolerance * max(scale,1);
boundaryFaces = double(mesh.topology.boundaryFaceIds);
labels = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"];

for index = 1:numel(boundaryFaces)
    faceId = boundaryFaces(index);
    vertices = double(mesh.topology.faces(faceId,:));
    coordinates = mesh.nodes(vertices,:);
    matched = false;
    for dimension = 1:3
        lowLabel = labels(2*dimension-1);
        highLabel = labels(2*dimension);
        if all(abs(coordinates(:,dimension)-bounds(dimension,1)) <= absoluteTolerance)
            mesh.faceTags(faceId) = lowLabel;
            matched = true;
        elseif all(abs(coordinates(:,dimension)-bounds(dimension,2)) <= absoluteTolerance)
            mesh.faceTags(faceId) = highLabel;
            matched = true;
        end
    end
    if ~matched
        error('tdgl:mesh:UnclassifiedBoxFace', ...
            'Boundary face %d could not be assigned to a box side.', faceId);
    end
end
end
