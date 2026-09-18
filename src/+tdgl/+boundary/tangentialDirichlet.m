function condition = tangentialDirichlet(mesh, faceIds, vectorPotential)
%TANGENTIALDIRICHLET Interpolate tangential A data into edge DOFs.
%   The edge degree of freedom is the oriented line integral of A along the
%   global edge. Two-point Gauss quadrature is exact for linear potentials.

faceIds = unique(double(faceIds(:)));
if isempty(faceIds) || any(faceIds < 1) || ...
        any(faceIds > double(mesh.topology.nFaces))
    error('tdgl:boundary:InvalidFaceIds', ...
        'faceIds must select at least one valid mesh face.');
end

edgeIds = unique(double(mesh.topology.faceEdges(faceIds,:)));
edges = double(mesh.topology.edges(edgeIds,:));
startPoint = mesh.nodes(edges(:,1),:);
endPoint = mesh.nodes(edges(:,2),:);
tangent = endPoint-startPoint;

gaussCoordinate = [-1/sqrt(3), 1/sqrt(3)];
values = zeros(numel(edgeIds),1);
for q = 1:2
    parameter = 0.5*(1+gaussCoordinate(q));
    points = startPoint + parameter*tangent;
    if isa(vectorPotential, 'function_handle')
        potentialValues = vectorPotential(points);
    else
        potentialValues = vectorPotential;
    end
    if isequal(size(potentialValues), [1 3])
        potentialValues = repmat(potentialValues,numel(edgeIds),1);
    end
    if ~isequal(size(potentialValues), size(points))
        error('tdgl:boundary:InvalidVectorPotential', ...
            'Vector potential must return one three-component row per point.');
    end
    values = values + 0.5*sum(potentialValues.*tangent,2);
end

condition = struct('edgeIds',edgeIds,'values',values,'faceIds',faceIds);
end
