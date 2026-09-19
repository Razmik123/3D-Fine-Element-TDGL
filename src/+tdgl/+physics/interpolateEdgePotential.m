function edgePotential = interpolateEdgePotential(mesh,vectorPotential)
%INTERPOLATEEDGEPOTENTIAL Interpolate A into all oriented edge integrals.

nEdges = double(mesh.topology.nEdges);
edges = double(mesh.topology.edges);
startPoint = mesh.nodes(edges(:,1),:);
endPoint = mesh.nodes(edges(:,2),:);
tangent = endPoint-startPoint;
edgePotential = zeros(nEdges,1);
gaussCoordinate = [-1/sqrt(3),1/sqrt(3)];
for q = 1:2
    parameter = 0.5*(1+gaussCoordinate(q));
    points = startPoint+parameter*tangent;
    if isa(vectorPotential,'function_handle')
        values = vectorPotential(points);
    else
        values = vectorPotential;
    end
    if isequal(size(values),[1 3])
        values = repmat(values,nEdges,1);
    end
    if ~isequal(size(values),size(points))
        error('tdgl:physics:InvalidVectorPotential', ...
            'Vector potential must return one three-component row per point.');
    end
    edgePotential = edgePotential+0.5*sum(values.*tangent,2);
end
end
