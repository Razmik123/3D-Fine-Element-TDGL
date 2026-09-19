function state = create(mesh,orderParameter,edgePotential,scalarPotential)
%CREATE Construct and validate one coupled TDGL-MQS state.

arguments
    mesh struct
    orderParameter (:,1) double
    edgePotential (:,1) double
    scalarPotential (:,1) double
end

nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
if numel(orderParameter) ~= nNodes || numel(scalarPotential) ~= nNodes || ...
        numel(edgePotential) ~= nEdges
    error('tdgl:state:FieldSizeMismatch', ...
        'State fields must match the mesh node and edge counts.');
end
if ~isreal(edgePotential) || ~isreal(scalarPotential)
    error('tdgl:state:ComplexElectromagneticPotential', ...
        'Vector and scalar electromagnetic potentials must be real.');
end

state = struct( ...
    'orderParameter',orderParameter(:), ...
    'edgePotential',edgePotential(:), ...
    'scalarPotential',scalarPotential(:));
end
