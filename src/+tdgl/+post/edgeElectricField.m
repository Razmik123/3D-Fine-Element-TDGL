function electricField = edgeElectricField(mesh,previousEdgePotential, ...
        currentEdgePotential,scalarPotential,dt)
%EDGEELECTRICFIELD Return oriented edge integrals of E=-dA/dt-grad(phi).

arguments
    mesh struct
    previousEdgePotential (:,1) double
    currentEdgePotential (:,1) double
    scalarPotential (:,1) double
    dt (1,1) double {mustBePositive}
end
nEdges = double(mesh.topology.nEdges);
nNodes = size(mesh.nodes,1);
if numel(previousEdgePotential) ~= nEdges || ...
        numel(currentEdgePotential) ~= nEdges || ...
        numel(scalarPotential) ~= nNodes
    error('tdgl:post:FieldSizeMismatch', ...
        'Potential arrays do not match the mesh degrees of freedom.');
end
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
electricField = -(currentEdgePotential-previousEdgePotential)/dt - ...
    G*scalarPotential;
end
