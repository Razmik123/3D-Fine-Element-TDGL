function result = vortexFaces(mesh,orderParameter,edgePotential,options)
%VORTEXFACES Detect integer phase winding through oriented mesh faces.
%   The branch decision uses gauge-covariant edge phase differences. The
%   returned integer is invariant under nodal discrete gauge transforms.

arguments
    mesh struct
    orderParameter (:,1) double
    edgePotential (:,1) double
    options.AmplitudeTolerance (1,1) double {mustBeNonnegative} = 1e-10
end
nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
if numel(orderParameter) ~= nNodes || numel(edgePotential) ~= nEdges
    error('tdgl:post:FieldSizeMismatch', ...
        'Order parameter and vector potential do not match the mesh.');
end

edges = double(mesh.topology.edges);
phaseDifference = angle(conj(orderParameter(edges(:,1))).* ...
    orderParameter(edges(:,2)));
covariantDifference = wrapAngle(phaseDifference-edgePotential);
faceEdges = double(mesh.topology.faceEdges);
faceSigns = double(mesh.topology.faceEdgeSigns);
gaugeCirculation = sum(faceSigns.*covariantDifference(faceEdges),2);
magneticFlux = sum(faceSigns.*edgePotential(faceEdges),2);
winding = round((gaugeCirculation+magneticFlux)/(2*pi));

faces = double(mesh.topology.faces);
valid = all(abs(orderParameter(faces)) > options.AmplitudeTolerance,2);
winding(~valid) = 0;
ids = find(winding ~= 0);
result = struct( ...
    'faceIds',ids, ...
    'winding',winding(ids), ...
    'allWinding',winding, ...
    'gaugeCirculation',gaugeCirculation, ...
    'magneticFlux',magneticFlux, ...
    'valid',valid);
end

function angleValue = wrapAngle(angleValue)
angleValue = atan2(sin(angleValue),cos(angleValue));
end
