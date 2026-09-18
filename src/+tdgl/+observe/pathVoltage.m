function voltage = pathVoltage(mesh, scalarPotential, edgePotentialPrevious, edgePotentialCurrent, dt, signedEdgePath)
%PATHVOLTAGE Gauge-invariant transient voltage along a mesh-edge path.
%   The signed path is oriented from terminal/node A to B. The returned
%   voltage V_A - V_B is
%       phi(A)-phi(B) - d/dt integral_A^B A dot dl.

arguments
    mesh struct
    scalarPotential (:,1) double
    edgePotentialPrevious (:,1) double
    edgePotentialCurrent (:,1) double
    dt (1,1) double {mustBePositive}
    signedEdgePath (:,1) {mustBeInteger}
end

nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
if numel(scalarPotential) ~= nNodes || ...
        numel(edgePotentialPrevious) ~= nEdges || ...
        numel(edgePotentialCurrent) ~= nEdges
    error('tdgl:observe:FieldSizeMismatch', ...
        'Potential arrays do not match the mesh degrees of freedom.');
end
[startNode,endNode] = validatePath(mesh,signedEdgePath);
edgeIds = abs(double(signedEdgePath));
signs = sign(double(signedEdgePath));
linePrevious = sum(signs.*edgePotentialPrevious(edgeIds));
lineCurrent = sum(signs.*edgePotentialCurrent(edgeIds));
voltage = scalarPotential(startNode)-scalarPotential(endNode) - ...
    (lineCurrent-linePrevious)/dt;
end

function [startNode,endNode] = validatePath(mesh,signedPath)
if isempty(signedPath)
    error('tdgl:observe:EmptyPath','Voltage path cannot be empty.');
end
edges = double(mesh.topology.edges);
startNode = NaN;
previousHead = NaN;
for index = 1:numel(signedPath)
    signedId = double(signedPath(index));
    edgeId = abs(signedId);
    if edgeId < 1 || edgeId > size(edges,1)
        error('tdgl:observe:InvalidPathEdge','Path contains an invalid edge ID.');
    end
    oriented = edges(edgeId,:);
    if signedId < 0, oriented = fliplr(oriented); end
    if index == 1
        startNode = oriented(1);
    elseif oriented(1) ~= previousHead
        error('tdgl:observe:DisconnectedPath', ...
            'Signed path edges do not form a continuous oriented path.');
    end
    previousHead = oriented(2);
end
endNode = previousHead;
end
