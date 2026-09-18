function signedEdgeIds = shortestEdgePath(mesh, startNode, endNode)
%SHORTESTEDGEPATH Return an oriented shortest path through mesh edges.

arguments
    mesh struct
    startNode (1,1) double {mustBeInteger,mustBePositive}
    endNode (1,1) double {mustBeInteger,mustBePositive}
end
if startNode > size(mesh.nodes,1) || endNode > size(mesh.nodes,1)
    error('tdgl:topology:InvalidPathEndpoint', ...
        'Path endpoints must be valid mesh node IDs.');
end

edges = double(mesh.topology.edges);
network = graph(edges(:,1),edges(:,2),[],size(mesh.nodes,1));
nodePath = shortestpath(network,startNode,endNode);
if isempty(nodePath)
    error('tdgl:topology:DisconnectedPath', ...
        'No mesh-edge path connects nodes %d and %d.',startNode,endNode);
end

edgeLookup = sparse(edges(:,1),edges(:,2),1:size(edges,1), ...
    size(mesh.nodes,1),size(mesh.nodes,1));
edgeLookup = edgeLookup + edgeLookup.';
signedEdgeIds = zeros(numel(nodePath)-1,1,'int32');
for index = 1:numel(nodePath)-1
    tail = nodePath(index);
    head = nodePath(index+1);
    edgeId = full(edgeLookup(tail,head));
    direction = 1;
    if tail > head, direction = -1; end
    signedEdgeIds(index) = int32(direction*edgeId);
end
end
