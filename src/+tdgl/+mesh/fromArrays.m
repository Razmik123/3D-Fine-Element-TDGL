function mesh = fromArrays(nodes, cells, options)
%FROMARRAYS Compile node/tetrahedron arrays into the canonical mesh model.
%   MESH = tdgl.mesh.fromArrays(NODES, CELLS) accepts N-by-3 coordinates
%   and M-by-4 one-based tetrahedron connectivity. Inverted tetrahedra are
%   reoriented consistently; degenerate cells are rejected.
%
%   Name-value options:
%     RegionIds       M-by-1 positive material-region identifiers.
%     VolumeTolerance Relative degeneracy tolerance (default 1e-12).

arguments
    nodes double
    cells {mustBeNumeric}
    options.RegionIds = []
    options.VolumeTolerance (1,1) double {mustBePositive} = 1e-12
end

if size(nodes,2) ~= 3 && size(nodes,1) == 3
    nodes = nodes.';
end
if size(cells,2) ~= 4 && size(cells,1) == 4
    cells = cells.';
end
if size(nodes,2) ~= 3 || size(cells,2) ~= 4
    error('tdgl:mesh:InvalidArrayShape', ...
        'Nodes must be N-by-3 and cells must be M-by-4.');
end
if isempty(nodes) || isempty(cells) || any(~isfinite(nodes), 'all')
    error('tdgl:mesh:InvalidArrayData', ...
        'Mesh arrays must be nonempty and node coordinates must be finite.');
end
if any(cells(:) < 1) || any(cells(:) ~= round(cells(:))) || ...
        any(cells(:) > size(nodes,1))
    error('tdgl:mesh:InvalidConnectivity', ...
        'Cell connectivity must contain valid one-based integer node IDs.');
end
if any(diff(sort(double(cells), 2), 1, 2) == 0, 'all')
    error('tdgl:mesh:RepeatedCellVertex', ...
        'A tetrahedron cannot contain a repeated vertex.');
end

cells = uint32(cells);
nCells = size(cells,1);
if isempty(options.RegionIds)
    regionIds = ones(nCells,1,'int32');
else
    regionIds = options.RegionIds(:);
    if numel(regionIds) ~= nCells || any(regionIds < 1) || ...
            any(regionIds ~= round(regionIds))
        error('tdgl:mesh:InvalidRegionIds', ...
            'RegionIds must contain one positive integer per tetrahedron.');
    end
    regionIds = int32(regionIds);
end

x1 = nodes(double(cells(:,1)),:);
x2 = nodes(double(cells(:,2)),:);
x3 = nodes(double(cells(:,3)),:);
x4 = nodes(double(cells(:,4)),:);
signedSixVolume = dot(cross(x2-x1, x3-x1, 2), x4-x1, 2);

coordinateScale = max(max(nodes,[],1) - min(nodes,[],1));
coordinateScale = max(coordinateScale, 1);
absoluteTolerance = options.VolumeTolerance * coordinateScale^3;
if any(abs(signedSixVolume) <= absoluteTolerance)
    bad = find(abs(signedSixVolume) <= absoluteTolerance, 1);
    error('tdgl:mesh:DegenerateTetrahedron', ...
        'Tetrahedron %d is degenerate or below the volume tolerance.', bad);
end

inverted = signedSixVolume < 0;
cells(inverted,[3 4]) = cells(inverted,[4 3]);

mesh = struct();
mesh.nodes = nodes;
mesh.cells = cells;
mesh.regionIds = regionIds;
mesh.topology = tdgl.mesh.buildTopology(cells, size(nodes,1));
mesh.geometry = tdgl.elements.tetraGeometry(nodes, cells);
mesh.faceRegionIds = zeros(mesh.topology.nFaces, 2, 'int32');
mesh.faceRegionIds(:,1) = regionIds(double(mesh.topology.faceCells(:,1)));
secondCell = mesh.topology.faceCells(:,2) > 0;
mesh.faceRegionIds(secondCell,2) = ...
    regionIds(double(mesh.topology.faceCells(secondCell,2)));
mesh.interfaceFaceIds = find(secondCell & ...
    mesh.faceRegionIds(:,1) ~= mesh.faceRegionIds(:,2));
mesh.faceTags = strings(mesh.topology.nFaces,1);
mesh.metadata = struct('generator', "tdgl.mesh.fromArrays");

mesh = addFaceGeometry(mesh);
tdgl.mesh.validate(mesh);
end

function mesh = addFaceGeometry(mesh)
faces = double(mesh.topology.faces);
p1 = mesh.nodes(faces(:,1),:);
p2 = mesh.nodes(faces(:,2),:);
p3 = mesh.nodes(faces(:,3),:);
rawNormal = cross(p2-p1, p3-p1, 2);
twiceArea = vecnorm(rawNormal, 2, 2);
mesh.faceAreas = 0.5 * twiceArea;
mesh.faceNormals = rawNormal ./ twiceArea;

mesh.boundaryNormals = zeros(size(mesh.faceNormals));
boundary = double(mesh.topology.boundaryFaceIds);
for index = 1:numel(boundary)
    faceId = boundary(index);
    cellId = double(mesh.topology.faceCells(faceId,1));
    faceCenter = mean(mesh.nodes(faces(faceId,:),:),1);
    cellCenter = mean(mesh.nodes(double(mesh.cells(cellId,:)),:),1);
    normal = mesh.faceNormals(faceId,:);
    if dot(normal, faceCenter-cellCenter) < 0
        normal = -normal;
    end
    mesh.boundaryNormals(faceId,:) = normal;
end
end
