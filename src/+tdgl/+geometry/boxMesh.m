function mesh = boxMesh(subdivisions, bounds)
%BOXMESH Create a conforming tetrahedral mesh of an axis-aligned box.
%   MESH = tdgl.geometry.boxMesh([NX NY NZ]) creates a unit-box mesh.
%   MESH = tdgl.geometry.boxMesh(..., [XMIN XMAX; YMIN YMAX; ZMIN ZMAX])
%   uses the supplied bounds. Each hexahedral grid cell is split into six
%   tetrahedra sharing its lower-to-upper main diagonal.
%
%   This generator is intended for verification and small examples. General
%   device geometries should be imported from a tagged tetrahedral mesh.

arguments
    subdivisions (1,3) double {mustBeInteger, mustBePositive}
    bounds (3,2) double = [0 1; 0 1; 0 1]
end

if any(bounds(:,2) <= bounds(:,1))
    error('tdgl:geometry:InvalidBounds', ...
        'Each upper box bound must be greater than its lower bound.');
end

nx = subdivisions(1);
ny = subdivisions(2);
nz = subdivisions(3);
x = linspace(bounds(1,1), bounds(1,2), nx + 1);
y = linspace(bounds(2,1), bounds(2,2), ny + 1);
z = linspace(bounds(3,1), bounds(3,2), nz + 1);
[X, Y, Z] = ndgrid(x, y, z);
nodes = [X(:), Y(:), Z(:)];

nodeId = reshape(uint32(1:size(nodes,1)), nx + 1, ny + 1, nz + 1);
cells = zeros(6 * nx * ny * nz, 4, 'uint32');
cursor = 0;

for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            v000 = nodeId(i,   j,   k);
            v100 = nodeId(i+1, j,   k);
            v010 = nodeId(i,   j+1, k);
            v110 = nodeId(i+1, j+1, k);
            v001 = nodeId(i,   j,   k+1);
            v101 = nodeId(i+1, j,   k+1);
            v011 = nodeId(i,   j+1, k+1);
            v111 = nodeId(i+1, j+1, k+1);

            localCells = [ ...
                v000 v100 v110 v111; ...
                v000 v110 v010 v111; ...
                v000 v010 v011 v111; ...
                v000 v011 v001 v111; ...
                v000 v001 v101 v111; ...
                v000 v101 v100 v111];
            cells(cursor + (1:6), :) = localCells;
            cursor = cursor + 6;
        end
    end
end

mesh = tdgl.mesh.fromArrays(nodes, cells, 'RegionIds', ones(size(cells,1),1));
mesh = tdgl.mesh.tagBoxBoundary(mesh, bounds);
mesh.metadata.generator = "tdgl.geometry.boxMesh";
mesh.metadata.bounds = bounds;
mesh.metadata.subdivisions = subdivisions;
end
