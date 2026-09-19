function mesh = cylinderInVacuum(sampleRadius,sampleHeight, ...
        outerRadius,outerHeight,options)
%CYLINDERINVACUUM Mesh a finite cylinder inside a larger vacuum cylinder.
%   Region 1 is the superconducting sample and region 2 is vacuum. The
%   outer boundary is split into outer-side, outer-top, and outer-bottom.

arguments
    sampleRadius (1,1) double {mustBePositive}
    sampleHeight (1,1) double {mustBePositive}
    outerRadius (1,1) double {mustBePositive}
    outerHeight (1,1) double {mustBePositive}
    options.Hmax (1,1) double {mustBePositive} = sampleRadius/3
    options.Hmin (1,1) double {mustBeNonnegative} = 0
end
if outerRadius <= sampleRadius || outerHeight <= sampleHeight
    error('tdgl:geometry:InsufficientVacuumPadding', ...
        'The outer cylinder must exceed the sample radius and height.');
end
if ~license('test','PDE_Toolbox')
    error('tdgl:geometry:PDEToolboxRequired', ...
        'cylinderInVacuum requires Partial Differential Equation Toolbox.');
end

geometry = multicylinder(outerRadius,outerHeight, ...
    'ZOffset',-outerHeight/2);
sampleGeometry = multicylinder(sampleRadius,sampleHeight, ...
    'ZOffset',-sampleHeight/2);
geometry = addCell(geometry,sampleGeometry);
model = createpde;
model.Geometry = geometry;
meshArguments = {'Hmax',options.Hmax,'GeometricOrder','linear'};
if options.Hmin > 0
    meshArguments = [meshArguments,{'Hmin',options.Hmin}];
end
pdeMesh = generateMesh(model,meshArguments{:});

nodes = pdeMesh.Nodes.';
cells = double(pdeMesh.Elements.');
centers = (nodes(cells(:,1),:)+nodes(cells(:,2),:)+ ...
    nodes(cells(:,3),:)+nodes(cells(:,4),:))/4;
radius = hypot(centers(:,1),centers(:,2));
tolerance = 1e-9*max(outerRadius,outerHeight);
insideSample = radius <= sampleRadius+tolerance & ...
    abs(centers(:,3)) <= sampleHeight/2+tolerance;
regionIds = 2*ones(size(cells,1),1);
regionIds(insideSample) = 1;
if ~any(insideSample) || all(insideSample)
    error('tdgl:geometry:InvalidCylinderPartition', ...
        'Generated mesh does not contain both sample and vacuum cells.');
end

mesh = tdgl.mesh.fromArrays(nodes,cells,'RegionIds',regionIds);
boundaryFaces = double(mesh.topology.boundaryFaceIds);
boundaryNodes = double(mesh.topology.faces(boundaryFaces,:));
faceZ = reshape(mesh.nodes(boundaryNodes,3),size(boundaryNodes));
top = all(abs(faceZ-outerHeight/2) <= tolerance,2);
bottom = all(abs(faceZ+outerHeight/2) <= tolerance,2);
mesh.faceTags(boundaryFaces(~top & ~bottom)) = "outer-side";
mesh.faceTags(boundaryFaces(top)) = "outer-top";
mesh.faceTags(boundaryFaces(bottom)) = "outer-bottom";
mesh.metadata.generator = "tdgl.geometry.cylinderInVacuum";
mesh.metadata.geometry = struct( ...
    'sampleRadius',sampleRadius,'sampleHeight',sampleHeight, ...
    'outerRadius',outerRadius,'outerHeight',outerHeight, ...
    'hmax',options.Hmax);
end
