function setup = cylinderVortexEntryCase(options)
%CYLINDERVORTEXENTRYCASE Configure a finite cylinder in remote uniform field.
%   The homogeneous field is prescribed only on the surrounding vacuum's
%   far boundary. Screening, demagnetization, and vortices are solved in
%   the complete sample-plus-vacuum domain.

arguments
    options.SampleRadius (1,1) double {mustBePositive} = 2
    options.SampleHeight (1,1) double {mustBePositive} = 4
    options.OuterRadius (1,1) double {mustBePositive} = 6
    options.OuterHeight (1,1) double {mustBePositive} = 10
    options.Hmax (1,1) double {mustBePositive} = 0.4
    options.Kappa (1,1) double {mustBePositive} = 5
    options.MaximumField (1,1) double {mustBeNonnegative} = 0.8
    options.RampTime (1,1) double {mustBePositive} = 10
    options.NoiseAmplitude (1,1) double {mustBeNonnegative} = 1e-3
    options.RandomSeed (1,1) double {mustBeInteger,mustBeNonnegative} = 1729
    options.Mesh = []
end

if isempty(options.Mesh)
    mesh = tdgl.geometry.cylinderInVacuum( ...
        options.SampleRadius,options.SampleHeight, ...
        options.OuterRadius,options.OuterHeight,'Hmax',options.Hmax);
else
    mesh = options.Mesh;
    tdgl.mesh.validate(mesh);
end
sampleCells = double(mesh.regionIds) == 1;
nCells = size(mesh.cells,1);
activeNodes = unique(double(mesh.cells(sampleCells,:)));

previousRandomState = rng;
cleanup = onCleanup(@() rng(previousRandomState));
rng(options.RandomSeed,'twister');
psi = complex(zeros(size(mesh.nodes,1),1));
psi(activeNodes) = 1+options.NoiseAmplitude*randn(numel(activeNodes),1);
clear cleanup;

zeroPotential = tdgl.physics.uniformFieldPotential([0 0 0]);
initialA = tdgl.physics.interpolateEdgePotential(mesh,zeroPotential);
initialState = tdgl.state.create(mesh,psi,initialA, ...
    zeros(size(mesh.nodes,1),1));

model = struct();
model.u = double(sampleCells);
model.a = double(sampleCells);
model.b = double(sampleCells);
model.K = double(sampleCells);
model.glActive = sampleCells;
model.conductivity = double(sampleCells);
model.muInv = ones(nCells,1);
model.kappa = options.Kappa;
model.sourceCurrent = [0 0 0];

outerFaces = double(mesh.topology.boundaryFaceIds);
boundaryProvider = @(time,~,~) remoteUniformBoundary( ...
    mesh,outerFaces,time,options.MaximumField,options.RampTime);
parameters = rmfield(options,'Mesh');
configuration = struct( ...
    'name',"cylinder-vortex-entry", ...
    'physics',"magnetoquasistatic TDGL", ...
    'appliedField',"homogeneous axial field on remote vacuum boundary", ...
    'parameters',parameters);
setup = struct( ...
    'mesh',mesh, ...
    'sampleCellMask',sampleCells, ...
    'model',model, ...
    'initialState',initialState, ...
    'boundaryProvider',boundaryProvider, ...
    'configuration',configuration);
end

function boundary = remoteUniformBoundary(mesh,faceIds,time,maximumField,rampTime)
amplitude = maximumField*min(max(time/rampTime,0),1);
potential = tdgl.physics.uniformFieldPotential([0 0 amplitude]);
boundary = tdgl.boundary.tangentialDirichlet(mesh,faceIds,potential);
end
