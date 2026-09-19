function analysis = analyzeCylinderRun(runDirectory,snapshotIndex)
%ANALYZECYLINDERRUN Reconstruct fields and vortices from a saved snapshot.

arguments
    runDirectory (1,1) string
    snapshotIndex (1,1) double {mustBeInteger,mustBePositive}
end
meshRecord = load(fullfile(runDirectory,'mesh.mat'),'mesh');
mesh = meshRecord.mesh;
snapshot = tdgl.io.readSnapshot(runDirectory,snapshotIndex);
magneticFluxDensity = tdgl.post.cellMagneticFluxDensity( ...
    mesh,snapshot.state.edgePotential);
vortexFaces = tdgl.post.vortexFaces( ...
    mesh,snapshot.state.orderParameter,snapshot.state.edgePotential);
segments = tdgl.post.vortexSegments(mesh,vortexFaces);
cellCenters = (mesh.nodes(double(mesh.cells(:,1)),:)+ ...
    mesh.nodes(double(mesh.cells(:,2)),:)+ ...
    mesh.nodes(double(mesh.cells(:,3)),:)+ ...
    mesh.nodes(double(mesh.cells(:,4)),:))/4;
sampleCells = double(mesh.regionIds) == 1;

analysis = struct();
analysis.runDirectory = runDirectory;
analysis.snapshotIndex = snapshotIndex;
analysis.time = snapshot.time;
analysis.stepIndex = snapshot.stepIndex;
analysis.mesh = mesh;
analysis.state = snapshot.state;
analysis.cellCenters = cellCenters;
analysis.sampleCellMask = sampleCells;
analysis.magneticFluxDensity = magneticFluxDensity;
analysis.magneticMagnitude = vecnorm(magneticFluxDensity,2,2);
analysis.vortexFaces = vortexFaces;
analysis.vortexSegments = segments;
analysis.metrics = struct( ...
    'meanSampleBz',mean(magneticFluxDensity(sampleCells,3)), ...
    'maximumSampleB',max(vecnorm(magneticFluxDensity(sampleCells,:),2,2)), ...
    'minimumSampleOrderMagnitude',minimumSampleOrder(mesh,sampleCells, ...
        snapshot.state.orderParameter), ...
    'piercedFaceCount',numel(vortexFaces.faceIds), ...
    'segmentCount',size(segments.startPoint,1));
end

function value = minimumSampleOrder(mesh,sampleCells,orderParameter)
nodes = unique(double(mesh.cells(sampleCells,:)));
value = min(abs(orderParameter(nodes)));
end
