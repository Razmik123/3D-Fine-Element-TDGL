function snapshot = readSnapshot(runDirectory,index)
%READSNAPSHOT Load one stored state without reading the complete trajectory.

arguments
    runDirectory (1,1) string
    index (1,1) double {mustBeInteger,mustBePositive}
end
fileName = fullfile(runDirectory,"trajectory.h5");
meshRecord = load(fullfile(runDirectory,"mesh.mat"),'mesh');
nNodes = size(meshRecord.mesh.nodes,1);
nEdges = double(meshRecord.mesh.topology.nEdges);
information = h5info(fileName,'/time');
nSnapshots = information.Dataspace.Size(2);
if index > nSnapshots
    error('tdgl:io:SnapshotOutOfRange', ...
        'Snapshot %d exceeds the stored count %d.',index,nSnapshots);
end
realPsi = h5read(fileName,'/psi/real',[1 index],[nNodes 1]);
imagPsi = h5read(fileName,'/psi/imag',[1 index],[nNodes 1]);
snapshot = struct();
snapshot.time = h5read(fileName,'/time',[1 index],[1 1]);
snapshot.stepIndex = double(h5read(fileName,'/step_index',[1 index],[1 1]));
snapshot.state = tdgl.state.create(meshRecord.mesh,complex(realPsi,imagPsi), ...
    h5read(fileName,'/A',[1 index],[nEdges 1]), ...
    h5read(fileName,'/phi',[1 index],[nNodes 1]));
snapshot.diagnostics = struct( ...
    'couplingIterations',h5read(fileName, ...
        '/diagnostics/coupling_iterations',[1 index],[1 1]), ...
    'couplingResidual',h5read(fileName, ...
        '/diagnostics/coupling_residual',[1 index],[1 1]));
end
