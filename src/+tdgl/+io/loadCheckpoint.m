function checkpoint = loadCheckpoint(runDirectory)
%LOADCHECKPOINT Load and validate the latest atomic restart record.

arguments
    runDirectory (1,1) string
end
meshRecord = load(fullfile(runDirectory,"mesh.mat"),'mesh');
checkpoint = load(fullfile(runDirectory,"checkpoint.mat"));
expected = tdgl.io.meshFingerprint(meshRecord.mesh);
if string(checkpoint.meshFingerprint) ~= expected
    error('tdgl:io:CheckpointMeshMismatch', ...
        'Checkpoint and stored mesh fingerprints differ.');
end
checkpoint.mesh = meshRecord.mesh;
end
