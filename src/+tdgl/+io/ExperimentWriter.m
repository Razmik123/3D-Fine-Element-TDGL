classdef ExperimentWriter < handle
    %EXPERIMENTWRITER Append-only HDF5 trajectory and atomic restart state.

    properties (SetAccess = private)
        RunDirectory (1,1) string
        TrajectoryFile (1,1) string
        CheckpointFile (1,1) string
        ManifestFile (1,1) string
        SnapshotCount (1,1) double = 0
        CheckpointEvery (1,1) double = 1
        MeshFingerprint (1,1) string
        Status (1,1) string = "running"
    end

    properties (Access = private)
        Mesh struct
        Configuration
        Manifest struct
    end

    methods
        function writer = ExperimentWriter(runDirectory,mesh,configuration,options)
            arguments
                runDirectory (1,1) string
                mesh struct
                configuration
                options.CheckpointEvery (1,1) double ...
                    {mustBeInteger,mustBePositive} = 1
                options.Resume (1,1) logical = false
            end
            writer.RunDirectory = runDirectory;
            writer.TrajectoryFile = fullfile(runDirectory,"trajectory.h5");
            writer.CheckpointFile = fullfile(runDirectory,"checkpoint.mat");
            writer.ManifestFile = fullfile(runDirectory,"manifest.json");
            writer.CheckpointEvery = options.CheckpointEvery;
            writer.Mesh = mesh;
            writer.Configuration = configuration;
            writer.MeshFingerprint = tdgl.io.meshFingerprint(mesh);

            if options.Resume
                writer.resumeExisting();
            else
                writer.createNew();
            end
        end

        function append(writer,time,stepIndex,state,diagnostics)
            %APPEND Commit one primary-field snapshot; time is written last.
            arguments
                writer
                time (1,1) double
                stepIndex (1,1) double {mustBeInteger,mustBeNonnegative}
                state struct
                diagnostics struct = struct()
            end
            state = tdgl.state.create(writer.Mesh,state.orderParameter, ...
                state.edgePotential,state.scalarPotential);
            column = writer.SnapshotCount+1;
            h5write(writer.TrajectoryFile,'/psi/real', ...
                real(state.orderParameter),[1 column],[numel(state.orderParameter) 1]);
            h5write(writer.TrajectoryFile,'/psi/imag', ...
                imag(state.orderParameter),[1 column],[numel(state.orderParameter) 1]);
            h5write(writer.TrajectoryFile,'/A',state.edgePotential, ...
                [1 column],[numel(state.edgePotential) 1]);
            h5write(writer.TrajectoryFile,'/phi',state.scalarPotential, ...
                [1 column],[numel(state.scalarPotential) 1]);
            [iterations,residual] = diagnosticValues(diagnostics);
            h5write(writer.TrajectoryFile,'/diagnostics/coupling_iterations', ...
                iterations,[1 column],[1 1]);
            h5write(writer.TrajectoryFile,'/diagnostics/coupling_residual', ...
                residual,[1 column],[1 1]);
            h5write(writer.TrajectoryFile,'/step_index',uint64(stepIndex), ...
                [1 column],[1 1]);
            h5write(writer.TrajectoryFile,'/time',time,[1 column],[1 1]);
            writer.SnapshotCount = column;
        end

        function checkpointIfDue(writer,time,stepIndex,state,diagnostics)
            %CHECKPOINTIFDUE Save restart state independently of snapshots.
            if mod(stepIndex,writer.CheckpointEvery) == 0
                writer.writeCheckpoint(time,stepIndex,state,diagnostics);
            end
        end

        function writeCheckpoint(writer,time,stepIndex,state,diagnostics)
            %WRITECHECKPOINT Atomically replace the latest restart record.
            temporary = writer.CheckpointFile+".tmp";
            snapshotIndex = writer.SnapshotCount;
            meshFingerprint = writer.MeshFingerprint;
            save(temporary,'time','stepIndex','state','diagnostics', ...
                'snapshotIndex','meshFingerprint','-v7.3');
            [success,message] = movefile(temporary,writer.CheckpointFile,'f');
            if ~success
                error('tdgl:io:CannotCommitCheckpoint', ...
                    'Cannot replace checkpoint: %s',message);
            end
        end

        function complete(writer)
            writer.Status = "complete";
            writer.Manifest.status = writer.Status;
            writer.Manifest.snapshotCount = writer.SnapshotCount;
            writer.Manifest.completedUtc = string(datetime('now','TimeZone','UTC', ...
                'Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
            tdgl.io.writeJsonAtomic(writer.ManifestFile,writer.Manifest);
        end

        function markFailed(writer,exception)
            writer.Status = "failed";
            writer.Manifest.status = writer.Status;
            writer.Manifest.snapshotCount = writer.SnapshotCount;
            writer.Manifest.failure = struct('identifier',string(exception.identifier), ...
                'message',string(exception.message));
            tdgl.io.writeJsonAtomic(writer.ManifestFile,writer.Manifest);
        end
    end

    methods (Access = private)
        function createNew(writer)
            if isfolder(writer.RunDirectory) || isfile(writer.RunDirectory)
                error('tdgl:io:RunDirectoryExists', ...
                    'Refusing to overwrite existing run directory %s.', ...
                    writer.RunDirectory);
            end
            [success,message] = mkdir(writer.RunDirectory);
            if ~success
                error('tdgl:io:CannotCreateRun','Cannot create run: %s',message);
            end
            mesh = writer.Mesh;
            configuration = writer.Configuration;
            save(fullfile(writer.RunDirectory,'mesh.mat'),'mesh','-v7.3');
            save(fullfile(writer.RunDirectory,'configuration.mat'), ...
                'configuration','-v7.3');
            writer.createTrajectory();
            writer.Manifest = struct( ...
                'schemaVersion',"1.0", ...
                'status',writer.Status, ...
                'meshFingerprint',writer.MeshFingerprint, ...
                'snapshotCount',0, ...
                'mesh',struct('nodes',size(mesh.nodes,1), ...
                    'edges',double(mesh.topology.nEdges), ...
                    'cells',size(mesh.cells,1)), ...
                'configuration',tdgl.io.jsonSafe(configuration), ...
                'provenance',tdgl.io.provenance());
            tdgl.io.writeJsonAtomic(writer.ManifestFile,writer.Manifest);
        end

        function createTrajectory(writer)
            nNodes = size(writer.Mesh.nodes,1);
            nEdges = double(writer.Mesh.topology.nEdges);
            createUnlimited(writer.TrajectoryFile,'/time',1,'double');
            createUnlimited(writer.TrajectoryFile,'/step_index',1,'uint64');
            createUnlimited(writer.TrajectoryFile,'/psi/real',nNodes,'double');
            createUnlimited(writer.TrajectoryFile,'/psi/imag',nNodes,'double');
            createUnlimited(writer.TrajectoryFile,'/A',nEdges,'double');
            createUnlimited(writer.TrajectoryFile,'/phi',nNodes,'double');
            createUnlimited(writer.TrajectoryFile, ...
                '/diagnostics/coupling_iterations',1,'double');
            createUnlimited(writer.TrajectoryFile, ...
                '/diagnostics/coupling_residual',1,'double');
            h5writeatt(writer.TrajectoryFile,'/','schema_version','1.0');
            h5writeatt(writer.TrajectoryFile,'/','mesh_fingerprint', ...
                char(writer.MeshFingerprint));
        end

        function resumeExisting(writer)
            required = [writer.TrajectoryFile,writer.CheckpointFile, ...
                writer.ManifestFile,fullfile(writer.RunDirectory,"mesh.mat")];
            if ~isfolder(writer.RunDirectory) || ~all(isfile(required))
                error('tdgl:io:IncompleteRun','Run directory is not restartable.');
            end
            storedMesh = load(fullfile(writer.RunDirectory,'mesh.mat'),'mesh');
            if tdgl.io.meshFingerprint(storedMesh.mesh) ~= writer.MeshFingerprint
                error('tdgl:io:MeshFingerprintMismatch', ...
                    'Restart mesh differs from the stored mesh.');
            end
            writer.Manifest = jsondecode(fileread(writer.ManifestFile));
            timeInfo = h5info(writer.TrajectoryFile,'/time');
            writer.SnapshotCount = timeInfo.Dataspace.Size(2);
            writer.Status = "running";
            if isfield(writer.Manifest,'resumeCount')
                resumeCount = double(writer.Manifest.resumeCount)+1;
            else
                resumeCount = 1;
            end
            configuration = writer.Configuration;
            resumeFile = fullfile(writer.RunDirectory,sprintf( ...
                'configuration-resume-%03d.mat',resumeCount));
            save(resumeFile,'configuration','-v7.3');
            writer.Manifest.status = writer.Status;
            writer.Manifest.snapshotCount = writer.SnapshotCount;
            writer.Manifest.resumeCount = resumeCount;
            writer.Manifest.lastResumeConfiguration = ...
                tdgl.io.jsonSafe(writer.Configuration);
            writer.Manifest.resumedUtc = string(datetime('now','TimeZone','UTC', ...
                'Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
            tdgl.io.writeJsonAtomic(writer.ManifestFile,writer.Manifest);
        end
    end
end

function createUnlimited(fileName,path,rowCount,dataType)
chunkRows = min(rowCount,65536);
h5create(fileName,path,[rowCount Inf], ...
    'ChunkSize',[chunkRows 1],'Deflate',1,'Datatype',dataType);
end

function [iterations,residual] = diagnosticValues(diagnostics)
iterations = NaN;
residual = NaN;
if isfield(diagnostics,'iterations'), iterations = diagnostics.iterations; end
if isfield(diagnostics,'history') && ~isempty(diagnostics.history)
    residual = diagnostics.history(end,end);
elseif isfield(diagnostics,'couplingResidual')
    residual = diagnostics.couplingResidual;
end
end
