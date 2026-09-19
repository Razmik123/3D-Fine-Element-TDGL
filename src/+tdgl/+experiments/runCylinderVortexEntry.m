function output = runCylinderVortexEntry(configuration)
%RUNCYLINDERVORTEXENTRY Run or resume the remote-field cylinder experiment.

arguments
    configuration struct
end
required = {'name','geometry','physics','time','solver','execution','storage'};
if ~all(isfield(configuration,required))
    error('tdgl:experiments:IncompleteConfiguration', ...
        'Cylinder experiment configuration is incomplete.');
end

resumeDirectory = string(configuration.storage.ResumeDirectory);
caseOptions = [namedOptions(configuration.geometry), ...
               namedOptions(configuration.physics)];
checkpoint = struct();
if strlength(resumeDirectory) > 0
    checkpoint = tdgl.io.loadCheckpoint(resumeDirectory);
    caseOptions = [caseOptions,{'Mesh',checkpoint.mesh}];
end
setup = tdgl.benchmarks.cylinderVortexEntryCase(caseOptions{:});
configuration.case = setup.configuration;
if strlength(resumeDirectory) == 0
    runDirectory = tdgl.io.makeRunDirectory( ...
        configuration.storage.ResultsRoot,configuration.name);
    writer = tdgl.io.ExperimentWriter(runDirectory,setup.mesh,configuration, ...
        'CheckpointEvery',configuration.storage.CheckpointEvery);
    initialState = setup.initialState;
    startTime = configuration.time.StartTime;
    initialStep = 0;
    writeInitial = true;
else
    runDirectory = resumeDirectory;
    initialState = checkpoint.state;
    startTime = checkpoint.time;
    initialStep = checkpoint.stepIndex;
    writeInitial = false;
    writer = tdgl.io.ExperimentWriter(runDirectory,setup.mesh,configuration, ...
        'CheckpointEvery',configuration.storage.CheckpointEvery,'Resume',true);
end
if startTime >= configuration.time.StopTime
    error('tdgl:experiments:RunAlreadyComplete', ...
        'Restart time is not earlier than the requested stop time.');
end

stepOptions = configuration.solver;
stepOptions.Execution = configuration.execution;
vortexObserver = struct('Name',"vortexFaceCount",'Evaluate', ...
    @(mesh,~,current,~,~,~,~) countVortexFaces(mesh,current));
try
    result = tdgl.time.integrateCoupled( ...
        setup.mesh,initialState,setup.model,setup.boundaryProvider, ...
        'StartTime',startTime, ...
        'StopTime',configuration.time.StopTime, ...
        'TimeStep',configuration.time.TimeStep, ...
        'StoreEvery',configuration.time.StoreEvery, ...
        'StoreHistory',false, ...
        'Writer',writer, ...
        'InitialStepIndex',initialStep, ...
        'WriteInitialSnapshot',writeInitial, ...
        'Observers',vortexObserver, ...
        'StepOptions',stepOptions);
catch exception
    writer.markFailed(exception);
    rethrow(exception);
end

output = struct( ...
    'runDirectory',runDirectory, ...
    'mesh',setup.mesh, ...
    'finalState',result.finalState, ...
    'stepDiagnostics',result.stepDiagnostics, ...
    'observations',result.observations);
end

function count = countVortexFaces(mesh,state)
data = tdgl.post.vortexFaces(mesh,state.orderParameter,state.edgePotential);
count = numel(data.faceIds);
end

function pairs = namedOptions(options)
names = fieldnames(options);
pairs = cell(1,2*numel(names));
for index = 1:numel(names)
    pairs{2*index-1} = names{index};
    pairs{2*index} = options.(names{index});
end
end
