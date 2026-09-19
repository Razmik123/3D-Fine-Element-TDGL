function result = integrateCoupled(mesh,initialState,modelProvider, ...
        boundaryProvider,options)
%INTEGRATECOUPLED Integrate the staggered self-consistent TDGL-MQS system.
%   Providers may be static structs or function handles. A model provider
%   is called as MODEL(TIME,STATE); a boundary provider as
%   BOUNDARY(TIME,STATE,MODEL). Observer Evaluate functions receive
%   (mesh,previous,current,model,time,dt,step).

arguments
    mesh struct
    initialState struct
    modelProvider
    boundaryProvider
    options.StartTime (1,1) double = 0
    options.StopTime (1,1) double {mustBePositive}
    options.TimeStep (1,1) double {mustBePositive}
    options.StoreEvery (1,1) double {mustBeInteger,mustBePositive} = 1
    options.Observers struct = struct([])
    options.StepOptions struct = struct()
end
if options.StopTime <= options.StartTime
    error('tdgl:time:InvalidTimeInterval', ...
        'StopTime must be greater than StartTime.');
end
state = tdgl.state.create(mesh,initialState.orderParameter, ...
    initialState.edgePotential,initialState.scalarPotential);
[observers,observerNames] = validateObservers(options.Observers);

nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
estimatedSteps = ceil((options.StopTime-options.StartTime)/options.TimeStep);
nStores = floor(estimatedSteps/options.StoreEvery)+2;
storedPsi = complex(zeros(nNodes,nStores));
storedA = zeros(nEdges,nStores);
storedPhi = zeros(nNodes,nStores);
storedTimes = zeros(1,nStores);
observations = cell(numel(observers),estimatedSteps);
diagnostics = repmat(struct('iterations',0,'couplingResidual',0,'dt',0), ...
    estimatedSteps,1);

time = options.StartTime;
stepId = 0;
storeId = 1;
[storedPsi(:,1),storedA(:,1),storedPhi(:,1)] = stateFields(state);
storedTimes(1) = time;

while time < options.StopTime
    dt = min(options.TimeStep,options.StopTime-time);
    nextTime = time+dt;
    model = evaluateModel(modelProvider,nextTime,state);
    boundary = evaluateBoundary(boundaryProvider,nextTime,state,model);
    previous = state;
    stepOptions = namedOptions(options.StepOptions);
    step = tdgl.solvers.stepCoupledStaggered( ...
        mesh,previous,model,'TimeStep',dt,'Boundary',boundary, ...
        stepOptions{:});
    state = step.state;
    stepId = stepId+1;
    diagnostics(stepId) = struct( ...
        'iterations',step.diagnostics.iterations, ...
        'couplingResidual',step.diagnostics.history(end,4), ...
        'dt',dt);

    for observerId = 1:numel(observers)
        observations{observerId,stepId} = observers(observerId).Evaluate( ...
            mesh,previous,state,model,nextTime,dt,step);
    end
    time = nextTime;
    if mod(stepId,options.StoreEvery) == 0 || time >= options.StopTime
        storeId = storeId+1;
        [storedPsi(:,storeId),storedA(:,storeId),storedPhi(:,storeId)] = ...
            stateFields(state);
        storedTimes(storeId) = time;
    end
end

result = struct();
result.time = storedTimes(1:storeId);
result.orderParameter = storedPsi(:,1:storeId);
result.edgePotential = storedA(:,1:storeId);
result.scalarPotential = storedPhi(:,1:storeId);
result.finalState = state;
result.stepDiagnostics = diagnostics(1:stepId);
result.observations = packObservations(observerNames, ...
    observations(:,1:stepId));
end

function model = evaluateModel(provider,time,state)
if isa(provider,'function_handle')
    model = provider(time,state);
else
    model = provider;
end
end

function boundary = evaluateBoundary(provider,time,state,model)
if isa(provider,'function_handle')
    boundary = provider(time,state,model);
else
    boundary = provider;
end
end

function [observers,names] = validateObservers(observers)
names = strings(numel(observers),1);
for index = 1:numel(observers)
    if ~isfield(observers(index),'Name') || ...
            ~isfield(observers(index),'Evaluate') || ...
            ~isa(observers(index).Evaluate,'function_handle')
        error('tdgl:time:InvalidObserver', ...
            'Every observer requires Name and function-handle Evaluate fields.');
    end
    names(index) = string(observers(index).Name);
end
if numel(unique(names)) ~= numel(names)
    error('tdgl:time:DuplicateObserverName', ...
        'Observer names must be unique.');
end
end

function packed = packObservations(names,values)
packed = struct();
for index = 1:numel(names)
    row = values(index,:);
    field = matlab.lang.makeValidName(names(index));
    if all(cellfun(@(item) isnumeric(item) && isscalar(item),row))
        packed.(field) = cell2mat(row);
    else
        packed.(field) = row;
    end
end
end

function [psi,edgePotential,scalarPotential] = stateFields(state)
psi = state.orderParameter;
edgePotential = state.edgePotential;
scalarPotential = state.scalarPotential;
end

function pairs = namedOptions(options)
names = fieldnames(options);
pairs = cell(1,2*numel(names));
for index = 1:numel(names)
    pairs{2*index-1} = names{index};
    pairs{2*index} = options.(names{index});
end
end
