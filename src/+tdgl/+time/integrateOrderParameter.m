function result = integrateOrderParameter(mesh,initialState,modelProvider,options)
%INTEGRATEORDERPARAMETER Integrate prescribed-field TDGL in time.
%   MODELPROVIDER may be a static model struct or a function handle
%   MODEL = MODELPROVIDER(TIME, PSI). Observers are struct records with Name
%   and Evaluate fields. Evaluate is called as
%       value = Evaluate(mesh,psiPrevious,psiCurrent,model,time,dt).

arguments
    mesh struct
    initialState (:,1) double
    modelProvider
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
if numel(initialState) ~= size(mesh.nodes,1)
    error('tdgl:time:InitialStateSizeMismatch', ...
        'initialState must contain one value per mesh node.');
end

observers = options.Observers;
observerNames = strings(numel(observers),1);
for index = 1:numel(observers)
    if ~isfield(observers(index),'Name') || ...
            ~isfield(observers(index),'Evaluate') || ...
            ~isa(observers(index).Evaluate,'function_handle')
        error('tdgl:time:InvalidObserver', ...
            'Every observer requires Name and function-handle Evaluate fields.');
    end
    observerNames(index) = string(observers(index).Name);
end
if numel(unique(observerNames)) ~= numel(observerNames)
    error('tdgl:time:DuplicateObserverName', ...
        'Observer names must be unique.');
end

estimatedSteps = ceil((options.StopTime-options.StartTime)/options.TimeStep);
storedStates = complex(zeros(size(mesh.nodes,1),floor(estimatedSteps/options.StoreEvery)+2));
storedTimes = zeros(1,size(storedStates,2));
observationValues = cell(numel(observers),estimatedSteps);
stepDiagnostics = repmat(struct('iterations',0,'residualNorm',0,'dt',0),estimatedSteps,1);

psi = initialState;
time = options.StartTime;
stepId = 0;
storeId = 1;
storedStates(:,storeId) = psi;
storedTimes(storeId) = time;

while time < options.StopTime
    dt = min(options.TimeStep,options.StopTime-time);
    nextTime = time+dt;
    if isa(modelProvider,'function_handle')
        model = modelProvider(nextTime,psi);
    else
        model = modelProvider;
    end
    previous = psi;
    nameValues = namedOptions(options.StepOptions);
    step = tdgl.solvers.stepOrderParameter( ...
        mesh,previous,model,'TimeStep',dt,nameValues{:});
    psi = step.orderParameter;
    stepId = stepId+1;
    stepDiagnostics(stepId) = struct( ...
        'iterations',step.iterations,'residualNorm',step.residualNorm,'dt',dt);

    for observerId = 1:numel(observers)
        observationValues{observerId,stepId} = observers(observerId).Evaluate( ...
            mesh,previous,psi,model,nextTime,dt);
    end
    time = nextTime;
    if mod(stepId,options.StoreEvery) == 0 || time >= options.StopTime
        storeId = storeId+1;
        storedStates(:,storeId) = psi;
        storedTimes(storeId) = time;
    end
end

result = struct();
result.time = storedTimes(1:storeId);
result.orderParameter = storedStates(:,1:storeId);
result.finalState = psi;
result.stepDiagnostics = stepDiagnostics(1:stepId);
result.observations = struct();
for observerId = 1:numel(observers)
    values = observationValues(observerId,1:stepId);
    if all(cellfun(@(item) isnumeric(item) && isscalar(item),values))
        result.observations.(matlab.lang.makeValidName(observerNames(observerId))) = ...
            cell2mat(values);
    else
        result.observations.(matlab.lang.makeValidName(observerNames(observerId))) = values;
    end
end
end

function pairs = namedOptions(options)
names = fieldnames(options);
pairs = cell(1,2*numel(names));
for index = 1:numel(names)
    pairs{2*index-1} = names{index};
    pairs{2*index} = options.(names{index});
end
end
