function step = stepCoupledStaggered(mesh,previous,model,options)
%STEPCOUPLEDSTAGGERED Advance TDGL and MQS with block Gauss-Seidel.
%   This auditable reference coupling alternates a nonlinear implicit psi
%   solve and a linear MQS A-phi-p solve until all primary fields stop
%   changing. It is the robust baseline and future preconditioner for a
%   monolithic Newton solve; it is not an operator-splitting single pass.

arguments
    mesh struct
    previous struct
    model struct
    options.TimeStep (1,1) double {mustBePositive}
    options.Boundary struct
    options.FixedPsiNodeIds = []
    options.FixedPsiValues = []
    options.FixedPhiNodeIds = []
    options.FixedPhiValues = []
    options.InitialGuess struct = struct()
    options.RelativeTolerance (1,1) double {mustBePositive} = 1e-8
    options.AbsoluteTolerance (1,1) double {mustBePositive} = 1e-10
    options.MaximumIterations (1,1) double {mustBeInteger,mustBePositive} = 30
    options.Relaxation (1,1) double {mustBeGreaterThan(options.Relaxation,0), ...
        mustBeLessThanOrEqual(options.Relaxation,1)} = 1
    options.OrderStepOptions struct = struct()
end

previous = validateState(mesh,previous);
required = {'u','a','b','K','conductivity','muInv','kappa','sourceCurrent'};
if ~all(isfield(model,required))
    error('tdgl:solvers:IncompleteCoupledModel', ...
        'The coupled model is missing one or more TDGL-MQS coefficients.');
end

iterate = initialGuess(mesh,previous,options.InitialGuess);
history = zeros(options.MaximumIterations,4);
converged = false;
orderStep = struct();
electromagneticStep = struct();

for iteration = 1:options.MaximumIterations
    orderModel = struct( ...
        'u',model.u,'a',model.a,'b',model.b,'K',model.K, ...
        'edgePotential',iterate.edgePotential, ...
        'scalarPotential',iterate.scalarPotential);
    orderOptions = namedOptions(options.OrderStepOptions);
    orderStep = tdgl.solvers.stepOrderParameter( ...
        mesh,previous.orderParameter,orderModel, ...
        'TimeStep',options.TimeStep, ...
        'InitialGuess',iterate.orderParameter, ...
        'FixedNodeIds',options.FixedPsiNodeIds, ...
        'FixedValues',options.FixedPsiValues,orderOptions{:});
    candidatePsi = orderStep.orderParameter;

    electromagneticStep = tdgl.solvers.stepElectromagnetic( ...
        mesh,previous.edgePotential,candidatePsi,model, ...
        'TimeStep',options.TimeStep,'Boundary',options.Boundary, ...
        'FixedPhiNodeIds',options.FixedPhiNodeIds, ...
        'FixedPhiValues',options.FixedPhiValues);
    candidateA = electromagneticStep.edgePotential;
    candidatePhi = electromagneticStep.scalarPotential;

    changes = [relativeChange(candidatePsi,iterate.orderParameter), ...
               relativeChange(candidateA,iterate.edgePotential), ...
               relativeChange(candidatePhi,iterate.scalarPotential)];
    residual = max(changes);
    history(iteration,:) = [changes,residual];
    if residual <= options.AbsoluteTolerance+options.RelativeTolerance
        iterate = tdgl.state.create(mesh,candidatePsi,candidateA,candidatePhi);
        converged = true;
        break;
    end
    iterate = tdgl.state.create(mesh, ...
        blend(candidatePsi,iterate.orderParameter,options.Relaxation), ...
        blend(candidateA,iterate.edgePotential,options.Relaxation), ...
        blend(candidatePhi,iterate.scalarPotential,options.Relaxation));
end

if ~converged
    error('tdgl:solvers:CouplingDidNotConverge', ...
        'Staggered TDGL-MQS coupling did not converge in %d iterations.', ...
        options.MaximumIterations);
end

step = struct();
step.state = iterate;
step.state.gaugeMultiplier = electromagneticStep.gaugeMultiplier;
step.cellMagneticFluxDensity = electromagneticStep.cellMagneticFluxDensity;
step.edgeElectricField = tdgl.post.edgeElectricField( ...
    mesh,previous.edgePotential,iterate.edgePotential, ...
    iterate.scalarPotential,options.TimeStep);
step.weakCurrents = electromagneticStep.weakCurrents;
step.diagnostics = struct( ...
    'converged',converged, ...
    'iterations',iteration, ...
    'history',history(1:iteration,:), ...
    'orderParameter',rmfield(orderStep,'orderParameter'), ...
    'electromagnetic',electromagneticStep.diagnostics);
end

function state = validateState(mesh,state)
required = {'orderParameter','edgePotential','scalarPotential'};
if ~all(isfield(state,required))
    error('tdgl:solvers:IncompleteState', ...
        'State requires orderParameter, edgePotential, and scalarPotential.');
end
state = tdgl.state.create(mesh,state.orderParameter, ...
    state.edgePotential,state.scalarPotential);
end

function guess = initialGuess(mesh,previous,supplied)
guess = previous;
if isempty(fieldnames(supplied)), return; end
required = {'orderParameter','edgePotential','scalarPotential'};
if ~all(isfield(supplied,required))
    error('tdgl:solvers:IncompleteInitialGuess', ...
        'InitialGuess requires all three primary state fields.');
end
guess = tdgl.state.create(mesh,supplied.orderParameter, ...
    supplied.edgePotential,supplied.scalarPotential);
end

function value = blend(candidate,current,relaxation)
value = relaxation*candidate+(1-relaxation)*current;
end

function value = relativeChange(current,previous)
value = norm(current-previous)/max([1,norm(current),norm(previous)]);
end

function pairs = namedOptions(options)
names = fieldnames(options);
pairs = cell(1,2*numel(names));
for index = 1:numel(names)
    pairs{2*index-1} = names{index};
    pairs{2*index} = options.(names{index});
end
end
