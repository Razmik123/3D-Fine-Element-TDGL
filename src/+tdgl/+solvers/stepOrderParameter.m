function step = stepOrderParameter(mesh,previous,model,options)
%STEPORDERPARAMETER Fully implicit backward-Euler TDGL order-parameter step.
%   This solver advances psi for prescribed A and phi. The temporal link
%   exp(-i*phi*dt) preserves the nodal time-gauge transformation exactly.
%   Coupling A and phi back to Maxwell is implemented in a later module.

arguments
    mesh struct
    previous (:,1) double
    model struct
    options.TimeStep (1,1) double {mustBePositive}
    options.InitialGuess = []
    options.FixedNodeIds = []
    options.FixedValues = []
    options.RelativeTolerance (1,1) double {mustBePositive} = 1e-10
    options.AbsoluteTolerance (1,1) double {mustBePositive} = 1e-12
    options.MaximumIterations (1,1) double {mustBeInteger,mustBePositive} = 25
    options.MinimumDamping (1,1) double {mustBePositive} = 2^-12
end

nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
if numel(previous) ~= nNodes
    error('tdgl:solvers:NodalFieldSizeMismatch', ...
        'previous must contain one complex value per node.');
end
required = {'u','a','b','K','edgePotential','scalarPotential'};
if ~all(isfield(model,required))
    error('tdgl:solvers:IncompleteOrderParameterModel', ...
        'Model must define u, a, b, K, edgePotential, and scalarPotential.');
end
if numel(model.edgePotential) ~= nEdges
    error('tdgl:solvers:EdgeFieldSizeMismatch', ...
        'model.edgePotential must contain one value per edge.');
end

phi = model.scalarPotential;
if isscalar(phi), phi = repmat(phi,nNodes,1); end
phi = phi(:);
if numel(phi) ~= nNodes
    error('tdgl:solvers:ScalarPotentialSizeMismatch', ...
        'scalarPotential must be scalar or contain one value per node.');
end

dt = options.TimeStep;
[massU,~] = tdgl.assembly.p1MassStiffness(mesh,model.u,0);
[massA,~] = tdgl.assembly.p1MassStiffness(mesh,model.a,0);
covariant = tdgl.assembly.covariantP1(mesh,model.edgePotential,model.K);
linearMatrix = massU/dt + covariant - massA;
transportedPrevious = exp(-1i*dt*phi).*previous;
rightHandSide = massU*transportedPrevious/dt;

fixedNodes = unique(double(options.FixedNodeIds(:)));
if isempty(fixedNodes)
    fixedValues = complex(zeros(0,1));
else
    fixedValues = options.FixedValues(:);
    if numel(fixedNodes) ~= numel(fixedValues)
        error('tdgl:solvers:InvalidOrderParameterBoundary', ...
            'FixedNodeIds and FixedValues must have equal length.');
    end
end
freeNodes = setdiff((1:nNodes).',fixedNodes);
realFree = [freeNodes; nNodes+freeNodes];

if isempty(options.InitialGuess)
    psi = transportedPrevious;
else
    psi = options.InitialGuess(:);
    if numel(psi) ~= nNodes
        error('tdgl:solvers:InitialGuessSizeMismatch', ...
            'InitialGuess must contain one value per node.');
    end
end
psi(fixedNodes) = fixedValues;

converged = false;
history = zeros(options.MaximumIterations,2);
for iteration = 1:options.MaximumIterations
    [nonlinear,jacobianNonlinear] = ...
        tdgl.assembly.nonlinearOrderParameter(mesh,psi,model.b);
    complexResidual = linearMatrix*psi+nonlinear-rightHandSide;
    residual = [real(complexResidual);imag(complexResidual)];
    residualNorm = norm(residual(realFree));
    tolerance = options.AbsoluteTolerance + options.RelativeTolerance* ...
        max(1,norm([real(rightHandSide(freeNodes));imag(rightHandSide(freeNodes))]));
    history(iteration,:) = [residualNorm,1];
    if residualNorm <= tolerance
        converged = true;
        break;
    end

    realLinear = [real(linearMatrix),-imag(linearMatrix); ...
                  imag(linearMatrix), real(linearMatrix)];
    jacobian = realLinear+jacobianNonlinear;
    increment = zeros(2*nNodes,1);
    increment(realFree) = -jacobian(realFree,realFree)\residual(realFree);

    damping = 1;
    accepted = false;
    while damping >= options.MinimumDamping
        candidate = psi+damping*(increment(1:nNodes)+1i*increment(nNodes+1:end));
        candidate(fixedNodes) = fixedValues;
        candidateNonlinear = tdgl.assembly.nonlinearOrderParameter( ...
            mesh,candidate,model.b);
        candidateResidual = linearMatrix*candidate+candidateNonlinear-rightHandSide;
        candidateNorm = norm([real(candidateResidual(freeNodes)); ...
                              imag(candidateResidual(freeNodes))]);
        if candidateNorm < residualNorm
            psi = candidate;
            accepted = true;
            history(iteration,2) = damping;
            break;
        end
        damping = damping/2;
    end
    if ~accepted
        error('tdgl:solvers:NewtonLineSearchFailed', ...
            'Order-parameter Newton line search failed at iteration %d.',iteration);
    end
end

if ~converged
    error('tdgl:solvers:NewtonDidNotConverge', ...
        'Order-parameter solve did not converge in %d iterations.', ...
        options.MaximumIterations);
end

step = struct();
step.orderParameter = psi;
step.iterations = iteration;
step.residualNorm = residualNorm;
step.tolerance = tolerance;
step.history = history(1:iteration,:);
step.converged = converged;
end
