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
    options.ActiveCellMask = []
    options.Execution = []
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
activeCells = activeCellMask(options.ActiveCellMask,size(mesh.cells,1));
if any(activeCells)
    activeNodes = unique(double(mesh.cells(activeCells,:)));
else
    step = emptyDomainStep(nNodes);
    return;
end
inactiveNodes = setdiff((1:nNodes).',activeNodes);
u = tdgl.assembly.cellCoefficient(model.u,size(mesh.cells,1),'u').*activeCells;
a = tdgl.assembly.cellCoefficient(model.a,size(mesh.cells,1),'a').*activeCells;
b = tdgl.assembly.cellCoefficient(model.b,size(mesh.cells,1),'b').*activeCells;
K = tdgl.assembly.cellCoefficient(model.K,size(mesh.cells,1),'K').*activeCells;
[massU,~] = tdgl.assembly.p1MassStiffness(mesh,u,0);
[massA,~] = tdgl.assembly.p1MassStiffness(mesh,a,0);
covariant = tdgl.assembly.covariantP1(mesh,model.edgePotential,K);
linearMatrix = massU/dt + covariant - massA;
previous(inactiveNodes) = 0;
transportedPrevious = exp(-1i*dt*phi).*previous;
rightHandSide = massU*transportedPrevious/dt;

fixedNodes = double(options.FixedNodeIds(:));
if isempty(fixedNodes)
    fixedValues = complex(zeros(0,1));
else
    fixedValues = options.FixedValues(:);
    if numel(fixedNodes) ~= numel(fixedValues) || ...
            numel(unique(fixedNodes)) ~= numel(fixedNodes) || ...
            any(~ismember(fixedNodes,activeNodes))
        error('tdgl:solvers:InvalidOrderParameterBoundary', ...
            ['Fixed psi nodes must be unique nodes in the active GL domain ', ...
             'and match FixedValues.']);
    end
end
freeNodes = setdiff(activeNodes,fixedNodes);
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
psi(inactiveNodes) = 0;
psi(fixedNodes) = fixedValues;

converged = false;
history = zeros(options.MaximumIterations,2);
linearSolve = struct();
for iteration = 1:options.MaximumIterations
    [nonlinear,jacobianNonlinear] = ...
        tdgl.assembly.nonlinearOrderParameter(mesh,psi,b);
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
    [increment(realFree),linearSolve] = tdgl.compute.solveLinear( ...
        jacobian(realFree,realFree),-residual(realFree),options.Execution);

    damping = 1;
    accepted = false;
    while damping >= options.MinimumDamping
        candidate = psi+damping*(increment(1:nNodes)+1i*increment(nNodes+1:end));
        candidate(fixedNodes) = fixedValues;
        candidate(inactiveNodes) = 0;
        candidateNonlinear = tdgl.assembly.nonlinearOrderParameter( ...
            mesh,candidate,b);
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
step.linearSolve = linearSolve;
end

function mask = activeCellMask(raw,nCells)
if isempty(raw)
    mask = true(nCells,1);
else
    mask = logical(raw(:));
    if numel(mask) ~= nCells
        error('tdgl:solvers:ActiveCellMaskSizeMismatch', ...
            'ActiveCellMask must contain one value per tetrahedron.');
    end
end
end

function step = emptyDomainStep(nNodes)
step = struct();
step.orderParameter = complex(zeros(nNodes,1));
step.iterations = 0;
step.residualNorm = 0;
step.tolerance = 0;
step.history = zeros(0,2);
step.converged = true;
step.linearSolve = struct();
end
