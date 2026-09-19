function step = stepElectromagnetic(mesh,previousEdgePotential,orderParameter,model,options)
%STEPELECTROMAGNETIC Solve one linearized MQS Maxwell/current step.
%   For fixed psi, this routine solves the coupled edge/nodal system for A,
%   physical scalar potential phi, and a mixed Coulomb-gauge multiplier.
%
%   The weak electromagnetic equation is
%     sigma*(A-A_prev)/dt + sigma*grad(phi)
%     + kappa^2*curl(muInv*curl(A)) - j_s = j_ext.

arguments
    mesh struct
    previousEdgePotential (:,1) double
    orderParameter (:,1) double
    model struct
    options.TimeStep (1,1) double {mustBePositive}
    options.Boundary struct
    options.FixedPhiNodeIds = []
    options.FixedPhiValues = []
end

nEdges = double(mesh.topology.nEdges);
nNodes = size(mesh.nodes,1);
nCells = size(mesh.cells,1);
if numel(previousEdgePotential) ~= nEdges || numel(orderParameter) ~= nNodes
    error('tdgl:solvers:FieldSizeMismatch', ...
        'Previous A and psi sizes must match mesh edge and node counts.');
end
required = {'conductivity','muInv','kappa','K','sourceCurrent'};
if ~all(isfield(model,required))
    error('tdgl:solvers:IncompleteElectromagneticModel', ...
        'Model must define conductivity, muInv, kappa, K, and sourceCurrent.');
end

conductivity = tdgl.assembly.cellCoefficient( ...
    model.conductivity,nCells,'conductivity');
[conductivityMass,~] = tdgl.assembly.edgeMassCurlCurl(mesh,conductivity,0);
[gaugeMass,curlCurl] = tdgl.assembly.edgeMassCurlCurl(mesh,1,model.muInv);
curlCurl = model.kappa^2*curlCurl;
[phaseLoad,densityMass] = tdgl.assembly.supercurrent(mesh,orderParameter,model.K);
externalLoad = tdgl.assembly.edgeLoad(mesh,model.sourceCurrent);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
dt = options.TimeStep;

[fixedEdges,fixedEdgeValues] = normalizedBoundary(options.Boundary,nEdges);
allBoundaryEdges = unique(double(mesh.topology.faceEdges( ...
    double(mesh.topology.boundaryFaceIds),:)));
if ~isequal(fixedEdges,allBoundaryEdges)
    error('tdgl:solvers:IncompleteOuterBoundary', ...
        'The MQS reference step requires A data on every exterior edge.');
end
freeEdges = setdiff((1:nEdges).',fixedEdges);

conductingCells = conductivity > 0;
components = tdgl.topology.conductingComponents(mesh,conductingCells);
if isempty(components)
    conductingNodes = zeros(0,1);
else
    conductingNodes = unique(vertcat(components{:}));
end
[fixedPhiNodes,fixedPhiValues] = normalizedScalarBoundary( ...
    options.FixedPhiNodeIds,options.FixedPhiValues,components,nNodes);
freePhiNodes = setdiff(conductingNodes,fixedPhiNodes);

boundaryNodes = unique(double(mesh.topology.faces( ...
    double(mesh.topology.boundaryFaceIds),:)));
gaugeNodes = setdiff((1:nNodes).',boundaryNodes);

Haa = conductivityMass/dt+curlCurl+densityMass;
baseRightHandSide = conductivityMass*previousEdgePotential/dt + ...
    phaseLoad+externalLoad;
GfreePhi = G(:,freePhiNodes);
GfixedPhi = G(:,fixedPhiNodes);
Ggauge = G(:,gaugeNodes);
gaugeCoupling = gaugeMass*Ggauge;

rhsA = baseRightHandSide(freeEdges) - ...
    Haa(freeEdges,fixedEdges)*fixedEdgeValues;
if ~isempty(fixedPhiNodes)
    rhsA = rhsA-conductivityMass(freeEdges,:)*GfixedPhi*fixedPhiValues;
end

Aphi = conductivityMass(freeEdges,:)*GfreePhi;
Ap = gaugeCoupling(freeEdges,:);

currentOperator = conductivityMass/dt+densityMass;
phiA = GfreePhi.'*currentOperator(:,freeEdges);
phiPhi = GfreePhi.'*conductivityMass*GfreePhi;
rhsPhi = GfreePhi.'*baseRightHandSide - ...
    GfreePhi.'*currentOperator(:,fixedEdges)*fixedEdgeValues;
if ~isempty(fixedPhiNodes)
    rhsPhi = rhsPhi-GfreePhi.'*conductivityMass*GfixedPhi*fixedPhiValues;
end

pA = Ggauge.'*gaugeMass(:,freeEdges);
rhsP = -Ggauge.'*gaugeMass(:,fixedEdges)*fixedEdgeValues;

nPhi = numel(freePhiNodes);
nP = numel(gaugeNodes);
systemMatrix = [ ...
    Haa(freeEdges,freeEdges), Aphi, Ap; ...
    phiA, phiPhi, sparse(nPhi,nP); ...
    pA, sparse(nP,nPhi), sparse(nP,nP)];
rightHandSide = [rhsA;rhsPhi;rhsP];
unknown = systemMatrix\rightHandSide;

edgePotential = zeros(nEdges,1);
edgePotential(fixedEdges) = fixedEdgeValues;
edgePotential(freeEdges) = unknown(1:numel(freeEdges));
offset = numel(freeEdges);
scalarPotential = zeros(nNodes,1);
scalarPotential(fixedPhiNodes) = fixedPhiValues;
scalarPotential(freePhiNodes) = unknown(offset+(1:nPhi));
offset = offset+nPhi;
gaugeMultiplier = zeros(nNodes,1);
gaugeMultiplier(gaugeNodes) = unknown(offset+(1:nP));

normalCurrent = -conductivityMass * ( ...
    (edgePotential-previousEdgePotential)/dt+G*scalarPotential);
supercurrent = phaseLoad-densityMass*edgePotential;
totalCurrent = supercurrent+normalCurrent+externalLoad;
currentContinuity = G(:,freePhiNodes).'*totalCurrent;
gaugeResidual = Ggauge.'*gaugeMass*edgePotential;
linearResidual = systemMatrix*unknown-rightHandSide;

step = struct();
step.edgePotential = edgePotential;
step.scalarPotential = scalarPotential;
step.gaugeMultiplier = gaugeMultiplier;
step.cellMagneticFluxDensity = ...
    tdgl.post.cellMagneticFluxDensity(mesh,edgePotential);
step.weakCurrents = struct( ...
    'superconducting',supercurrent, ...
    'normal',normalCurrent, ...
    'external',externalLoad, ...
    'total',totalCurrent);
step.diagnostics = struct( ...
    'linearResidual',norm(linearResidual)/max(1,norm(rightHandSide)), ...
    'currentContinuityResidual',norm(currentContinuity)/max(1,norm(totalCurrent)), ...
    'gaugeResidual',norm(gaugeResidual)/max(1,norm(edgePotential)), ...
    'nFreeEdges',numel(freeEdges), ...
    'nScalarPotentialDofs',nPhi, ...
    'nGaugeDofs',nP);
step.phiReaction = G(:,fixedPhiNodes).'*totalCurrent;
step.fixedPhiNodeIds = fixedPhiNodes;
end

function [edgeIds,values] = normalizedBoundary(boundary,nEdges)
rawIds = double(boundary.edgeIds(:));
rawValues = boundary.values(:);
if numel(rawIds) ~= numel(rawValues) || ...
        any(rawIds < 1) || any(rawIds > nEdges) || ...
        numel(unique(rawIds)) ~= numel(rawIds)
    error('tdgl:solvers:InvalidBoundaryCondition', ...
        'Boundary edge IDs must be unique, valid, and match values.');
end
[edgeIds,order] = sort(rawIds);
values = rawValues(order);
end

function [fixedNodes,fixedValues] = normalizedScalarBoundary( ...
        rawNodes,rawValues,components,nNodes)
rawNodes = double(rawNodes(:));
rawValues = rawValues(:);
if isempty(components)
    if ~isempty(rawNodes)
        error('tdgl:solvers:PhiWithoutConductor', ...
            'Scalar-potential constraints require a conducting region.');
    end
    fixedNodes = zeros(0,1);
    fixedValues = zeros(0,1);
    return;
end
conductingNodes = unique(vertcat(components{:}));
if numel(rawNodes) ~= numel(rawValues) || ...
        numel(unique(rawNodes)) ~= numel(rawNodes) || ...
        any(rawNodes < 1) || any(rawNodes > nNodes) || ...
        any(rawNodes ~= fix(rawNodes)) || ...
        ~all(ismember(rawNodes,conductingNodes))
    error('tdgl:solvers:InvalidScalarPotentialBoundary', ...
        ['Fixed phi node IDs must be unique conducting-node indices ', ...
         'and match FixedPhiValues.']);
end

fixedNodes = rawNodes;
fixedValues = rawValues;
for componentId = 1:numel(components)
    component = components{componentId};
    if isempty(intersect(component,fixedNodes))
        fixedNodes(end+1,1) = component(1); %#ok<AGROW>
        fixedValues(end+1,1) = 0; %#ok<AGROW>
    end
end
[fixedNodes,order] = sort(fixedNodes);
fixedValues = fixedValues(order);
end
