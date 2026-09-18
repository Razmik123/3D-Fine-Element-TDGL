function solution = solveMagnetostatic(mesh, model)
%SOLVEMAGNETOSTATIC Solve a mixed edge-element magnetostatic problem.
%   This verified foundation solves
%       curl(muInv curl(A)) = J
%   with prescribed tangential A on the complete exterior boundary and a
%   weak Coulomb constraint. It is the static Maxwell kernel used before
%   adding TDGL current and transient conductivity blocks.
%
%   MODEL fields:
%     muInv          scalar or one value per tetrahedron (default 1)
%     sourceCurrent  vector-field function or 1-by-3 vector (default zero)
%     boundary       output of tdgl.boundary.tangentialDirichlet

arguments
    mesh struct
    model struct
end

if ~isfield(model,'boundary')
    error('tdgl:solvers:MissingBoundaryCondition', ...
        'model.boundary must prescribe tangential A on the outer boundary.');
end
if ~isfield(model,'muInv'), model.muInv = 1; end
if ~isfield(model,'sourceCurrent'), model.sourceCurrent = [0 0 0]; end

[edgeMass, curlCurl] = tdgl.assembly.edgeMassCurlCurl(mesh,1,model.muInv);
loadVector = tdgl.assembly.edgeLoad(mesh,model.sourceCurrent);
[G,~,~] = tdgl.topology.incidenceMatrices(mesh);

rawFixedEdges = double(model.boundary.edgeIds(:));
rawFixedValues = model.boundary.values(:);
if numel(rawFixedEdges) ~= numel(rawFixedValues) || ...
        numel(unique(rawFixedEdges)) ~= numel(rawFixedEdges)
    error('tdgl:solvers:InvalidBoundaryCondition', ...
        'Boundary edge IDs and values must be unique and have equal length.');
end
[fixedEdges,order] = sort(rawFixedEdges);
fixedValues = rawFixedValues(order);

allBoundaryEdges = unique(double(mesh.topology.faceEdges( ...
    double(mesh.topology.boundaryFaceIds),:)));
if ~isequal(fixedEdges, allBoundaryEdges)
    error('tdgl:solvers:IncompleteOuterBoundary', ...
        ['The current magnetostatic kernel requires tangential A on every ', ...
         'exterior boundary edge.']);
end

nEdges = double(mesh.topology.nEdges);
freeEdges = setdiff((1:nEdges).', fixedEdges);
edgePotential = zeros(nEdges,1);
edgePotential(fixedEdges) = fixedValues;

boundaryNodes = unique(double(mesh.topology.faces( ...
    double(mesh.topology.boundaryFaceIds),:)));
interiorNodes = setdiff((1:size(mesh.nodes,1)).', boundaryNodes);
gaugeCoupling = edgeMass * G(:,interiorNodes);

rightHandSide = loadVector(freeEdges) - ...
    curlCurl(freeEdges,fixedEdges)*fixedValues;
constraintRightHandSide = -gaugeCoupling(fixedEdges,:).'*fixedValues;

if isempty(interiorNodes)
    systemMatrix = curlCurl(freeEdges,freeEdges);
    reducedSolution = systemMatrix \ rightHandSide;
    gaugeMultiplier = zeros(size(mesh.nodes,1),1);
else
    coupling = gaugeCoupling(freeEdges,:);
    systemMatrix = [curlCurl(freeEdges,freeEdges), coupling; ...
                    coupling.', sparse(numel(interiorNodes),numel(interiorNodes))];
    completeRightHandSide = [rightHandSide; constraintRightHandSide];
    reducedSolution = systemMatrix \ completeRightHandSide;
    gaugeMultiplier = zeros(size(mesh.nodes,1),1);
    gaugeMultiplier(interiorNodes) = reducedSolution(numel(freeEdges)+1:end);
    reducedSolution = reducedSolution(1:numel(freeEdges));
end

edgePotential(freeEdges) = reducedSolution;
algebraicResidual = curlCurl*edgePotential + ...
    edgeMass*G*gaugeMultiplier-loadVector;
freeResidual = norm(algebraicResidual(freeEdges)) / ...
    max(1,norm(loadVector(freeEdges)));
gaugeResidual = norm(G(:,interiorNodes).'*edgeMass*edgePotential) / ...
    max(1,norm(edgePotential));

solution = struct();
solution.edgePotential = edgePotential;
solution.gaugeMultiplier = gaugeMultiplier;
solution.cellMagneticFluxDensity = ...
    tdgl.post.cellMagneticFluxDensity(mesh,edgePotential);
solution.diagnostics = struct( ...
    'freeResidual',freeResidual, ...
    'gaugeResidual',gaugeResidual, ...
    'nFreeEdges',numel(freeEdges), ...
    'nGaugeDofs',numel(interiorNodes));
end
