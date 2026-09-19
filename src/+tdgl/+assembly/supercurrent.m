function [phaseLoad,densityMass] = supercurrent(mesh,orderParameter,coefficient)
%SUPERCURRENT Assemble the two gauge-current contributions.
%   For j_s = K*Im(conj(psi)*grad(psi)) - K*|psi|^2*A, PHASELOAD
%   represents the first term tested against Nedelec functions and
%   DENSITYMASS represents integral K*|psi|^2*N_i dot N_j.

arguments
    mesh struct
    orderParameter (:,1) double
    coefficient = 1
end

nNodes = size(mesh.nodes,1);
nCells = size(mesh.cells,1);
nEdges = double(mesh.topology.nEdges);
if numel(orderParameter) ~= nNodes
    error('tdgl:assembly:NodalFieldSizeMismatch', ...
        'orderParameter must contain one value per mesh node.');
end
coefficient = tdgl.assembly.cellCoefficient(coefficient,nCells,'coefficient');
[barycentric,weights] = tdgl.elements.tetraQuadrature(4);

loadRows = zeros(6*nCells,1);
loadValues = zeros(6*nCells,1);
matrixRows = zeros(36*nCells,1);
matrixColumns = zeros(36*nCells,1);
matrixValues = zeros(36*nCells,1);
loadCursor = 0;
matrixCursor = 0;

for cellId = 1:nCells
    nodeIds = double(mesh.cells(cellId,:));
    edgeIds = double(mesh.topology.cellEdges(cellId,:));
    signs = double(mesh.topology.cellEdgeSigns(cellId,:));
    gradients = mesh.geometry.gradLambda(:,:,cellId);
    [edgeBasis,~] = tdgl.elements.nedelec1(barycentric,gradients);
    edgeBasis = edgeBasis.*reshape(signs,1,1,6);
    localPsi = orderParameter(nodeIds);
    gradientPsi = localPsi.'*gradients;
    localLoad = zeros(6,1);
    localMass = zeros(6);
    scale = coefficient(cellId)*mesh.geometry.detJ(cellId);

    for q = 1:numel(weights)
        psi = barycentric(q,:)*localPsi;
        phaseCurrent = imag(conj(psi)*gradientPsi);
        density = abs(psi)^2;
        vectors = squeeze(edgeBasis(q,:,:)).';
        weight = weights(q)*scale;
        localLoad = localLoad + weight*(vectors*phaseCurrent.');
        localMass = localMass + weight*density*(vectors*vectors.');
    end

    loadIndices = loadCursor+(1:6);
    loadRows(loadIndices) = edgeIds;
    loadValues(loadIndices) = localLoad;
    loadCursor = loadCursor+6;

    matrixIndices = matrixCursor+(1:36);
    [localRows,localColumns] = ndgrid(edgeIds,edgeIds);
    matrixRows(matrixIndices) = localRows(:);
    matrixColumns(matrixIndices) = localColumns(:);
    matrixValues(matrixIndices) = localMass(:);
    matrixCursor = matrixCursor+36;
end

phaseLoad = accumarray(loadRows,loadValues,[nEdges 1],@sum,0);
densityMass = sparse(matrixRows,matrixColumns,matrixValues,nEdges,nEdges);
end
