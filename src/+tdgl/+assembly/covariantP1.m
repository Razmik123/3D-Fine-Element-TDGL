function matrix = covariantP1(mesh, edgePotential, coefficient)
%COVARIANTP1 Assemble the Hermitian covariant-gradient matrix.
%   MATRIX represents integral K (grad+iA)v* dot (grad-iA)u over the
%   GL-active mesh. A is represented by first-order Nedelec edge DOFs.

arguments
    mesh struct
    edgePotential (:,1) double
    coefficient = 1
end

nCells = size(mesh.cells,1);
nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
if numel(edgePotential) ~= nEdges
    error('tdgl:assembly:EdgeFieldSizeMismatch', ...
        'edgePotential must contain one value per global edge.');
end
coefficient = tdgl.assembly.cellCoefficient(coefficient,nCells,'coefficient');
[barycentric,weights] = tdgl.elements.tetraQuadrature(4);

rows = zeros(16*nCells,1);
columns = zeros(16*nCells,1);
values = complex(zeros(16*nCells,1));
cursor = 0;

for cellId = 1:nCells
    nodeIds = double(mesh.cells(cellId,:));
    edgeIds = double(mesh.topology.cellEdges(cellId,:));
    signs = double(mesh.topology.cellEdgeSigns(cellId,:));
    gradients = mesh.geometry.gradLambda(:,:,cellId);
    [edgeBasis,~] = tdgl.elements.nedelec1(barycentric,gradients);
    localEdgeDofs = signs(:).*edgePotential(edgeIds);
    localMatrix = complex(zeros(4));

    for q = 1:numel(weights)
        basisVectors = squeeze(edgeBasis(q,:,:)).';
        A = localEdgeDofs.'*basisVectors;
        lambda = barycentric(q,:);
        covariantTrial = gradients - 1i*(lambda.'*A);
        localMatrix = localMatrix + weights(q)*mesh.geometry.detJ(cellId) * ...
            coefficient(cellId) * (conj(covariantTrial)*covariantTrial.');
        % The expression above is equivalent to
        % (grad(lambda_i)+i*A*lambda_i) dot
        % (grad(lambda_j)-i*A*lambda_j).
    end

    indices = cursor+(1:16);
    [localRows,localColumns] = ndgrid(nodeIds,nodeIds);
    rows(indices) = localRows(:);
    columns(indices) = localColumns(:);
    values(indices) = localMatrix(:);
    cursor = cursor+16;
end

matrix = sparse(rows,columns,values,nNodes,nNodes);
end
