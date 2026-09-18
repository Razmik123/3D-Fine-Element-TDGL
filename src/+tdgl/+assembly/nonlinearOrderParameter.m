function [residual,jacobian] = nonlinearOrderParameter(mesh,orderParameter,coefficient)
%NONLINEARORDERPARAMETER Assemble b*|psi|^2*psi and its real Jacobian.
%   JACOBIAN acts on [delta(real(psi)); delta(imag(psi))].

arguments
    mesh struct
    orderParameter (:,1) double
    coefficient = 1
end

nNodes = size(mesh.nodes,1);
nCells = size(mesh.cells,1);
if numel(orderParameter) ~= nNodes
    error('tdgl:assembly:NodalFieldSizeMismatch', ...
        'orderParameter must contain one value per mesh node.');
end
coefficient = tdgl.assembly.cellCoefficient(coefficient,nCells,'coefficient');
[barycentric,weights] = tdgl.elements.tetraQuadrature(4);

residual = complex(zeros(nNodes,1));
rows = zeros(16*nCells,1);
columns = zeros(16*nCells,1);
jxx = zeros(16*nCells,1);
jxy = zeros(16*nCells,1);
jyx = zeros(16*nCells,1);
jyy = zeros(16*nCells,1);
cursor = 0;

for cellId = 1:nCells
    nodeIds = double(mesh.cells(cellId,:));
    localPsi = orderParameter(nodeIds);
    localResidual = complex(zeros(4,1));
    localJxx = zeros(4);
    localJxy = zeros(4);
    localJyx = zeros(4);
    localJyy = zeros(4);
    scale = coefficient(cellId)*mesh.geometry.detJ(cellId);

    for q = 1:numel(weights)
        shape = barycentric(q,:).';
        value = shape.'*localPsi;
        x = real(value);
        y = imag(value);
        weight = weights(q)*scale;
        localResidual = localResidual + weight*shape*(abs(value)^2*value);
        shapeProduct = shape*shape.';
        localJxx = localJxx + weight*(3*x^2+y^2)*shapeProduct;
        localJxy = localJxy + weight*(2*x*y)*shapeProduct;
        localJyx = localJyx + weight*(2*x*y)*shapeProduct;
        localJyy = localJyy + weight*(x^2+3*y^2)*shapeProduct;
    end
    residual(nodeIds) = residual(nodeIds)+localResidual;

    indices = cursor+(1:16);
    [localRows,localColumns] = ndgrid(nodeIds,nodeIds);
    rows(indices) = localRows(:);
    columns(indices) = localColumns(:);
    jxx(indices) = localJxx(:);
    jxy(indices) = localJxy(:);
    jyx(indices) = localJyx(:);
    jyy(indices) = localJyy(:);
    cursor = cursor+16;
end

Jxx = sparse(rows,columns,jxx,nNodes,nNodes);
Jxy = sparse(rows,columns,jxy,nNodes,nNodes);
Jyx = sparse(rows,columns,jyx,nNodes,nNodes);
Jyy = sparse(rows,columns,jyy,nNodes,nNodes);
jacobian = [Jxx Jxy; Jyx Jyy];
end
