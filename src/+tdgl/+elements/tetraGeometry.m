function geometry = tetraGeometry(nodes, cells)
%TETRAGEOMETRY Compute affine maps and barycentric gradients per cell.
%   The reference tetrahedron has vertices (0,0,0), (1,0,0), (0,1,0),
%   and (0,0,1). GradLambda(:,:,K) contains physical gradients of the four
%   barycentric basis functions for cell K as row vectors.

nCells = size(cells,1);
jacobian = zeros(3,3,nCells);
inverseJacobian = zeros(3,3,nCells);
detJ = zeros(nCells,1);
gradLambda = zeros(4,3,nCells);
referenceGradients = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];

for cellId = 1:nCells
    vertexIds = double(cells(cellId,:));
    coordinates = nodes(vertexIds,:);
    J = [coordinates(2,:)-coordinates(1,:); ...
         coordinates(3,:)-coordinates(1,:); ...
         coordinates(4,:)-coordinates(1,:)].';
    determinant = det(J);
    jacobian(:,:,cellId) = J;
    inverseJacobian(:,:,cellId) = inv(J);
    detJ(cellId) = determinant;
    gradLambda(:,:,cellId) = referenceGradients / J;
end

geometry = struct();
geometry.jacobian = jacobian;
geometry.inverseJacobian = inverseJacobian;
geometry.detJ = detJ;
geometry.volume = detJ / 6;
geometry.gradLambda = gradLambda;
end
