function [values, curls] = nedelec1(barycentric, gradLambda)
%NEDELEC1 Evaluate lowest-order first-family tetrahedral edge functions.
%   Local edges are (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).
%   The local basis for oriented edge i->j is
%       N_ij = lambda_i grad(lambda_j) - lambda_j grad(lambda_i).
%   Its line integral along that local edge is one.

if size(barycentric,2) ~= 4 || ~isequal(size(gradLambda), [4 3])
    error('tdgl:elements:InvalidNedelecInput', ...
        'Expected Q-by-4 barycentric coordinates and 4-by-3 gradients.');
end

localEdges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nPoints = size(barycentric,1);
values = zeros(nPoints,3,6);
curls = zeros(6,3);

for edgeId = 1:6
    i = localEdges(edgeId,1);
    j = localEdges(edgeId,2);
    values(:,:,edgeId) = ...
        barycentric(:,i) * gradLambda(j,:) - ...
        barycentric(:,j) * gradLambda(i,:);
    curls(edgeId,:) = 2 * cross(gradLambda(i,:), gradLambda(j,:));
end
end
