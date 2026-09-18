function [values, gradients] = lagrangeP1(barycentric, gradLambda)
%LAGRANGEP1 Evaluate first-order tetrahedral Lagrange basis functions.

if size(barycentric,2) ~= 4
    error('tdgl:elements:InvalidBarycentricCoordinates', ...
        'Barycentric coordinates must be Q-by-4.');
end
if any(abs(sum(barycentric,2)-1) > 1e-12)
    error('tdgl:elements:InvalidBarycentricCoordinates', ...
        'Each barycentric coordinate row must sum to one.');
end

values = barycentric;
if nargin > 1
    if ~isequal(size(gradLambda), [4 3])
        error('tdgl:elements:InvalidGradientShape', ...
            'gradLambda must be 4-by-3 for one tetrahedron.');
    end
    gradients = gradLambda;
else
    gradients = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
end
end
