function [barycentric, weights] = tetraQuadrature(order)
%TETRAQUADRATURE Symmetric quadrature on the reference tetrahedron.
%   Weights integrate over the reference tetrahedron and therefore sum to
%   1/6. Orders 1, 2, and 4 are currently provided.

arguments
    order (1,1) double {mustBeInteger, mustBePositive} = 2
end

switch order
    case 1
        barycentric = [1 1 1 1] / 4;
        weights = 1/6;
    case 2
        a = 0.5854101966249685;
        b = 0.1381966011250105;
        barycentric = [ ...
            a b b b; ...
            b a b b; ...
            b b a b; ...
            b b b a];
        weights = repmat(1/24, 4, 1);
    case 4
        a = 11/14;
        b = 1/14;
        c = 0.3994035761667992;
        d = 0.1005964238332008;
        barycentric = [ ...
            1/4 1/4 1/4 1/4; ...
            a b b b; b a b b; b b a b; b b b a; ...
            c c d d; c d c d; c d d c; ...
            d c c d; d c d c; d d c c];
        weights = [ ...
            -0.0131555555555556; ...
            repmat(0.00762222222222222,4,1); ...
            repmat(0.0248888888888889,6,1)];
    otherwise
        error('tdgl:elements:UnsupportedQuadratureOrder', ...
            'Tetrahedral quadrature order %d is not implemented.', order);
end
end
