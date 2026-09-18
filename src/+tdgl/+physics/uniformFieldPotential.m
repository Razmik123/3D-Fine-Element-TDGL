function potential = uniformFieldPotential(magneticFluxDensity, origin)
%UNIFORMFIELDPOTENTIAL Symmetric-gauge vector potential for uniform B.
%   A(x) = 1/2 B cross (x-origin), so curl(A) = B.

arguments
    magneticFluxDensity (1,3) double
    origin (1,3) double = [0 0 0]
end

B = magneticFluxDensity;
x0 = origin;
potential = @(points) 0.5 * cross(repmat(B,size(points,1),1), ...
    points-x0, 2);
end
