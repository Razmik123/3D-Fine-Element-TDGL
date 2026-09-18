function value = electrochemicalPotential(orderPrevious, orderCurrent, scalarPotential, dt)
%ELECTROCHEMICALPOTENTIAL Compute phi + time derivative of order phase.
%   The phase increment is evaluated without an explicit unwrap using the
%   argument of conj(psi_n)*psi_{n+1}. Values near a vortex core are marked
%   NaN because phase is undefined there.

arguments
    orderPrevious (:,1) double
    orderCurrent (:,1) double
    scalarPotential (:,1) double
    dt (1,1) double {mustBePositive}
end
if ~isequal(numel(orderPrevious),numel(orderCurrent),numel(scalarPotential))
    error('tdgl:observe:FieldSizeMismatch', ...
        'Order parameter and scalar potential arrays must have equal size.');
end

phaseRate = angle(conj(orderPrevious).*orderCurrent)/dt;
value = scalarPotential + phaseRate;
undefined = min(abs(orderPrevious),abs(orderCurrent)) < 100*eps;
value(undefined) = NaN;
end
