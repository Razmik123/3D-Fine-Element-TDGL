function voltage = terminalElectrochemicalVoltage(mesh,faceIdsA,faceIdsB, ...
        orderPrevious,orderCurrent,scalarPotential,dt)
%TERMINALELECTROCHEMICALVOLTAGE Difference of surface-averaged mu_gi.
%   The returned value is mean_A(phi+d(arg psi)/dt) minus the corresponding
%   mean on B. Nodes where the order-parameter phase is undefined are
%   excluded. This observable is appropriate for superconducting probes.

electrochemical = tdgl.observe.electrochemicalPotential( ...
    orderPrevious,orderCurrent,scalarPotential,dt);
valueA = finiteSurfaceAverage(mesh,faceIdsA,electrochemical);
valueB = finiteSurfaceAverage(mesh,faceIdsB,electrochemical);
voltage = valueA-valueB;
end

function value = finiteSurfaceAverage(mesh,faceIds,nodalField)
faceIds = double(faceIds(:));
faces = double(mesh.topology.faces(faceIds,:));
areas = mesh.faceAreas(faceIds);
faceValues = zeros(numel(faceIds),1);
faceWeights = zeros(numel(faceIds),1);
for index = 1:numel(faceIds)
    values = nodalField(faces(index,:));
    valid = isfinite(values);
    if any(valid)
        faceValues(index) = mean(values(valid));
        faceWeights(index) = areas(index)*nnz(valid)/3;
    end
end
if sum(faceWeights) == 0
    error('tdgl:observe:UndefinedTerminalPhase', ...
        'The order-parameter phase is undefined on the complete terminal.');
end
value = sum(faceWeights.*faceValues)/sum(faceWeights);
end
