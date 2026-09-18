function energy = orderParameterEnergy(mesh,orderParameter,model)
%ORDERPARAMETERENERGY Compute kinetic and local GL energy components.

[massA,~] = tdgl.assembly.p1MassStiffness(mesh,model.a,0);
covariant = tdgl.assembly.covariantP1(mesh,model.edgePotential,model.K);
[barycentric,weights] = tdgl.elements.tetraQuadrature(4);
quartic = 0;
b = tdgl.assembly.cellCoefficient(model.b,size(mesh.cells,1),'b');

for cellId = 1:size(mesh.cells,1)
    nodeIds = double(mesh.cells(cellId,:));
    values = barycentric*orderParameter(nodeIds);
    quartic = quartic + b(cellId)*mesh.geometry.detJ(cellId) * ...
        sum(weights.*abs(values).^4)/2;
end

energy = struct();
energy.kinetic = real(orderParameter'*covariant*orderParameter);
energy.quadratic = -real(orderParameter'*massA*orderParameter);
energy.quartic = quartic;
energy.total = energy.kinetic+energy.quadratic+energy.quartic;
end
