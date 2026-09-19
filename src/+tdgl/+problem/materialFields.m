function fields = materialFields(problem)
%MATERIALFIELDS Expand region materials into cellwise coefficient arrays.
%   GLActive is deliberately returned separately. A caller must choose a
%   no-proximity restricted domain or a coefficient-based proximity model;
%   this routine does not silently activate GL physics in a normal region.

mesh = problem.mesh;
nCells = size(mesh.cells,1);
fields = struct( ...
    'conductivity',zeros(nCells,1), ...
    'muInv',zeros(nCells,1), ...
    'glActive',false(nCells,1), ...
    'a',zeros(nCells,1), ...
    'b',zeros(nCells,1), ...
    'K',zeros(nCells,1), ...
    'u',zeros(nCells,1));

for regionId = unique(double(mesh.regionIds(:))).'
    cells = double(mesh.regionIds) == regionId;
    material = problem.materialsByRegion{regionId};
    fields.conductivity(cells) = material.conductivity;
    fields.muInv(cells) = 1/material.relativePermeability;
    fields.glActive(cells) = material.glActive;
    fields.a(cells) = material.a;
    fields.b(cells) = material.b;
    fields.K(cells) = material.K;
    fields.u(cells) = material.u;
end
end
