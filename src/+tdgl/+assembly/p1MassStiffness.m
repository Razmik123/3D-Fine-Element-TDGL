function [massMatrix, stiffnessMatrix] = p1MassStiffness(mesh, massCoefficient, stiffnessCoefficient)
%P1MASSSTIFFNESS Assemble scalar P1 mass and diffusion matrices.

arguments
    mesh struct
    massCoefficient = 1
    stiffnessCoefficient = 1
end

nCells = size(mesh.cells,1);
nNodes = size(mesh.nodes,1);
massCoefficient = tdgl.assembly.cellCoefficient( ...
    massCoefficient, nCells, 'massCoefficient');
stiffnessCoefficient = tdgl.assembly.cellCoefficient( ...
    stiffnessCoefficient, nCells, 'stiffnessCoefficient');

nEntries = 16*nCells;
rows = zeros(nEntries,1);
columns = zeros(nEntries,1);
massValues = zeros(nEntries,1);
stiffnessValues = zeros(nEntries,1);
cursor = 0;
referenceMassPattern = ones(4) + eye(4);

for cellId = 1:nCells
    nodeIds = double(mesh.cells(cellId,:));
    volume = mesh.geometry.volume(cellId);
    gradients = mesh.geometry.gradLambda(:,:,cellId);
    localMass = massCoefficient(cellId) * volume/20 * referenceMassPattern;
    localStiffness = stiffnessCoefficient(cellId) * volume * ...
        (gradients * gradients.');
    indices = cursor + (1:16);
    [localRows, localColumns] = ndgrid(nodeIds,nodeIds);
    rows(indices) = localRows(:);
    columns(indices) = localColumns(:);
    massValues(indices) = localMass(:);
    stiffnessValues(indices) = localStiffness(:);
    cursor = cursor + 16;
end

massMatrix = sparse(rows, columns, massValues, nNodes, nNodes);
stiffnessMatrix = sparse(rows, columns, stiffnessValues, nNodes, nNodes);
end
