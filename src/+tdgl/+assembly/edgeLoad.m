function loadVector = edgeLoad(mesh, vectorField)
%EDGELOAD Assemble integral of a vector field against Nedelec test functions.
%   VECTORFIELD is a function handle accepting Q-by-3 coordinates and
%   returning Q-by-3 vector values, or a constant 1-by-3 vector.

nCells = size(mesh.cells,1);
nEdges = double(mesh.topology.nEdges);
[barycentric, referenceWeights] = tdgl.elements.tetraQuadrature(2);
rows = zeros(6*nCells,1);
values = zeros(6*nCells,1);
cursor = 0;

for cellId = 1:nCells
    nodeIds = double(mesh.cells(cellId,:));
    coordinates = barycentric * mesh.nodes(nodeIds,:);
    if isa(vectorField, 'function_handle')
        fieldValues = vectorField(coordinates);
    else
        fieldValues = vectorField;
    end
    if isequal(size(fieldValues), [1 3])
        fieldValues = repmat(fieldValues, size(coordinates,1), 1);
    end
    if ~isequal(size(fieldValues), size(coordinates))
        error('tdgl:assembly:InvalidVectorField', ...
            'Vector field must return one three-component row per point.');
    end

    gradients = mesh.geometry.gradLambda(:,:,cellId);
    [basis, ~] = tdgl.elements.nedelec1(barycentric, gradients);
    signs = double(mesh.topology.cellEdgeSigns(cellId,:));
    basis = basis .* reshape(signs,1,1,6);
    localLoad = zeros(6,1);
    for q = 1:numel(referenceWeights)
        vectors = squeeze(basis(q,:,:)).';
        localLoad = localLoad + referenceWeights(q)*mesh.geometry.detJ(cellId) * ...
            (vectors*fieldValues(q,:).');
    end

    indices = cursor + (1:6);
    rows(indices) = double(mesh.topology.cellEdges(cellId,:));
    values(indices) = localLoad;
    cursor = cursor + 6;
end

loadVector = accumarray(rows, values, [nEdges 1], @sum, 0);
end
