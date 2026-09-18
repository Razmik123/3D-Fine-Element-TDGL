function problem = compile(mesh, assignments, experiment)
%COMPILE Resolve materials, interfaces, boundaries, and terminals.
%   This function performs no PDE solve. It converts user-facing physical
%   definitions into explicit mesh entities and rejects ambiguous models.
%
%   ASSIGNMENTS is a struct array with fields regionId and material.
%   EXPERIMENT fields:
%     name, outerBoundaryTags, boundaryConditions, terminals,
%     interfaceModels, time, initialCondition.

arguments
    mesh struct
    assignments struct
    experiment struct
end

tdgl.mesh.validate(mesh);
regionIds = unique(double(mesh.regionIds(:))).';
assignedIds = arrayfun(@(item) double(item.regionId),assignments);
if numel(unique(assignedIds)) ~= numel(assignedIds)
    error('tdgl:problem:DuplicateMaterialAssignment', ...
        'Every material region must be assigned exactly once.');
end
if ~isequal(sort(regionIds),sort(assignedIds))
    error('tdgl:problem:IncompleteMaterialAssignment', ...
        'Material assignments must exactly cover mesh region IDs.');
end

materialsByRegion = cell(max(regionIds),1);
for index = 1:numel(assignments)
    material = assignments(index).material;
    required = {'kind','conductivity','relativePermeability','glActive'};
    if ~all(isfield(material,required))
        error('tdgl:problem:InvalidMaterial', ...
            'Region %d has an incomplete material record.',assignedIds(index));
    end
    materialsByRegion{assignedIds(index)} = material;
end

if ~isfield(experiment,'outerBoundaryTags') || ...
        isempty(experiment.outerBoundaryTags)
    error('tdgl:problem:MissingOuterBoundary', ...
        'The surrounding electromagnetic outer boundary must be explicit.');
end
outerFaces = tdgl.boundary.selectFaces(mesh,experiment.outerBoundaryTags);
if ~isequal(sort(outerFaces),sort(double(mesh.topology.boundaryFaceIds)))
    error('tdgl:problem:IncompleteOuterBoundary', ...
        'outerBoundaryTags must cover every exterior face exactly.');
end

if ~isfield(experiment,'boundaryConditions')
    experiment.boundaryConditions = struct([]);
end
compiledBoundary = compileBoundaryConditions(mesh,experiment.boundaryConditions);

if ~isfield(experiment,'terminals')
    experiment.terminals = struct([]);
end
compiledTerminals = compileTerminals(mesh,experiment.terminals);
validateTerminalExcitations(compiledTerminals);

if ~isfield(experiment,'interfaceModels')
    experiment.interfaceModels = struct([]);
end
compiledInterfaces = compileInterfaces( ...
    mesh,materialsByRegion,experiment.interfaceModels);

problem = struct();
problem.mesh = mesh;
problem.materialsByRegion = materialsByRegion;
problem.outerFaceIds = outerFaces;
problem.boundaryConditions = compiledBoundary;
problem.terminals = compiledTerminals;
problem.interfaces = compiledInterfaces;
problem.experiment = experiment;
problem.metadata = struct('compilerVersion',tdgl.version());
end

function compiled = compileBoundaryConditions(mesh,conditions)
compiled = conditions;
occupied = containers.Map('KeyType','char','ValueType','any');
for index = 1:numel(conditions)
    faceIds = tdgl.boundary.selectFaces(mesh,conditions(index).tags);
    nodeIds = unique(double(mesh.topology.faces(faceIds,:)));
    edgeIds = unique(double(mesh.topology.faceEdges(faceIds,:)));
    compiled(index).faceIds = faceIds;
    compiled(index).nodeIds = nodeIds;
    compiled(index).edgeIds = edgeIds;

    if ismember(conditions(index).type,["dirichlet","tangential-dirichlet"])
        key = char(conditions(index).field);
        if ~isKey(occupied,key), occupied(key) = []; end
        dofs = nodeIds;
        if conditions(index).field == "A", dofs = edgeIds; end
        overlap = intersect(occupied(key),dofs);
        if ~isempty(overlap)
            error('tdgl:problem:ConflictingBoundaryConditions', ...
                'Strong conditions for field %s overlap.',key);
        end
        occupied(key) = union(occupied(key),dofs);
    end
end
end

function compiled = compileTerminals(mesh,terminals)
compiled = terminals;
usedFaces = [];
names = strings(numel(terminals),1);
for index = 1:numel(terminals)
    names(index) = terminals(index).name;
    faceIds = tdgl.boundary.selectFaces(mesh,terminals(index).tags);
    if ~isempty(intersect(usedFaces,faceIds))
        error('tdgl:problem:OverlappingTerminals', ...
            'Terminal surfaces must not overlap.');
    end
    usedFaces = union(usedFaces,faceIds);
    compiled(index).faceIds = faceIds;
    compiled(index).nodeIds = unique(double(mesh.topology.faces(faceIds,:)));
    compiled(index).area = sum(mesh.faceAreas(faceIds));
end
if numel(unique(names)) ~= numel(names)
    error('tdgl:problem:DuplicateTerminalName', ...
        'Terminal names must be unique.');
end
end

function validateTerminalExcitations(terminals)
if isempty(terminals), return; end
modes = string({terminals.mode});
if any(modes == "voltage") && ~any(modes == "ground")
    error('tdgl:problem:MissingPotentialReference', ...
        'A voltage-driven experiment requires a ground terminal.');
end

currentIds = find(modes == "current");
if isempty(currentIds), return; end
values = zeros(numel(currentIds),1);
allConstant = true;
for index = 1:numel(currentIds)
    excitation = terminals(currentIds(index)).excitation;
    if isnumeric(excitation) && isscalar(excitation)
        values(index) = excitation;
    else
        allConstant = false;
    end
end
if allConstant && abs(sum(values)) > 1e-12*max(1,sum(abs(values)))
    error('tdgl:problem:UnbalancedTerminalCurrent', ...
        'Constant terminal currents must sum to zero.');
end
end

function compiled = compileInterfaces(mesh,materialsByRegion,models)
faceIds = double(mesh.interfaceFaceIds(:));
compiled = repmat(struct('faceIds',[],'regionIds',[],'electromagneticType', ...
    "transmission",'orderParameterType',"",'gammaB',0),0,1);
if isempty(faceIds), return; end

pairs = sort(double(mesh.faceRegionIds(faceIds,:)),2);
[uniquePairs,~,pairMap] = unique(pairs,'rows');
for pairId = 1:size(uniquePairs,1)
    pair = uniquePairs(pairId,:);
    matchingModel = find(arrayfun(@(item) ...
        isequal(sort(double(item.regionIds)),pair),models));
    glActive = [materialsByRegion{pair(1)}.glActive, ...
                materialsByRegion{pair(2)}.glActive];
    needsOrderLaw = any(glActive);
    if needsOrderLaw && numel(matchingModel) ~= 1
        error('tdgl:problem:MissingInterfaceModel', ...
            ['Interface [%d %d] touches a GL-active region and requires ', ...
             'exactly one explicit interface model.'],pair(1),pair(2));
    end
    record = struct('faceIds',faceIds(pairMap == pairId), ...
        'regionIds',pair,'electromagneticType',"transmission", ...
        'orderParameterType',"no-order-parameter",'gammaB',0);
    if needsOrderLaw
        record.orderParameterType = models(matchingModel).orderParameterType;
        record.gammaB = models(matchingModel).gammaB;
    end
    compiled(end+1,1) = record; %#ok<AGROW>
end
end
