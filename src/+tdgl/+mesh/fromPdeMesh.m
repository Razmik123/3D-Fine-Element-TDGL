function mesh = fromPdeMesh(pdeMesh, options)
%FROMPDEMESH Convert a PDE Toolbox FEMesh to the canonical TDGL mesh.
%   Boundary semantics are intentionally not guessed. Use BoundaryTagger or
%   tdgl.mesh.tagBoundary after conversion to assign physical face roles.

arguments
    pdeMesh
    options.RegionIds = []
    options.BoundaryTagger = []
end

if ~isprop(pdeMesh,'Nodes') || ~isprop(pdeMesh,'Elements')
    error('tdgl:mesh:InvalidPdeMesh', ...
        'Input must expose PDE Toolbox Nodes and Elements properties.');
end

regionIds = options.RegionIds;
if isempty(regionIds) && isprop(pdeMesh,'ElementIDToRegionID')
    candidate = pdeMesh.ElementIDToRegionID;
    if ~isempty(candidate)
        regionIds = candidate;
    end
end
if isempty(regionIds)
    regionIds = ones(size(pdeMesh.Elements,2),1);
end

mesh = tdgl.mesh.fromArrays(pdeMesh.Nodes,pdeMesh.Elements, ...
    'RegionIds',regionIds);
mesh.metadata.generator = "tdgl.mesh.fromPdeMesh";
if ~isempty(options.BoundaryTagger)
    if ~isa(options.BoundaryTagger,'function_handle')
        error('tdgl:mesh:InvalidBoundaryTagger', ...
            'BoundaryTagger must be a function handle.');
    end
    mesh = options.BoundaryTagger(mesh);
end
end
