function space = scalarPotentialSpace(mesh,conductingCellMask,terminals,options)
%SCALARPOTENTIALSPACE Build constrained nodal coordinates for physical phi.
%   Ordinary conducting nodes receive individual unknowns. Current and
%   floating metal terminals receive one equipotential unknown each.
%   Voltage/ground terminals enter the fixed offset. One harmless reference
%   is removed from every otherwise floating conducting component.

arguments
    mesh struct
    conductingCellMask (:,1) logical
    terminals struct = struct([])
    options.FixedNodeIds = []
    options.FixedValues = []
end

nNodes = size(mesh.nodes,1);
nEdges = double(mesh.topology.nEdges);
components = tdgl.topology.conductingComponents(mesh,conductingCellMask);
if isempty(components)
    if ~isempty(terminals) || ~isempty(options.FixedNodeIds)
        error('tdgl:boundary:PhiConstraintWithoutConductor', ...
            'Scalar-potential constraints require a conducting region.');
    end
    space = emptySpace(nNodes,nEdges);
    return;
end
conductingNodes = unique(vertcat(components{:}));

[fixedNodes,fixedValues] = validateFixed( ...
    options.FixedNodeIds,options.FixedValues,conductingNodes,nNodes);
fixedMask = false(nNodes,1);
fixedMask(fixedNodes) = true;
offset = zeros(nNodes,1);
offset(fixedNodes) = fixedValues;

nTerminals = numel(terminals);
terminalNodes = cell(nTerminals,1);
terminalModes = strings(nTerminals,1);
terminalTargets = zeros(nTerminals,1);
terminalNames = strings(nTerminals,1);
terminalComponent = zeros(nTerminals,1);
occupied = false(nNodes,1);
groupTerminalIds = zeros(0,1);

for terminalId = 1:nTerminals
    required = {'name','mode','excitation','nodeIds'};
    if ~all(isfield(terminals(terminalId),required))
        error('tdgl:boundary:IncompleteTerminal', ...
            'Every terminal requires name, mode, excitation, and nodeIds.');
    end
    nodes = unique(double(terminals(terminalId).nodeIds(:)));
    mode = string(terminals(terminalId).mode);
    if isempty(nodes) || any(nodes < 1) || any(nodes > nNodes) || ...
            ~all(ismember(nodes,conductingNodes))
        error('tdgl:boundary:InvalidTerminalNodes', ...
            'Terminal nodes must be valid nodes of a conducting region.');
    end
    if any(occupied(nodes))
        error('tdgl:boundary:OverlappingTerminals', ...
            'Electrical terminal node sets must not overlap.');
    end
    componentIds = find(cellfun(@(item) all(ismember(nodes,item)),components));
    if numel(componentIds) ~= 1
        error('tdgl:boundary:TerminalSpansComponents', ...
            'A terminal must belong to exactly one conducting component.');
    end
    terminalNodes{terminalId} = nodes;
    terminalModes(terminalId) = mode;
    terminalNames(terminalId) = string(terminals(terminalId).name);
    terminalComponent(terminalId) = componentIds;
    occupied(nodes) = true;

    switch mode
        case {"voltage","ground"}
            value = scalarExcitation(terminals(terminalId).excitation,mode);
            if any(fixedMask(nodes) & abs(offset(nodes)-value) > 100*eps)
                error('tdgl:boundary:ConflictingTerminalPotential', ...
                    'A terminal voltage conflicts with a fixed nodal value.');
            end
            fixedMask(nodes) = true;
            offset(nodes) = value;
        case {"current","floating"}
            if mode == "current"
                terminalTargets(terminalId) = scalarExcitation( ...
                    terminals(terminalId).excitation,mode);
            end
            groupTerminalIds(end+1,1) = terminalId; %#ok<AGROW>
        case "probe"
            occupied(nodes) = false;
        otherwise
            error('tdgl:boundary:UnsupportedTerminalMode', ...
                'Unsupported terminal mode "%s".',mode);
    end
end

referenceTerminal = false(nTerminals,1);
for componentId = 1:numel(components)
    component = components{componentId};
    componentTerminals = find(terminalComponent == componentId & ...
        ismember(terminalModes,["current","floating"]));
    hasFixedPotential = any(fixedMask(component));
    hasVoltageTerminal = any(terminalComponent == componentId & ...
        ismember(terminalModes,["voltage","ground"]));
    currentTerminals = find(terminalComponent == componentId & ...
        terminalModes == "current");
    if ~hasVoltageTerminal && ~isempty(currentTerminals)
        imbalance = sum(terminalTargets(currentTerminals));
        if abs(imbalance) > 1e-12*max(1,sum(abs(terminalTargets(currentTerminals))))
            error('tdgl:boundary:UnbalancedComponentCurrent', ...
                'Prescribed currents must balance on each isolated component.');
        end
    end
    if ~hasFixedPotential
        if ~isempty(componentTerminals)
            referenceId = componentTerminals(1);
            referenceTerminal(referenceId) = true;
            nodes = terminalNodes{referenceId};
            fixedMask(nodes) = true;
            offset(nodes) = 0;
        else
            fixedMask(component(1)) = true;
            offset(component(1)) = 0;
        end
    end
end

activeGroupIds = groupTerminalIds(~referenceTerminal(groupTerminalIds));
groupNodes = unique(vertcat(terminalNodes{groupTerminalIds}));
ordinaryNodes = setdiff(conductingNodes,union(find(fixedMask),groupNodes));
nOrdinary = numel(ordinaryNodes);
nGroups = numel(activeGroupIds);
rows = [ordinaryNodes; vertcat(terminalNodes{activeGroupIds})];
columns = [(1:nOrdinary).'; groupColumns(activeGroupIds,terminalNodes,nOrdinary)];
basis = sparse(rows,columns,1,nNodes,nOrdinary+nGroups);
targetDivergence = [zeros(nOrdinary,1);terminalTargets(activeGroupIds)];

[G,~,~] = tdgl.topology.incidenceMatrices(mesh);
terminalTestColumns = cell(1,nTerminals);
for terminalId = 1:nTerminals
    if terminalModes(terminalId) ~= "probe"
        indicator = sparse(terminalNodes{terminalId},1,1,nNodes,1);
        terminalTestColumns{terminalId} = G*indicator;
    else
        terminalTestColumns{terminalId} = sparse(nEdges,1);
    end
end
if isempty(terminalTestColumns)
    terminalTests = sparse(nEdges,0);
else
    terminalTests = horzcat(terminalTestColumns{:});
end

space = struct();
space.basis = basis;
space.offset = offset;
space.targetDivergence = targetDivergence;
space.conductingNodes = conductingNodes;
space.fixedNodeIds = find(fixedMask);
space.ordinaryNodeIds = ordinaryNodes;
space.terminal = struct( ...
    'names',terminalNames, ...
    'modes',terminalModes, ...
    'nodeIds',{terminalNodes}, ...
    'componentIds',terminalComponent, ...
    'reference',referenceTerminal, ...
    'tests',terminalTests);
end

function space = emptySpace(nNodes,nEdges)
space = struct( ...
    'basis',sparse(nNodes,0), ...
    'offset',zeros(nNodes,1), ...
    'targetDivergence',zeros(0,1), ...
    'conductingNodes',zeros(0,1), ...
    'fixedNodeIds',zeros(0,1), ...
    'ordinaryNodeIds',zeros(0,1), ...
    'terminal',struct('names',strings(0,1),'modes',strings(0,1), ...
        'nodeIds',{cell(0,1)},'componentIds',zeros(0,1), ...
        'reference',false(0,1),'tests',sparse(nEdges,0)));
end

function [nodes,values] = validateFixed(nodes,values,conductingNodes,nNodes)
nodes = double(nodes(:));
values = values(:);
if numel(nodes) ~= numel(values) || numel(unique(nodes)) ~= numel(nodes) || ...
        any(nodes < 1) || any(nodes > nNodes) || any(nodes ~= fix(nodes)) || ...
        ~all(ismember(nodes,conductingNodes))
    error('tdgl:boundary:InvalidFixedPotential', ...
        'Fixed phi data require unique conducting-node IDs and matching values.');
end
end

function value = scalarExcitation(excitation,mode)
if ~isnumeric(excitation) || ~isscalar(excitation) || ...
        ~isreal(excitation) || ~isfinite(excitation)
    error('tdgl:boundary:UnevaluatedTerminalExcitation', ...
        'The %s terminal excitation must be a finite real scalar.',mode);
end
value = double(excitation);
end

function columns = groupColumns(groupIds,terminalNodes,offset)
columns = zeros(sum(cellfun(@numel,terminalNodes(groupIds))),1);
cursor = 0;
for groupIndex = 1:numel(groupIds)
    count = numel(terminalNodes{groupIds(groupIndex)});
    columns(cursor+(1:count)) = offset+groupIndex;
    cursor = cursor+count;
end
end
