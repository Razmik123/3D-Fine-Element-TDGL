function condition = boundaryCondition(name, field, type, tags, value)
%BOUNDARYCONDITION Create a declarative physical boundary condition.

arguments
    name (1,1) string
    field (1,1) string {mustBeMember(field,["A","psi","phi"])}
    type (1,1) string
    tags string
    value = []
end

condition = struct( ...
    'name',name, ...
    'field',field, ...
    'type',type, ...
    'tags',string(tags(:)), ...
    'value',value);
end
