function safe = jsonSafe(value)
%JSONSAFE Convert MATLAB values into portable JSON-compatible metadata.

if isa(value,'function_handle')
    safe = struct('kind',"function_handle",'definition',string(func2str(value)));
elseif isstruct(value)
    safe = value;
    names = fieldnames(value);
    for element = 1:numel(value)
        for index = 1:numel(names)
            safe(element).(names{index}) = ...
                tdgl.io.jsonSafe(value(element).(names{index}));
        end
    end
elseif iscell(value)
    safe = cellfun(@tdgl.io.jsonSafe,value,'UniformOutput',false);
elseif isnumeric(value) || islogical(value) || ischar(value) || isstring(value)
    safe = value;
elseif isa(value,'datetime') || isa(value,'duration')
    safe = string(value);
else
    safe = struct('kind',string(class(value)),'display',string(evalc('disp(value)')));
end
end
