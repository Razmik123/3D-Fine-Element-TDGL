function definition = terminal(name, tags, mode, excitation)
%TERMINAL Create a current, voltage, ground, floating, or probe terminal.

arguments
    name (1,1) string
    tags string
    mode (1,1) string {mustBeMember(mode,["current","voltage","ground","floating","probe"])}
    excitation = []
end

if mode == "ground" && isempty(excitation)
    excitation = 0;
end
if ismember(mode,["current","voltage"]) && isempty(excitation)
    error('tdgl:problem:MissingTerminalExcitation', ...
        'Current and voltage terminals require an excitation.');
end

definition = struct( ...
    'name',name, ...
    'tags',string(tags(:)), ...
    'mode',mode, ...
    'excitation',excitation);
end
