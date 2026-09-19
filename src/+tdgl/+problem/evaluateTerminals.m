function evaluated = evaluateTerminals(terminals,time,state)
%EVALUATETERMINALS Evaluate time-dependent terminal excitations.
%   Excitation functions use the explicit contract VALUE = F(TIME,STATE).

arguments
    terminals struct
    time (1,1) double
    state struct
end
evaluated = terminals;
for index = 1:numel(terminals)
    excitation = terminals(index).excitation;
    if isa(excitation,'function_handle')
        value = excitation(time,state);
        if ~isnumeric(value) || ~isscalar(value) || ...
                ~isreal(value) || ~isfinite(value)
            error('tdgl:problem:InvalidTerminalExcitationValue', ...
                'Terminal excitation functions must return a finite real scalar.');
        end
        evaluated(index).excitation = double(value);
    end
end
end
