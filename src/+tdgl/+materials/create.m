function material = create(kind, options)
%CREATE Define one dimensionless material record.
%   The defaults describe reference materials for verification. Real runs
%   must replace them with calibrated values and record the scale mapping.

arguments
    kind (1,1) string {mustBeMember(kind,["superconductor","normal-metal","insulator","vacuum"])}
    options.Name (1,1) string = ""
    options.Conductivity = []
    options.RelativePermeability (1,1) double {mustBePositive} = 1
    options.GLActive = []
    options.A (1,1) double = NaN
    options.B (1,1) double {mustBePositive} = 1
    options.K (1,1) double {mustBePositive} = 1
    options.U (1,1) double {mustBePositive} = 1
end

if options.Name == ""
    options.Name = kind;
end

switch kind
    case "superconductor"
        defaultConductivity = 1;
        defaultGLActive = true;
        defaultA = 1;
    case "normal-metal"
        defaultConductivity = 1;
        defaultGLActive = false;
        defaultA = -1;
    otherwise
        defaultConductivity = 0;
        defaultGLActive = false;
        defaultA = -1;
end

if isempty(options.Conductivity)
    options.Conductivity = defaultConductivity;
elseif ~isnumeric(options.Conductivity) || ~isscalar(options.Conductivity) || ...
        ~isfinite(options.Conductivity) || options.Conductivity < 0
    error('tdgl:materials:InvalidConductivity', ...
        'Conductivity must be a finite nonnegative scalar.');
end
if isnan(options.A)
    options.A = defaultA;
end
if isempty(options.GLActive)
    options.GLActive = defaultGLActive;
elseif ~islogical(options.GLActive) || ~isscalar(options.GLActive)
    error('tdgl:materials:InvalidGLActive', ...
        'GLActive must be a logical scalar.');
end

material = struct( ...
    'name',options.Name, ...
    'kind',kind, ...
    'conductivity',options.Conductivity, ...
    'relativePermeability',options.RelativePermeability, ...
    'glActive',options.GLActive, ...
    'a',options.A, ...
    'b',options.B, ...
    'K',options.K, ...
    'u',options.U);
end
