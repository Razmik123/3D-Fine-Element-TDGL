function policy = execution(mode,options)
%EXECUTION Create an explicit CPU/GPU execution policy.
%   AUTO selects the GPU only for sufficiently large systems and when the
%   estimated sparse storage fits comfortably in currently available GPU
%   memory. Explicit GPU mode fails if no suitable device is available.

arguments
    mode (1,1) string {mustBeMember(mode,["auto","cpu","gpu"])} = "auto"
    options.MinimumUnknowns (1,1) double {mustBeInteger,mustBeNonnegative} = 20000
    options.MemoryFraction (1,1) double {mustBePositive,mustBeLessThanOrEqual( ...
        options.MemoryFraction,0.8)} = 0.35
    options.AllowFallback (1,1) logical = true
end

policy = struct( ...
    'mode',mode, ...
    'minimumUnknowns',options.MinimumUnknowns, ...
    'memoryFraction',options.MemoryFraction, ...
    'allowFallback',options.AllowFallback);
end
