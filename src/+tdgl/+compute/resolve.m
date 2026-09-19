function backend = resolve(policy,matrix)
%RESOLVE Select a concrete execution backend for one linear system.

if nargin < 1 || isempty(policy)
    policy = tdgl.compute.execution("auto");
elseif ischar(policy) || isstring(policy)
    policy = tdgl.compute.execution(string(policy));
end
required = {'mode','minimumUnknowns','memoryFraction','allowFallback'};
if ~isstruct(policy) || ~all(isfield(policy,required))
    error('tdgl:compute:InvalidExecutionPolicy', ...
        'Execution must be created by tdgl.compute.execution.');
end

available = gpuDeviceCount("available") > 0;
estimatedBytes = sparseStorageEstimate(matrix);
useGPU = false;
reason = "CPU requested";
switch string(policy.mode)
    case "cpu"
        useGPU = false;
    case "gpu"
        if ~available
            error('tdgl:compute:GPUUnavailable', ...
                'GPU execution was requested but no supported GPU is available.');
        end
        device = gpuDevice;
        if ~device.SupportsDouble
            error('tdgl:compute:GPUDoubleUnsupported', ...
                'The selected GPU does not support double precision.');
        end
        useGPU = true;
        reason = "GPU explicitly requested";
    case "auto"
        if available
            device = gpuDevice;
            fits = estimatedBytes <= policy.memoryFraction*device.AvailableMemory;
            largeEnough = size(matrix,1) >= policy.minimumUnknowns;
            useGPU = device.SupportsDouble && fits && largeEnough;
            if useGPU
                reason = "automatic size and memory criteria passed";
            elseif ~largeEnough
                reason = "system below GPU crossover threshold";
            elseif ~fits
                reason = "estimated storage exceeds GPU memory budget";
            else
                reason = "GPU double precision unavailable";
            end
        else
            reason = "no supported GPU available";
        end
    otherwise
        error('tdgl:compute:InvalidMode','Unknown execution mode.');
end

backend = struct( ...
    'name',string(ternary(useGPU,"gpu","cpu")), ...
    'useGPU',useGPU, ...
    'requestedMode',string(policy.mode), ...
    'allowFallback',logical(policy.allowFallback), ...
    'estimatedMatrixBytes',estimatedBytes, ...
    'reason',reason);
if useGPU
    device = gpuDevice;
    backend.deviceName = string(device.Name);
    backend.availableMemory = double(device.AvailableMemory);
else
    backend.deviceName = "";
    backend.availableMemory = NaN;
end
end

function bytes = sparseStorageEstimate(matrix)
if issparse(matrix)
    bytesPerValue = 8;
    if ~isreal(matrix), bytesPerValue = 16; end
    bytes = double(nnz(matrix))*(bytesPerValue+8) + ...
        double(size(matrix,2)+1)*8;
else
    bytesPerValue = 8;
    if ~isreal(matrix), bytesPerValue = 16; end
    bytes = double(numel(matrix))*bytesPerValue;
end
% Sparse direct factors can be much larger than the input matrix.
bytes = 6*bytes;
end

function value = ternary(condition,yes,no)
if condition, value = yes; else, value = no; end
end
