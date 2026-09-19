function [solution,information] = solveLinear(matrix,rightHandSide,policy)
%SOLVELINEAR Solve on CPU or GPU according to an execution policy.

arguments
    matrix
    rightHandSide
    policy = []
end
backend = tdgl.compute.resolve(policy,matrix);
started = tic;
fellBack = false;
failure = "";
if backend.useGPU
    try
        solution = gather(gpuArray(matrix)\gpuArray(rightHandSide));
    catch exception
        if ~backend.allowFallback
            rethrow(exception);
        end
        failure = string(exception.message);
        solution = matrix\rightHandSide;
        fellBack = true;
    end
else
    solution = matrix\rightHandSide;
end

information = backend;
information.elapsedSeconds = toc(started);
information.fellBack = fellBack;
information.failure = failure;
if fellBack
    information.name = "cpu";
    information.reason = "GPU solve failed; CPU fallback used";
end
end
