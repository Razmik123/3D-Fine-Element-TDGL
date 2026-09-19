function record = provenance
%PROVENANCE Capture code, MATLAB, CPU, and GPU identity for a run.

classFile = mfilename('fullpath');
repositoryRoot = fileparts(fileparts(fileparts(fileparts(classFile))));
[status,commit] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
if status ~= 0, commit = "unavailable"; end
record = struct();
record.createdUtc = string(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
record.gitCommit = strtrim(string(commit));
record.solver = tdgl.version();
record.computeEnvironment = tdgl.compute.environment();
record.computer = string(computer);
end
