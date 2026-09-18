function results = run_all_tests
%RUN_ALL_TESTS Run the complete local MATLAB verification suite.

repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repositoryRoot);
startup;
results = runtests(fullfile(repositoryRoot,'tests','unit'), ...
    'IncludeSubfolders',true);
assertSuccess(results);
end
