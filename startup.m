function startup
%STARTUP Add the TDGL source tree to the MATLAB path.
%   Run this function once after opening the repository in MATLAB.

repositoryRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot, 'src'));
end
