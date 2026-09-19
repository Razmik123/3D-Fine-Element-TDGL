function tests = testComputeBackend
%TESTCOMPUTEBACKEND Verify CPU, automatic, and optional GPU linear solves.
tests = functiontests(localfunctions);
end

function setupOnce(~)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
startup;
end

function testExplicitCpuSolve(testCase)
matrix = sparse([4 -1 0;-1 4 -1;0 -1 3]);
rhs = [1;2;3];
[solution,information] = tdgl.compute.solveLinear( ...
    matrix,rhs,tdgl.compute.execution("cpu"));
verifyLessThan(testCase,norm(matrix*solution-rhs),1e-13);
verifyEqual(testCase,information.name,"cpu");
verifyFalse(testCase,information.fellBack);
end

function testAutoKeepsSmallSystemOnCpu(testCase)
matrix = speye(5);
backend = tdgl.compute.resolve( ...
    tdgl.compute.execution("auto",'MinimumUnknowns',100),matrix);
verifyEqual(testCase,backend.name,"cpu");
verifyEqual(testCase,backend.reason,"system below GPU crossover threshold");
end

function testExplicitGpuSolveWhenAvailable(testCase)
assumeGreaterThan(testCase,gpuDeviceCount("available"),0);
matrix = gallery('poisson',8);
rhs = ones(size(matrix,1),1);
[solution,information] = tdgl.compute.solveLinear( ...
    matrix,rhs,tdgl.compute.execution("gpu",'AllowFallback',false));
verifyEqual(testCase,information.name,"gpu");
verifyLessThan(testCase,norm(matrix*solution-rhs)/norm(rhs),2e-12);
end
