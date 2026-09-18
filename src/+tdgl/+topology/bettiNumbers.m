function betti = bettiNumbers(mesh, maximumDenseSize)
%BETTINUMBERS Compute mesh Betti numbers from chain-complex ranks.
%   This diagnostic uses numerical dense rank and is intended for geometry
%   validation on modest meshes. A scalable cohomology-basis algorithm will
%   replace it before large multiply-connected production runs.

arguments
    mesh struct
    maximumDenseSize (1,1) double {mustBePositive} = 5000
end

[G, C, D] = tdgl.topology.incidenceMatrices(mesh);
if max([numel(G), numel(C), numel(D)]) > maximumDenseSize^2
    error('tdgl:topology:RankProblemTooLarge', ...
        ['The diagnostic dense-rank calculation is too large. ', ...
         'Use a smaller topology-only mesh.']);
end

rankG = rank(full(G));
rankC = rank(full(C));
rankD = rank(full(D));
nV = size(G,2);
nE = size(G,1);
nF = size(C,1);
nT = size(D,1);

betti = [nV-rankG, nE-rankG-rankC, nF-rankC-rankD, nT-rankD];
end
