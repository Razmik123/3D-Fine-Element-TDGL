function figureHandle = plotCylinderSnapshot(analysis)
%PLOTCYLINDERSNAPSHOT Plot magnetic distribution, condensate, and vortices.

figureHandle = figure('Color','white','Name',sprintf( ...
    'Cylinder TDGL snapshot t=%g',analysis.time));
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact');

nexttile(layout,1);
sample = analysis.sampleCellMask;
scatter3(analysis.cellCenters(:,1),analysis.cellCenters(:,2), ...
    analysis.cellCenters(:,3),14,analysis.magneticFluxDensity(:,3),'filled');
axis equal tight;
xlabel('x/\xi'); ylabel('y/\xi'); zlabel('z/\xi');
title(sprintf('B_z, sample mean %.4g',analysis.metrics.meanSampleBz));
colorbar;
hold on;
scatter3(analysis.cellCenters(sample,1),analysis.cellCenters(sample,2), ...
    analysis.cellCenters(sample,3),5,'k');
hold off;

nexttile(layout,2);
sampleNodes = unique(double(analysis.mesh.cells(sample,:)));
scatter3(analysis.mesh.nodes(sampleNodes,1),analysis.mesh.nodes(sampleNodes,2), ...
    analysis.mesh.nodes(sampleNodes,3),18, ...
    abs(analysis.state.orderParameter(sampleNodes)),'filled');
axis equal tight;
xlabel('x/\xi'); ylabel('y/\xi'); zlabel('z/\xi');
title(sprintf('|\psi| and vortex segments (%d)', ...
    analysis.metrics.segmentCount));
colorbar;
hold on;
segments = analysis.vortexSegments;
for index = 1:size(segments.startPoint,1)
    points = [segments.startPoint(index,:);segments.endPoint(index,:)];
    plot3(points(:,1),points(:,2),points(:,3),'k-','LineWidth',2);
end
hold off;
end
