%RUN_CYLINDER_VORTEX_ENTRY Execute and visualize the production experiment.
repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repositoryRoot);
startup;
configuration = cylinder_vortex_entry_config;
output = tdgl.experiments.runCylinderVortexEntry(configuration);
trajectory = h5info(fullfile(output.runDirectory,'trajectory.h5'),'/time');
lastSnapshot = trajectory.Dataspace.Size(2);
analysis = tdgl.experiments.analyzeCylinderRun( ...
    output.runDirectory,lastSnapshot);
figureHandle = tdgl.visualization.plotCylinderSnapshot(analysis);
exportgraphics(figureHandle,fullfile(output.runDirectory,'final-state.png'), ...
    'Resolution',180);
fprintf('Run saved to %s\n',output.runDirectory);
fprintf('Final vortex segments: %d\n',analysis.metrics.segmentCount);
