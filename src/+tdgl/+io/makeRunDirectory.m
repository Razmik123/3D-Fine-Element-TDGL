function runDirectory = makeRunDirectory(resultsRoot,experimentName)
%MAKERUNDIRECTORY Return a unique, not-yet-created experiment directory.

arguments
    resultsRoot (1,1) string
    experimentName (1,1) string
end
safeName = regexprep(lower(experimentName),'[^a-z0-9_-]+','-');
safeName = regexprep(safeName,'(^-+|-+$)','');
if strlength(safeName) == 0, safeName = "experiment"; end
stamp = string(datetime('now','Format','yyyyMMdd-HHmmss'));
base = fullfile(resultsRoot,safeName,stamp);
runDirectory = base;
suffix = 1;
while isfolder(runDirectory) || isfile(runDirectory)
    runDirectory = base+"-"+string(suffix);
    suffix = suffix+1;
end
end
