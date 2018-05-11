function count_vid_frames(parentDir, sid)

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

logFile = fopen(fullfile(parentDir, ['sid_', num2str(sid), '_frameCounts.txt']), 'w');
vidFiles = dir(fullfile(parentDir, ['sid_', num2str(sid), '_tid*.avi']));

for iFile = 1:numel(vidFiles)
    nFrames = count_frames(fullfile(parentDir, vidFiles(iFile).name));
    fprintf(logFile, [num2str(nFrames), ',', vidFiles(iFile).name, '\n']);
end
fclose(logFile);