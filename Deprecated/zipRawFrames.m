function zipRawFrames(sid, parentDir)
%============================================================================================================================
% ADD RAW VIDEO FRAMES TO A ZIPPED ARCHIVE
%
% Inputs:
%       sid = the session ID of the videos you want to process
%       parentDir = the directory containing the raw data for the video frames you want to combine
%               e.g. 'U:\2P Behavior Video\2017_07_30\_Movies'
%============================================================================================================================

    % Make sure this zip file doesn't already exist
    zipFileName = ['sid_', num2str(sid), '_RawFrames'];
    assert(exist(fullfile(parentDir, [zipFileName, '.zip']), 'file')==0, 'Error: a zip file with that name already exists in this directory');

    % Zip images
    tic
    disp('Zipping raw video data...');
    sessionDirs = dir(fullfile(parentDir, ['*sid_', num2str(sid), '_tid*']));
    zipPaths = strcat([parentDir, '\'], {sessionDirs.name});
    zip(fullfile(parentDir, zipFileName), zipPaths);
    disp(['Zipping complete. Duration = ', num2str(toc), ' sec'])

%     % Delete raw video frames if zipping was successful
%     disp('Deleting raw video frames...')
%     delDirs = sessionDirs;
%     zipDir = dir(fullfile(parentDir, [zipFileName, '.zip']));
%     assert(~isempty(zipDir), 'Error: no zipped folder was found for these files');
%     for iFile = 1:length(delDirs)
%         rmdir(fullfile(parentDir, delDirs(iFile).name), 's');
%     end
%     disp('Raw video frames deleted')
end