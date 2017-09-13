function archiveRawFrames(sid, parentDir)
%==========================================================================================================
% ADD RAW VIDEO FRAMES TO A COMPRESSED ARCHIVE
% Calls 7zip from the Windows command line to archive raw video file
%
% Inputs:
%       sid = the session ID of the videos you want to process
%       parentDir = the directory containing the raw data for the video frames you want to combine
%               e.g. 'U:\2P Behavior Video\2017_07_30\_Movies'
%==========================================================================================================

    % Make sure this .7z file doesn't already exist
    zipFileName = ['sid_', num2str(sid), '_RawFrames'];
    assert(exist(fullfile(parentDir, [zipFileName, '.7z']), 'file')==0, 'Error: a .7z file with that name already exists in this directory');

    % Zip images
    tic
    disp('Zipping raw video data...');
    
    cmdStr_1 = 'D: & ' % Switch to the data drive
    cmdStr_2 = ['cd ', parentDir, ' & '] % Set working directory to parentDir
    cmdStr_3 = ['7z a sid_', num2str(sid), '_RawFrames.7z
    system(['cd ', parentDir]);
    
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