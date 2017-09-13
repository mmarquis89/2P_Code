function output = system7zip(parentDir, archiveName, archiveType, filterString)
%==================================================================================================================
% ADD FILES TO A COMPRESSED ARCHIVE
% Calls 7zip from the Windows command line to archive files in one of several formats.
%
% Inputs:
%       parentDir = the directory containing the files you want to archive (e.g. 'D:\2p Behavior Video\2017_07_30')
%       archiveName = the name (minus file extension) of the archive file (e.g. 'RawFrames')
%       archiveType = string containing the mane of any 7zip-supported archive format (e.g. 'zip')
%                     The valid values for archiveType are: 'zip', '7z', 'tar', and 'wim'
%       filterString = string to specify the files in parentDir to add to the archive (e.g. *sid_0_tid*_mp4)
%==================================================================================================================

% Make sure archiveType is valid
validTypes = {'zip', '7z', 'tar', 'wim'};
assert(sum(strcmp(archiveType, validTypes)) == 1, ['Error: "', archiveType, ...
        '" is not a ...supported archive type. The valid types are: "zip", "7z", "tar", and "wim"']);

% Make sure a file with this name does not already exist
assert(exist(fullfile(parentDir,[archiveName, '.', archiveType]), 'file')== 0, 'Error: an archive with that name already exists in this directory');

% Add files to archive
disp(['Adding files to ', archiveName, '.', archiveType, '...'])
cmdStr_1 = [parentDir(1:2), ' & ']; % Switch working directory to the correct drive if necessary
cmdStr_2 = ['cd ', parentDir, ' & ']; % Set working directory to parentDir
cmdStr_3 = '"C:\Program Files\7-Zip\7z.exe"'; % Need to provide the entire path for some reason
cmdStr_4 = [' a ', archiveName, '.', archiveType, ' ', filterString]; % Create 7zip command
output = system([cmdStr_1, cmdStr_2, cmdStr_3, cmdStr_4]);
disp('Archival complete')

end