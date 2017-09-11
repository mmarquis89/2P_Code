function preRegRoutine_MM(parentDir)
%% Load ScanImage 5.1 imaging data

% parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05_exp2'; % Update - new data directory that you want to process
% Note that the files should be sorted in chronological order due to the timestamp at the beginning of the filename

% Identify data files for each session
myFiles = dir(fullfile([parentDir, '\*sid*tid*']));
fileNames = sort({myFiles.name})';
sidLocs = strfind(fileNames, 'sid_');
for iFile = 1:length(fileNames)
    sessionNums(iFile) = str2double(fileNames{iFile}(sidLocs{iFile}+4));
end

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

disp('Processing...')
for iSession = unique(sessionNums)
    
    % Make session folder for new files if necessary
    sessionDir=[parentDir sprintf('\\sid_%.0f',iSession)];
    if ~isdir(sessionDir)
        mkdir(sessionDir)
    end
    
    % Separate names of files from the current session
    currFiles = fileNames(sessionNums == iSession);
    
    for iFile = 1:length(currFiles)
        
        % Load a .tif file, subract baseline, and crop edges
        workingFile=readTif(fullfile(parentDir, currFiles{iFile}));
        workingFile=squeeze(workingFile); % [Lines Pixels Planes Volumes]
        workingFile=workingFile-min(workingFile(:));
        workingFile=workingFile+1;
        yrange=2:size(workingFile,1)-2;
        xrange=2:size(workingFile,2)-2;
        croppedFile=workingFile(yrange,xrange,:,:);
        
        % Discard flyback frames
        croppedFile(:,:,1:4,:)=[]; % ScanImage should tell you
        
        % Median filter each plane
        filtFile=uint16(zeros(size(croppedFile)));
        for iVol=1:size(croppedFile,4)
            for iPlane=1:size(croppedFile,3)
                working_frame=medfilt2(croppedFile(:,:,iPlane,iVol),[1 1]);
                filtFile(:,:,iPlane,iVol)=working_frame+1; % Offset 0
            end
        end
        
        
        %----- Save to session structure -----
        fName=currFiles{iFile};
        trialNum=iFile;

        % Get trial type (LeftWind, RightWind, CenterWind, OdorNoWind)
        if ~isempty(strfind(fName, 'LeftWind'))
            trType = 'LeftWind';
        elseif ~isempty(strfind(fName, 'RightWind'))
            trType = 'RightWind';
        elseif ~isempty(strfind(fName, 'CenterWind'))
            trType = 'CenterWind';
        elseif ~isempty(strfind(fName, 'OdorNoWind'))
            trType = 'OdorNoWind';
        end
        
        % Create session file on first loop
        if iFile==1
            arraySize=size(filtFile);
            arraySize(end+1)=size(currFiles,1);
            wholeSession=uint16(zeros(arraySize)); % [x, y ,z, volume, trialNum]
        end
        
        wholeSession(:,:,:,:,trialNum)=filtFile;
        trialType{trialNum}=trType;
        origFileNames{trialNum}=fName;
        
        disp(['Session ', num2str(iSession), ', Trial #', num2str(trialNum)]);
        
    end%for
    
    % Save session data ( [x, y , plane, volume, trialNum] )
    save(fullfile(sessionDir, sprintf('sid_%.0f_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames');
    clear wholeSession trialTypes origFileNames
    
    % Make sure zipped archive doesn't already exist
    if exist(fullfile(parentDir,['TrialData_sid_', num2str(iSession), '.zip']), 'file') == 0
        
        % Zip raw data files
        disp('All session data saved')
        zipFiles = dir(fullfile(parentDir, ['*sid_', num2str(iSession), '_tid*']));
        disp('Zipping individual trial data files...')
        zipPaths = strcat([parentDir, '\'], {zipFiles.name});
        zip(fullfile(parentDir,['TrialData_sid_', num2str(iSession)]), zipPaths);
        disp('Zipping complete')
        
        % Delete raw trial data files after if zipping was successful
        disp('Deleting raw trial data files...')
        delFiles = zipFiles;
        zipDir = dir(fullfile(parentDir, ['TrialData_sid_', num2str(iSession), '.zip']));
        assert(~isempty(zipDir), 'Error: no zipped folder was found for these files');
        for iFile = 1:length(delFiles)
            delete(fullfile(parentDir, delFiles(iFile).name));
        end
        disp('Raw stacks deleted')
    else
        disp('Error: zip archive already exists for session ', num2str(iSession), '...skipping archival and deletion')
    end%if
    
end%for
end%function



