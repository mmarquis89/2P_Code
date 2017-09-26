function preRegRoutine_MM(parentDir, sids, expDate)
%===================================================================================================
%% Load ScanImage 5.1 imaging data
% 
% Inputs:
%       parentDir   = directory containing the raw data that you want to process
%       sids        = array containing the session IDs you want to process. Pass [] to default to 
%                     detecting and processing all sids present in the parent directory
%
% NOTE: the files should be sorted in chronological order due to the timestamp at the beginning 
% of the filename. If filenames do not sort in this pattern, they must be renamed before processing.
%===================================================================================================

% Identify data files for each session
myFiles = dir(fullfile(parentDir, '*sid*tid*'));
fileNames = sort({myFiles.name})';
sidLocs = strfind(fileNames, 'sid_');
for iFile = 1:length(fileNames)
    sessionNums(iFile) = str2double(fileNames{iFile}(sidLocs{iFile}+4));
end

% Identify session numbers to be processed
if isempty(sids)
    mySessions = unique(sessionNums);
else
    mySessions = sids;
end

disp('Processing...')
for iSession = mySessions
    
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
    save(fullfile(sessionDir, sprintf('sid_%.0f_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate');
    clear wholeSession trialTypes origFileNames

end%for

end%function



