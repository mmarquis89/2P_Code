function preReg_routine_MM(parentDir, sids, expDate)
%===================================================================================================
%% Load ScanImage 5.1 imaging data
% 
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
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
    sessionDir = [parentDir sprintf('\\sid_%.0f',iSession)];
    if ~isdir(sessionDir)
        mkdir(sessionDir)
    end
    
    % Separate names of files from the current session
    currFiles = fileNames(sessionNums == iSession);
    
    for iFile = 1:length(currFiles)
        
        % Load a .tif file
        rawFile = read_tif(fullfile(parentDir, currFiles{iFile}));
        rawFile = squeeze(rawFile); % [Lines Pixels Planes Volumes Channels]
        chanData = zeros(size(rawFile)); chanData = chanData(1:end-4,1:end-4,1:end-4,:,:);
        nChannels = size(chanData, 5);
        
        for iChannel = 1:nChannels
            % Offset so minimum value = 1
            workingFile = squeeze(rawFile(:,:,:,:,iChannel));
            workingFile = workingFile-min(workingFile(:));
            workingFile = workingFile+1;
%             % Crop edges (why do I even do this here?? Need to ask AKM what this is for)
%             yrange = 3:size(workingFile,1)-2;
%             xrange = 3:size(workingFile,2)-2;
%             croppedFile = workingFile(yrange,xrange,:,:);

            % Discard flyback frames
            workingFile(:,:,1:4,:) = [];  % --> [x, y, plane, volume]
            
            % Median filter each plane
%             filtFile = uint16(zeros(size(croppedFile)));
%             for iVol = 1:size(croppedFile,4)
%                 for iPlane = 1:size(croppedFile,3)
%                     working_frame = medfilt2(croppedFile(:,:,iPlane,iVol),[1 1]);
%                     filtFile(:,:,iPlane,iVol) = working_frame;
%                 end
%             end
            
            chanData(:,:,:,:,iChannel) = workingFile(:,:,:,:,1); % --> [x, y, plane, volume, channel]
        end
        
        % Separate channels if file includes data from both PMTs
        if nChannels == 1
            chanData_1 = squeeze(chanData);
        else
            chanData_1 = squeeze(chanData(:,:,:,:,1));
            chanData_2 = squeeze(chanData(:,:,:,:,2));
        end
        
        %----- Save to session structure -----
        fName = currFiles{iFile};

        % Get trial type (LeftWind, RightWind, CenterWind, OdorNoWind)
        if ~isempty(strfind(fName, 'LeftWind'))
            trType = 'LeftWind';
        elseif ~isempty(strfind(fName, 'RightWind'))
            trType = 'RightWind';
        elseif ~isempty(strfind(fName, 'CenterWind'))
            trType = 'CenterWind';
        elseif ~isempty(strfind(fName, 'OdorNoWind'))
            trType = 'OdorNoWind';
        elseif ~isempty(strfind(fName, 'NoWind'))
            trType = 'NoWind';
        else
            error(['Error: "', trType, '" is not a valid trial type']);
        end
        
        % Create session array(s) on first loop
        if iFile == 1
            disp('Creating session data array(s)...')
            arraySize = size(squeeze(chanData(:,:,:,:,1))); % --> [x, y, plane, volume]
            arraySize(end+1) = size(currFiles,1);           % --> [x, y, plane
            wholeSession_1 = uint16(zeros(arraySize));      % --> [x, y, plane, volume, trialNum]
            if nChannels > 1
                % Make a second array if there are two channels of imaging data
                wholeSession_2 = uint16(zeros(arraySize));  % --> [x, y, plane, volume, trialNum]
            end
        end
        
        % Save data to session array(s)
        wholeSession_1(:,:,:,:,iFile) = squeeze(chanData(:,:,:,:,1));     % --> [x, y, plane, volume, trial]
        if nChannels > 1
            wholeSession_2(:,:,:,:,iFile) = squeeze(chanData(:,:,:,:,2)); % --> [x, y, plane, volume, trial]
        end
        trialType{iFile} = trType;
        origFileNames{iFile} = fName;
        
        disp(['Session ', num2str(iSession), ', Trial #', num2str(iFile), ' of ', num2str(length(currFiles))]);
        
    end%for
    
    % Save session data
    if nChannels > 1
        
        % Save each channel as a separate file
        disp('Saving channel 1...')
        wholeSession = wholeSession_1; % --> [x, y , plane, volume, trial])
        save(fullfile(sessionDir, sprintf('sid_%.0f_Chan_1_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate');
        
        disp('Saving channel 2...')
        wholeSession = wholeSession_2; % --> [x, y , plane, volume, trial])
        save(fullfile(sessionDir, sprintf('sid_%.0f_Chan_2_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate');
        
        % Make and save another array for ratiometric analysis if necessary
        disp('Saving ratio channel...')
        wholeSession = wholeSession_1 ./ wholeSession_2; % --> [x, y, plane, volume, trial])
        save(fullfile(sessionDir, sprintf('sid_%.0f_ChanRatio_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate');     
        
    else
        disp('Saving session data...')
        wholeSession = wholeSession_1; % --> [x, y , plane, volume, trial])
        save(fullfile(sessionDir, sprintf('sid_%.0f_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate');
    end
end%for

end%function



