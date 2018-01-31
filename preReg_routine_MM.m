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
        chanData = zeros(size(rawFile)); chanData = chanData(:,:,1:end-4,:,:);
        nChannels = size(chanData, 5);
 
        % Save scanimage data from the first trial
        if iFile == 1
            tifObj = Tiff(fullfile(parentDir, currFiles{iFile}), 'r');
            infoStr = tifObj.getTag('ImageDescription');  
            scanimageInfo = strsplit(infoStr, '\n')';
        end
              
        for iChannel = 1:nChannels
            
            % Offset so minimum value = 1
            workingFile = squeeze(rawFile(:,:,:,:,iChannel));
            workingFile = workingFile-min(workingFile(:));
            workingFile = workingFile+1;
            
%             % Crop edges (why do I even do this here?? Need to ask AKM what this part was for)
%             yrange = 3:size(workingFile,1)-2;
%             xrange = 3:size(workingFile,2)-2;
%             croppedFile = workingFile(yrange,xrange,:,:);

            % Discard flyback frames
            workingFile(:,:,1:4,:) = [];                         % --> [y, x, plane, volume]            
            chanData(:,:,:,:,iChannel) = workingFile(:,:,:,:,1); % --> [y, x, plane, volume, channel]
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
        elseif ~isempty(strfind(fName, 'StopBall'))
            trType = 'StopBall';
        elseif ~isempty(strfind(fName, 'OdorA'))
            trType = 'OdorA';
        elseif ~isempty(strfind(fName, 'OdorB'))
            trType = 'OdorB';
        elseif ~isempty(strfind(fName, 'NoOdor'))
            trType = 'NoOdor';
        else
            error(['Error: "', trType, '" is not a valid trial type']);
        end
        
        % Create session array(s) on first loop
        if iFile == 1
            disp('Creating session data array(s)...')
            arraySize = size(squeeze(chanData(:,:,:,:,1))); % --> [y, x, plane, volume]
            arraySize(end+1) = size(currFiles,1);           % --> [y, x, plane]
            wholeSession_1 = uint16(zeros(arraySize));      % --> [y, x, plane, volume, trialNum]
            if nChannels > 1
                % Make a second array if there are two channels of imaging data
                wholeSession_2 = uint16(zeros(arraySize));  % --> [y, x, plane, volume, trialNum]
            end
        end
        
        % Save data to session array(s)
        wholeSession_1(:,:,:,:,iFile) = squeeze(chanData(:,:,:,:,1));     % --> [y, x, plane, volume, trial]
        if nChannels > 1
            wholeSession_2(:,:,:,:,iFile) = squeeze(chanData(:,:,:,:,2)); % --> [y, x, plane, volume, trial]
        end
        trialType{iFile} = trType;
        origFileNames{iFile} = fName;
        
        disp(['Session ', num2str(iSession), ', Trial #', num2str(iFile), ' of ', num2str(length(currFiles))]);
        
    end%for
    
    % Save session data
    if nChannels > 1
        
        % Use the red channel to create and save reference images
        channelNum = 2;
        refImages = [];
        for iPlane = 1:size(wholeSession_2, 3)
            refImages{iPlane} = squeeze(mean(mean(wholeSession_2(:,:,iPlane,:,:),4),5)); % --> [y, x]
        end
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_refImages.mat', iSession)), 'refImages', 'channelNum');
            
        % Save each channel as a separate file
        disp('Saving channel 1...')
        wholeSession = wholeSession_1; % --> [y, x , plane, volume, trial])
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_Chan_1_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate', 'scanimageInfo');
        
        disp('Saving channel 2...')
        wholeSession = wholeSession_2; % --> [y, x , plane, volume, trial])
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_Chan_2_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate', 'scanimageInfo');
        
        % Make and save another array for ratiometric analysis if necessary
        disp('Saving ratio channel...')
        wholeSession = wholeSession_1 ./ wholeSession_2; % --> [y, x, plane, volume, trial])
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_ChanRatio_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate', 'scanimageInfo');     
        
    else
        
        % Use the GCaMP channel to create and save reference images
        channelNum = 1;
        refImages = [];
        for iPlane = 1:size(wholeSession_1, 3)
            refImages{iPlane} = squeeze(mean(mean(wholeSession_1(:,:,iPlane,:,:),4),5)); % --> [y, x]
        end
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_refImages.mat', iSession)), 'refImages', 'channelNum');
        
        disp('Saving session data...')
        wholeSession = wholeSession_1; % --> [y, x , plane, volume, trial])
        savefast(fullfile(sessionDir, sprintf('sid_%.0f_sessionFile.mat',iSession)),'wholeSession','trialType','origFileNames', 'expDate', 'scanimageInfo');
    end
end%for

end%function



