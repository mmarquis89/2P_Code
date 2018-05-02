function preReg_routine_MM(parentDir, sid, expDate)
%===================================================================================================
%% Load ScanImage 5.1 imaging data
%
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
%       sid        = the session ID you want to process.
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

disp('Processing...')

% Make session folder for new files if necessary
sessionDir = [parentDir sprintf('\\sid_%.0f',sid)];
if ~isdir(sessionDir)
    mkdir(sessionDir)
end

% Separate names of files from the current session
currFiles = fileNames(sessionNums == sid);
nTrials = numel(currFiles);

for iFile = 1:nTrials
    
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
        workingFile = workingFile - min(workingFile(:));
        workingFile = workingFile + 1;
        
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
    
    % Create session array(s) on first loop
    if iFile == 1
        chDataSize = size(chanData);
        sz = [chDataSize(1:4), nTrials];
        nVolumes = sz(4);
        if nChannels > 1
            m = matfile(fullfile(sessionDir,sprintf('sid_%.0f_Chan_1_sessionFile.mat',sid)), 'Writable', true);
            m2 = matfile(fullfile(sessionDir,sprintf('sid_%.0f_Chan_2_sessionFile.mat',sid)), 'Writable', true);
            m.wholeSession(sz(1),sz(2),sz(3),(nVolumes * nTrials)) = 0;
            m2.wholeSession(sz(1),sz(2),sz(3),(nVolumes * nTrials)) = 0;
        else
            m = matfile(fullfile(sessionDir,sprintf('sid_%.0f_sessionFile.mat',sid)), 'Writable', true);
            m.wholeSession(sz(1),sz(2),sz(3),(nVolumes * nTrials)) = 0;
        end
    end
    
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
    elseif (~isempty(strfind(fName, 'NoOdor')) || ~isempty(strfind(fName, 'NoStim')))
        trType = 'NoStim';
    elseif ~isempty(strfind(fName, 'CarrierStreamStop'))
        trType = 'CarrierStreamStop';
    elseif ~isempty(strfind(fName, 'Laser'))
        trType = 'Laser';
    else
        error(['Error: "', trType, '" is not a valid trial type']);
    end
    
    % Save data to session array(s)    
    m.wholeSession(1:sz(1), 1:sz(2),1:sz(3), (((iFile-1) * nVolumes) +  1):(iFile * nVolumes)) = chanData(:,:,:,:,1);     % --> [y, x, plane, allVolumes]
    if nChannels > 1
        m2.wholeSession(1:sz(1), 1:sz(2),1:sz(3), (((iFile-1) * nVolumes) +  1):(iFile * nVolumes)) = chanData(:,:,:,:,2); % --> [y, x, plane, allVolumes]
    end
    trialType{iFile} = trType;
    origFileNames{iFile} = fName;
    
    disp(['Session ', num2str(sid), ', Trial #', num2str(iFile), ' of ', num2str(length(currFiles))]);
    
end%iFile

% Save other imaging metadata variables in separate file
savefast(fullfile(parentDir, ['imgMetadata']), 'trialType', 'origFileNames', 'expDate', 'scanimageInfo');

% Create reference images
if nChannels == 2
    channelNum = 2;
    refImages = [];
    for iPlane = 1:sz(3)
        currData = m2.wholeSession(:,:,iPlane,:,:);
        refImages{iPlane} = squeeze(mean(mean(currData,4),5)); % --> [y, x]
    end
    savefast(fullfile(sessionDir, sprintf('sid_%.0f_refImages.mat', iSession)), 'refImages', 'channelNum');
    
    clear m m2
else
    
    % Use the GCaMP channel to create and save reference images
    disp('Creating reference images...')
    channelNum = 1;
    refImages = [];
    tic
    for iPlane = 1:sz(3)
        disp(['Plane #', num2str(iPlane)]);
        currData = m.wholeSession(1:sz(1),1:sz(2),iPlane,1:sz(4),1:sz(5));
        refImages{iPlane} = squeeze(mean(mean(currData,4),5)); % --> [y, x]
    end
    savefast(fullfile(sessionDir, sprintf('sid_%.0f_refImages.mat', sid)), 'refImages', 'channelNum');
    disp(['Reference images saved in ', num2str(toc), ' seconds']);
    clear m
end

disp('Saving complete')

end%function



