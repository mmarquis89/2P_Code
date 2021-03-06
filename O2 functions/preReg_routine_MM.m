function preReg_routine_MM(parentDir, sid, expDate, varargin)
%===================================================================================================
% Load ScanImage 5.1 imaging data
%
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
%       sid        = the session ID you want to process.
%       
%       expDate    = the identifier of the experiment you want to process.
% 
% OPTIONAL NAME-VALUE PAIR INPUTS:
%       
%       'OutputDir' = (default: parentDir) the directory to save the processed data in
%
% NOTE: the files should be sorted in chronological order due to the timestamp at the beginning
% of the filename. If filenames do not sort in this pattern, they must be renamed before processing.
%===================================================================================================

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', parentDir);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;

sid = str2double(sid);

% Identify data files for each session
myFiles = dir(fullfile(parentDir, '*sid*tid*'));
fileNames = sort({myFiles.name})';
sidLocs = strfind(fileNames, 'sid_');
for iFile = 1:length(fileNames)
    sessionNums(iFile) = str2double(fileNames{iFile}(sidLocs{iFile}+4));
end

disp('Processing...')

% Make session folder for new files if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Separate names of files from the current session
currFiles = fileNames(sessionNums == sid);
nTrials = numel(currFiles);
disp(nTrials)
for iFile = 1:nTrials
        
    % Load a .tif file
    rawFile = read_tif(fullfile(parentDir, currFiles{iFile}));
    rawFile = squeeze(rawFile); % [Lines Pixels Planes Volumes Channels]
    chanData = zeros(size(rawFile)); chanData = chanData(:,:,1:end-4,:,:);
    nChannels = size(chanData, 5);
    
    % Extract scanimage data from the first trial
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
        chDataSize = size(chanData_1);
        sz = [chDataSize(1:4), nTrials];
        wholeSession = zeros(sz);
        if nChannels > 1
            wholeSession2 = zeros(sz);            
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
    wholeSession(:,:,:,:,iFile) = chanData_1(:,:,:,:,1);
    if nChannels > 1
       wholeSession2(:,:,:,:,iFile) = chanData_2(:,:,:,:,2); 
    end
    trialType{iFile} = trType;
    origFileNames{iFile} = fName;
        
end%iFile


disp(size(wholeSession))
disp('Data processing complete')


% Save session data file
save(fullfile(outputDir, ['sid_', num2str(sid), '_sessionFile.mat']), 'wholeSession', '-v7.3')
if nChannels == 2
   clear wholeSession
   wholeSession = wholeSession2;
   save(fullfile(outputDir, [sid_', num2str(sid), '_Chan_2_sessionFile.mat']), 'wholeSession', '-v7.3');
end

disp(outputDir)
disp('Session data file saved')

% Save other imaging metadata variables in separate file
save(fullfile(outputDir, 'imgMetadata'), 'trialType', 'origFileNames', 'expDate', 'scanimageInfo');

disp('Metadata saved')

% Create reference images
if nChannels == 2

    % Use red channel for reference images
    channelNum = 2;
    refImages = [];
    for iPlane = 1:sz(3)
        refImages{iPlane} = squeeze(mean(mean(wholeSession2(:,:,iPlane,:,:),4),5)); % --> [y, x]
    end
    save(fullfile(outputDir, sprintf('sid_%.0f_refImages.mat', iSession)), 'refImages', 'channelNum', '-v7.3');
else
    
    % Use the GCaMP channel to create and save reference images
    channelNum = 1;
    refImages = [];
    for iPlane = 1:sz(3)
        refImages{iPlane} = squeeze(mean(mean(wholeSession(:,:,iPlane,:,:),4),5)); % --> [y, x]
    end
    save(fullfile(outputDir, sprintf('sid_%.0f_refImages.mat', sid)), 'refImages', 'channelNum', '-v7.3');
end

disp('Reference images created and saved')

end%function