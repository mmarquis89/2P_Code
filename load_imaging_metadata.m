function [outputMetadata, dataFileObj] = load_imaging_metadata()
%=======================================================================================================================================================
%
%  Prompts user for input data file(s) containing 2P imaging data and loads everything except the actual imaging data array. Also does some basic 
%  pre-processing/metadata extraction before returning it all as a single data structure.
%
%  OUTPUT:
%       outputData = a structure with the following fields:
%               behaviorLabels   = the labels corresponding to each number in the behavior annotation data
%               expDate          = date of experiment in YYYY_MM_DD format
%               goodTrials       = logical vector specifying all the trials that are not missing any behavior video frames
%               laserPower       = power of the 2P laser during imaging
%               nFrames          = the number of frames/trial in the behavior video/annotation data
%               nPlanes          = the total number of imaging planes in the session
%               nTrials          = the total number of trials in the session
%               nVolumes         = the number of imaging volumes per trial in the session
%               origFileNames    = cell array containing the file names of all the raw .tif files that generated the input data structure
%               refImg           = 1 x nPlanes cell array, with each cell containing the average of one plane across all volumes and trials [y, x]
%               scanFrameRate    = the acquisition frame rate of the 2P data
%               sid              = the session ID of the data that was loaded
%               stimDurs         = 1 x nTrials numeric vector containing the duration of the stimulus in each trial
%               stimOnsetTimes   = 1 x nTrials numeric vector containing the onset time of the stimulus in seconds for each trial
%               stimSepTrials    = a structure with 1 x nTrials logical vectors for each stimType
%               stimTypes        = cell array containing the names of each unique trialType in the session
%               trialAnnotations = cell array with behavior annotation data (see process_anvil_annotations() for more info)
%               trialDuration    = total duration of trial in seconds
%               trialType        = cell array with the stimulus type for each trial that went in the data structure
%               volumeRate       = volume acquisition rate for the imaging data
%
%      dataFileObj = a writeable matfile object referencing the original .mat file that was selected. The actual imaging data can be accessed by 
%                    referencing or loading the variable "wholeSession" from this structure. wholeSession has dimensions [y, x, plane, volume, trial].
%
%========================================================================================================================================================

% Prompt user for imaging data file
[dataFile, sessionDataPath, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'B:\Dropbox (HMS)\2P Data\Imaging Data\');
if dataFile == 0
    disp('Initialization cancelled')
    outputMetadata = []; % Skip loading if user clicked "Cancel"
else
    imgData = [];
    disp(['Loading ' dataFile, '...'])
    dataFileObj = matfile([sessionDataPath, dataFile], 'Writable', true); % Fields are: 'wholeSession','trialType','origFileNames','tE_sec', 'scanimageInfo', 'expDate' (scanimageInfo not present in older exps)
    imgData.trialType = dataFileObj.trialType;
    imgData.origFileNames = dataFileObj.origFileNames;
    imgData.scanimageInfo = dataFileObj.scanimageInfo;
    imgData.expDate = dataFileObj.expDate;
    wholeSessionSize = size(dataFileObj, 'wholeSession'); % --> [y, x, plane, volume, trial]
    
    disp([dataFile, ' loaded'])
    
    % Prompt user for a stimulus computer metadata file
    [metadataFile, metaDataPath, ~] = uigetfile('*.mat', 'Select a metadata file', ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', imgData.expDate]);
    if metadataFile == 0
        errordlg('No metadata file selected');
        return
    end
    
    % Extract metadata
    metadata = load([metaDataPath, metadataFile]);
    imgData.sid = metadata.metaData.sid;
    imgData.trialDuration = metadata.metaData.trialDuration;
     
    % Get stimulus timing info from file names
    fileNames = imgData.origFileNames;
    imgData.stimOnsetTimes = []; imgData.stimDurs = [];
    for iTrial = 1:numel(fileNames)
        currName = strsplit(fileNames{iTrial}, {'-', '_'});
        imgData.stimOnsetTimes(iTrial)= str2double(currName(find(cellfun(@strcmp, currName, repmat({'Onset'}, size(currName)))) + 1 ));
        imgData.stimDurs(iTrial) = str2double(currName(find(cellfun(@strcmp, currName, repmat({'Dur'}, size(currName)))) + 1 ));
    end
    
    % Prompt user for annotation data file
    [annotDataFile, annotDataPath, ~] = uigetfile('*.mat', 'Select a behavioral annotation data file if desired', sessionDataPath);
    if annotDataFile == 0
        disp('No behavioral annotation data selected')
        annotDataPath = metaDataPath;
        annotData.behaviorLabels = [];
        annotData.ballStopLabels = [];
        annotData.trialAnnotations = [];
        if exist(fullfile(annotDataPath, ['sid_', num2str(imgData.sid), '_frameCountLog.mat']), 'file')
            % Check frame counts for behavior video
            parentDir = fullfile('B:\Dropbox (HMS)\2P Data\Behavior Vids', imgData.expDate, '_Movies')
            [annotData.goodTrials, ~, ~, ~] = frame_count_check(parentDir, imgData.sid, 25, sum(imgData.trialDuration));
        else
            annotData.goodTrials = [];
        end
    else
        disp(['Loading ' annotDataFile, '...'])
        annotData = load([annotDataPath, annotDataFile]);
        if ~isfield(annotData, 'ballStopLabels')
            annotData.ballStopLabels = [];
        end
        disp([annotDataFile, ' loaded'])
    end
    
    % Combine imaging data and annotation data into one structure   
    outputMetadata = setstructfields(imgData, annotData); 
    outputMetadata.nTrials = wholeSessionSize(5);
    outputMetadata.nPlanes = wholeSessionSize(3);
    outputMetadata.nVolumes = wholeSessionSize(4);
    outputMetadata.stimTypes = sort(unique(outputMetadata.trialType));
    if ~isempty(annotData.trialAnnotations)
        outputMetadata.nFrames = max(cellfun(@height, outputMetadata.trialAnnotations(annotData.goodTrials)));
    else
        outputMetadata.nFrames = [];
    end
    
    % Extract metadata from scanimage info
    if isfield(outputMetadata, 'scanimageInfo')
        outputMetadata.volumeRate = frameStringKeyLookup(strjoin(outputMetadata.scanimageInfo, '\n'), 'scanimage.SI.hRoiManager.scanVolumeRate');
        outputMetadata.laserPower = frameStringKeyLookup(strjoin(outputMetadata.scanimageInfo, '\n'), 'scanimage.SI.hBeams.powers');
        outputMetadata.scanFrameRate = frameStringKeyLookup(strjoin(outputMetadata.scanimageInfo, '\n'), 'scanimage.SI.hRoiManager.scanFrameRate');
    else
        % For backwards compatibility
        outputMetadata.volumeRate = 6.44; % This is true for all older experiments
        outputMetadata.laserPower = [];
        outputMetadata.scanFrameRate = [];
    end
    
    % For compatibility with early experiments
    if ~isfield(outputMetadata, 'expDate')
        cellDate = inputdlg('Please enter experiment date in YYYY_MM_DD format');
        outputMetadata.expDate = cellDate{1};
    end
        
    % Separate trials by stim type
    outputMetadata.stimSepTrials = [];
    for iStim = 1:length(outputMetadata.stimTypes)
        outputMetadata.stimSepTrials.(outputMetadata.stimTypes{iStim}) = logical(cellfun(@(x) ...
            strcmp(x, outputMetadata.stimTypes{iStim}), outputMetadata.trialType));
    end
    
    % Match frame times to volumes if annotation data was provided (contains the video frame that 
    % most closely matches the time of the volume for each volume in the trial)
    if ~isempty(outputMetadata.trialAnnotations)
        volTimes = (1:outputMetadata.nVolumes)' ./ outputMetadata.volumeRate;
        frameTimes = outputMetadata.trialAnnotations{find(outputMetadata.goodTrials, 1)}.frameTime;
        for iVol = 1:outputMetadata.nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
            volFrames = volFrames';
        end
        outputMetadata.volFrames = volFrames;
    else
        outputMetadata.volFrames = 1:outputMetadata.nVolumes;
    end
    
    % Prompt user for a reference images file
    [refImageFile, refImgPath, ~] = uigetfile('*.mat', 'Select a reference images file if desired', sessionDataPath);
    if refImageFile == 0
        disp('No reference images file selected - creating reference images from main data file')
        
        % Create mean reference image for each plane
        outputMetadata.refImg = [];
        for iPlane = 1:outputMetadata.nPlanes
            outputMetadata.refImg{iPlane} = squeeze(mean(mean(dataFileObj.wholeSession(:,:,iPlane,:,:),4),5)); % --> [y, x]
        end
    else
        refImages = load([refImgPath, refImageFile]);
        outputMetadata.refImg = refImages.refImages;
    end
    
    % Order fields alphabetically
    outputMetadata = orderfields(outputMetadata);
    
end