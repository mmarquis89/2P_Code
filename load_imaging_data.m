function outputData = load_imaging_data()
%=======================================================================================================================================================
%
%  Prompts user for input data file(s) containing 2P imaging data and loads the data. Can handle either 1 or 2 channel imaging. Also does some basic
%  pre-processing/metadata extraction before returning it all as a single data structure.
%
%  OUTPUT:
%       outputData = a structure with the following fields:
%               behaviorLabels   = the labels corresponding to each number in the behavior annotation data
%               expDate          = date of experiment in YYYY_MM_DD format
%               goodTrials       = logical vector specifying all the trials that are not missing any behavior video frames
%               nFrames          = the number of frames/trial in the behavior video/annotation data
%               nPlanes          = the total number of imaging planes in the session
%               nTrials          = the total number of trials in the session
%               nVolumes         = the number of imaging volumes per trial in the session
%               origFileNames    = cell array containing the file names of all the raw .tif files that generated the input data structure
%               refImg           = 1 x nPlanes cell array, with each cell containing the average of one plane across all volumes and trials [x, y]
%               sid              = the session ID of the data that was loaded
%               stimSepTrials    = a structure with 1 x nTrials logical vectors for each stimType, as well as one for wind trials (windTrials)
%               stimTypes        = cell array containing the names of each unique trialType in the session
%               tE_sec           = time in seconds it took to register the data
%               trialAnnotations = cell array with behavior annotation data (see process_anvil_annotations() for more info)
%               trialDuration    = [pre-stim, stimLength, postStim] trial time in seconds
%               trialType        = cell array with the stimulus type for each trial that went in the data structure
%               volumeRate       = volume acquisition rate for the imaging data
%               wholeSession     = the imaging data in an array with dimensions [x, y, plane, volume, trial]
%
%========================================================================================================================================================

% Prompt user for imaging data file
[dataFile, pathName, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
if dataFile == 0
    disp('Initialization cancelled')
    outputData = []; % Skip loading if user clicked "Cancel"
else
    disp(['Loading ' dataFile, '...'])
    imgData = load([pathName, dataFile]); % Fields are: 'regProduct','trialType','origFileNames','tE_sec', 'expDate'
    disp([dataFile, ' loaded'])
    
    % Prompt user for a metadata file
    [metadataFile, pathName, ~] = uigetfile('*.mat', 'Select a metadata file', 'D:\Dropbox (HMS)\2P Data\Behavior Vids\');
    if metadataFile == 0
        disp('No metadata file selected')
        imgData.trialDuration = [];
        imgData.volumeRate = [];
    else
        metadata = load([pathName, metadataFile]);
        imgData.trialDuration = metadata.metaData.trialDuration;
        imgData.volumeRate = metadata.metaData.run_obj_frameRate ./ metadata.metaData.run_obj_nVols;
    end
    
    % Prompt user for behavioral annotation data file
    [annotDataFile, pathName, ~] = uigetfile('*.mat', 'Select a behavioral annotation data file if desired', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
    if annotDataFile == 0
        disp('No behavioral annotation data selected')
        annotData.behaviorLabels = [];
        annotData.goodTrials = [];
        annotData.trialAnnotations = [];
    else
        disp(['Loading ' annotDataFile, '...'])
        annotData = load([pathName, annotDataFile]);
        disp([annotDataFile, ' loaded'])
    end
    
    outputData = setstructfields(imgData, annotData); % Combine imaging data and annotation data into one structure   
    
    % Process raw data structure
    if ~isfield(outputData, 'wholeSession')
        outputData.wholeSession = outputData.regProduct; % --> [x, y, plane, volume, trial]
        outputData = rmfield(outputData, 'regProduct');
    end
    outputData.nTrials = size(outputData.wholeSession, 5);
    singleTrial = squeeze(outputData.wholeSession(:,:,:,:,1));  % --> [x, y, plane, volume]
    outputData.nPlanes = size(singleTrial, 3);
    outputData.nVolumes = size(singleTrial, 4);
    outputData.stimTypes = sort(unique(outputData.trialType));
    if ~isempty(annotData.trialAnnotations)
        outputData.nFrames = max(cellfun(@height, outputData.trialAnnotations));
    else
        outputData.nFrames = [];
    end
    
    % For compatibility with early experiments
    if ~isfield(outputData, 'expDate')
        cellDate = inputdlg('Please enter experiment date in YYYY_MM_DD format');
        outputData.expDate = cellDate{1};
    end
    
    % Extract session number
    origFileName = outputData.origFileNames{1};
    sidLoc = strfind(origFileName, 'sid_');
    outputData.sid = origFileName(sidLoc+4);
    
    % Create mean reference image for each plane
    outputData.refImg = [];
    for iPlane = 1:outputData.nPlanes
        outputData.refImg{iPlane} = squeeze(mean(mean(outputData.wholeSession(:,:,iPlane,:,:),4),5)); % --> [x, y]
    end
    
    % Separate trials by stim type
    outputData.stimSepTrials = [];
    for iStim = 1:length(outputData.stimTypes)
        outputData.stimSepTrials.(outputData.stimTypes{iStim}) = logical(cellfun(@(x) ...
            strcmp(x, outputData.stimTypes{iStim}), outputData.trialType));
    end
    
    % Separate out all wind stim trials
    windTrials = zeros(1, outputData.nTrials);
    if isfield(outputData.stimSepTrials, 'LeftWind')
        windTrials = windTrials + logical(outputData.stimSepTrials.LeftWind);
    end
    if isfield(outputData.stimSepTrials, 'RightWind')
        windTrials = windTrials + logical(outputData.stimSepTrials.RightWind);
    end
    if isfield(outputData.stimSepTrials, 'CenterWind')
        windTrials = windTrials + logical(outputData.stimSepTrials.CenterWind);
    end
    outputData.stimSepTrials.windTrials = logical(windTrials);
   
    % Order fields alphabetically
    outputData = orderfields(outputData);
    
end