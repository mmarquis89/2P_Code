function [outputMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, varargin)
%=======================================================================================================================================================
%
%  Loads a data file containing 2P imaging data. Also does some basic
%  pre-processing/metadata extraction before returning it all as a single data structure.
%
%
%       parentDir                = the path to a directory containing all the data/metadata files
%
%       sessionDataFile          = a session data file containing just an array of imaging data
%                                  named 'wholeSession' with dimensions [y, x, plane, volume, trial]
%
%  OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%  These are collectively the default assumtions for filenames, which should all be in the same directory as the
%  selected session data file. Pass arguments to override, or [] to prompt a file selection dialog for the
%  indicated filename.
%
%       'StimMetadataFile'       = (default: 'stimMetadata.mat')
%       'AnnotFile'              = (default: 'Annotations.mat')
%       'imgMetadataFile'        = (default: 'imgMetadata.mat')
%       'refImgFile'             = (default: 'refImages_Reg.mat')
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
%      wholeSession =  the imaging data array with dimensions [y, x, plane, volume, trial].
%
%========================================================================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'StimMetadataFile', 'stimMetadata.mat');
addParameter(p, 'AnnotFile', 'Annotations.mat');
addParameter(p, 'imgMetadataFile', 'imgMetadata.mat');
addParameter(p, 'refImgFile', 'refImages_Reg.mat')
parse(p, varargin{:});
stimMetadataFileName = p.Results.StimMetadataFile;
annotFileName = p.Results.AnnotFile;
imgMetadataFileName = p.Results.imgMetadataFile;
refImgFileName = p.Results.refImgFile;

% Load imaging data session file
disp(['Loading ' sessionDataFile, '...'])
load(fullfile(parentDir, sessionDataFile)) % only variable is 'wholeSession'
wholeSessionSize = size(wholeSession);
disp('Session data loaded')

% Load imaging metadata file
imgMetadata = load(fullfile(parentDir, imgMetadataFileName)); % fields 'trialType','origFileNames','tE_sec', 'scanimageInfo', 'expDate' (scanimageInfo not present in older exps)   

% Load stimulus computer metadata file
load(fullfile(parentDir, stimMetadataFileName)); % variable "metaData" with fields 'trialDuration, 'interTrialInterval', 'stimType', 'sid', 'taskFile', 'outputData'

% Combine metadata
imgMetadata.sid = metaData.sid;
imgMetadata.trialDuration = metaData.trialDuration;

% Get stimulus timing info from file names
fileNames = imgMetadata.origFileNames;
imgMetadata.stimOnsetTimes = []; imgMetadata.stimDurs = [];
for iTrial = 1:numel(fileNames)
    currName = fileNames{iTrial};
    imgMetadata.stimOnsetTimes(iTrial) = str2double(regexp(currName, '(?<=Onset-).*(?=-Dur)', 'match'));
    imgMetadata.stimDurs(iTrial) = str2double(regexp(currName, '(?<=Dur-).*(?=_tt)', 'match'));
end

% Load annotation data file
if ~isempty(annotFileName)
    disp(['Loading ' annotFileName, '...'])
    annotData = load(fullfile(parentDir, annotFileName)); % variables 'ballStopLabels, 'behaviorLabels', 'goodTrials', 'trialAnnotations'
    if ~isfield(annotData, 'ballStopLabels')
        annotData.ballStopLabels = [];
    end
    disp([annotFileName, ' loaded'])
else
    disp('No behavioral annotation data loaded')
    annotData.behaviorLabels = [];
    annotData.ballStopLabels = [];
    annotData.trialAnnotations = [];
    annotData.goodTrials = [];
end

% Combine imaging data and annotation data into one structure
outputMetadata = setstructfields(imgMetadata, annotData);
outputMetadata.nTrials = size(wholeSession, 5);
outputMetadata.nPlanes = size(wholeSession, 3);
outputMetadata.nVolumes = size(wholeSession, 4);
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

% Load a reference images file or create new reference images
if ~isempty(refImgFileName)
    refImages = load(fullfile(parentDir, refImgFileName));
    outputMetadata.refImg = refImages.refImages;
else
    disp('No reference images file selected - creating reference images from main data file')
    
    % Create mean reference image for each plane
    outputMetadata.refImg = [];
    for iPlane = 1:outputMetadata.nPlanes
        disp(num2str(iPlane))
        outputMetadata.refImg{iPlane} = squeeze(mean(mean(wholeSession(:,:,iPlane,:,:), 4, 5))); % --> [y, x]
    end
end

% Order fields alphabetically
outputMetadata = orderfields(outputMetadata);
