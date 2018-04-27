function [outputMetadata, dataFileObj] = load_imaging_metadata(varargin)
%=======================================================================================================================================================
%
%  Prompts user for input data file(s) containing 2P imaging data and loads everything except the actual imaging data array. Also does some basic 
%  pre-processing/metadata extraction before returning it all as a single data structure.
%
%
%  OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%  These are collectively the default assumtions for filenames, which should all be in the same directory as the 
%  selected session data file. Pass arguments to override, or [] to prompt a file selection dialog for the 
%  indicated filename.
%
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
%      dataFileObj = a writeable matfile object referencing the original .mat file that was selected. The actual imaging data can be accessed by 
%                    referencing or loading the variable "wholeSession" from this structure. wholeSession has dimensions [y, x, plane, volume, trial].
%
%========================================================================================================================================================


% Parse optional arguments
p = inputParser;
addParameter(p, 'AnnotFile', 'Annotations.mat');
addParameter(p, 'imgMetadataFile', 'imgMetadata.mat');
addParameter(p, 'refImgFile', 'refImages_Reg.mat')
parse(p, varargin{:});
annotFileName = p.Results.AnnotFile;
imgMetadataFileName = p.Results.imgMetadataFile;
refImgFileName = p.Results.refImgFile;

% Prompt user for imaging data file
[dataFile, sessionDataPath, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'B:\Dropbox (HMS)\2P Data\Imaging Data\');
if dataFile == 0
    % Skip loading if user clicked "Cancel"
    disp('Initialization cancelled')
    outputMetadata = []; 
    dataFileObj = [];
else
    imgMetadata = [];
    disp(['Loading ' dataFile, '...'])
    dataFileObj = matfile([sessionDataPath, dataFile], 'Writable', true); % Only field is 'wholeSession'
    if isempty(imgMetadataFileName)
       [imgDataFile, imgDataFilePath, ~] = uigetfile('*.mat', 'Select an imaging metadata file', sessionDataPath);
       if imgDataFile == 0
           errordlg('No imaging metadata file selected!');
           return
       end
       load(fullfile(imgDataFilePath, imgDataFile)); % variable "imgMetadata" with fields 'trialType','origFileNames','tE_sec', 'scanimageInfo', 'expDate' (scanimageInfo not present in older exps)
    else
       load(fullfile(sessionDataPath, imgMetadataFileName)); % variable "imgMetadata" with fields 'trialType','origFileNames','tE_sec', 'scanimageInfo', 'expDate' (scanimageInfo not present in older exps) 
    end
    wholeSessionSize = size(dataFileObj, 'wholeSession'); % --> [y, x, plane, volume, trial]
    disp([dataFile, ' loaded'])
    
    % Prompt user for a stimulus computer metadata file
    [metadataFile, metaDataPath, ~] = uigetfile('*.mat', 'Select a metadata file', ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', imgMetadata.expDate]);
    if metadataFile == 0
        errordlg('No metadata file selected');
        return
    end
    
    % Extract metadata
    metadata = load([metaDataPath, metadataFile]);
    imgMetadata.sid = metadata.metaData.sid;
    imgMetadata.trialDuration = metadata.metaData.trialDuration;
     
    % Get stimulus timing info from file names
    fileNames = imgMetadata.origFileNames;
    imgMetadata.stimOnsetTimes = []; imgMetadata.stimDurs = [];
    for iTrial = 1:numel(fileNames)
        currName = strsplit(fileNames{iTrial}, {'-', '_'});
        imgMetadata.stimOnsetTimes(iTrial)= str2double(currName(find(cellfun(@strcmp, currName, repmat({'Onset'}, size(currName)))) + 1 ));
        imgMetadata.stimDurs(iTrial) = str2double(currName(find(cellfun(@strcmp, currName, repmat({'Dur'}, size(currName)))) + 1 ));
    end
    
    % Prompt user for annotation data file
    if isempty(annotFileName)
        [annotFileName, annotDataPath, ~] = uigetfile('*.mat', 'Select a behavioral annotation data file if desired', sessionDataPath);
    else
        annotDataPath = sessionDataPath;
        if ~exist(fullfile(sessionDataPath, annotFileName), 'file')
            annotFileName = 0;
        end
    end
    if annotFileName == 0
        disp('No behavioral annotation data loaded')
        annotDataPath = metaDataPath;
        annotData.behaviorLabels = [];
        annotData.ballStopLabels = [];
        annotData.trialAnnotations = [];
        if exist(fullfile(annotDataPath, ['sid_', num2str(imgMetadata.sid), '_frameCountLog.mat']), 'file')
            % Check frame counts for behavior video
            parentDir = fullfile('B:\Dropbox (HMS)\2P Data\Behavior Vids', imgMetadata.expDate, '_Movies')
            [annotData.goodTrials, ~, ~, ~] = frame_count_check(parentDir, imgMetadata.sid, 25, sum(imgMetadata.trialDuration));
        else
            annotData.goodTrials = [];
        end
    else
        disp(['Loading ' annotFileName, '...'])
        annotData = load(fullfile(annotDataPath, annotFileName));
        if ~isfield(annotData, 'ballStopLabels')
            annotData.ballStopLabels = [];
        end
        disp([annotFileName, ' loaded'])
    end

    % Combine imaging data and annotation data into one structure   
    outputMetadata = setstructfields(imgMetadata, annotData); 
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
    if isempty(refImgFileName)
        [refImgFileName, refImgPath, ~] = uigetfile('*.mat', 'Select a reference images file if desired', sessionDataPath);
    else
        refImgPath = sessionDataPath;
        if ~exist(fullfile(sessionDataPath, refImgFileName), 'file')
            refImgFileName = 0;
        end
    end
    if refImgFileName == 0
        disp('No reference images file selected - creating reference images from main data file')
        
        % Create mean reference image for each plane
        outputMetadata.refImg = [];
        for iPlane = 1:outputMetadata.nPlanes
            sessionData = dataFileObj.wholeSession(:,:,iPlane,:,:);
            outputMetadata.refImg{iPlane} = squeeze(mean(mean(sessionData,4),5)); % --> [y, x]
        end
        clear sessionData
    else
        refImages = load(fullfile(refImgPath, refImgFileName));
        outputMetadata.refImg = refImages.refImages;
    end
    
    % Order fields alphabetically
    outputMetadata = orderfields(outputMetadata);
    
end