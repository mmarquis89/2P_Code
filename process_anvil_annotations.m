function trialAnnotations = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, frameRate, trialDuration)
%========================================================================================================================= 
% READS AND PROCESSES A FRAME-BY-FRAME BEHAVIORAL ANNOTATION FILE CREATED IN ANVIL
%
% This function requires information about the frame counts for both the individual and concatenated video files for 
% the selected session, saved in the parent directory as 'sid_0_frameCountLog.mat' and 
% 'sid_0_AllTrials_frameCountLog.mat', respectively (with the appropriate session ID #).
%
% Processed annotation data will be saved as a .mat file in the specified directory containing a table with the 
% annotation data as well as a cell array with labels for the behavioral codes and a logical vector of valid trials.
%
% INPUTS:
%       sid  =  session ID of the data you want to process.
%       parentDir  =  path to the directory containing the annotation file and the frame count logs.
%       saveDir  =  directory to save processed annotation data in.
%       annotationFileName = name of the Anvil-exported .txt file with annotations.
%       frameRate = frame rate that the behavioral video was acquired at.
%       trial Duration = duration of each trial in seconds (must be the same for all trials).
%
% OUTPUS:
%       trialAnnotations = a 1xn cell array (where n is the number of trials in the session), with each cell containing an 
%                          mx3 table (where m is the number of video frames for that trial) with the following column 
%                          names: [frameNum, actionNums, frameTime]. [actions] contains a behavioral code for each frame, 
%                          corresponding to an entry in the "behaviorLabels" array. [frameTime] is the trial time in 
%                          seconds corresponding to each frame.
%
%       behaviorLabels   = an array of strings containing the various strings that correspond to numbers in the [actionNums]
%                          field of the trialAnnotations tables (0 = first entry in behaviorLabels, etc.) Note hardcoded
%                          values for this variable below.
%
%       goodTrials =       a 1 x n logical array (where n is the number of trials in the session) indicating which trials 
%                          are missing one or more video frames.
%                           
%==========================================================================================================================

% Remember to update this if annotation coding changes
behaviorLabels = {'None', 'AbdominalContraction', 'Locomotion', 'Grooming', 'IsolatedLegMovement'};

% Read and parse annotation table data
annotationData = readtable(fullfile(parentDir, annotationFileName));
frameNum = annotationData.Frame;
frameNum = frameNum + 1; % Use 1-indexing for frame numbers
actionNums = annotationData.ActionTypes_ActionTypes;
actionNums(actionNums == -1000) = 0;
annotationTable = table(frameNum, actionNums);

% Import video frame count data 
individualVidFrameCounts = load(fullfile(parentDir, ['sid_', num2str(sid), '_frameCountLog.mat']));
nTrials = length(individualVidFrameCounts.frameCounts);
frameCounts = [individualVidFrameCounts.frameCounts.nFrames];
concatenatedFrameCounts = load(fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials_frameCountLog.mat']));
allTrialsFrameCount = concatenatedFrameCounts.frameCount;

% Check for trials with expected number of frames
goodTrials = (frameCounts == frameRate * trialDuration);

% Calculate frame times in seconds for good trials
frameTimes = (1:frameRate*trialDuration) * (1/frameRate);

% Make sure frame counts are consistent with each other and annotation data
assert(sum(frameCounts) == allTrialsFrameCount, 'Error: sum of individual frame counts is not equal to concatenated video frame count');
assert(allTrialsFrameCount == length(frameNum), 'Error: frame number mismatch between annotation data and video file data');

% Separate annotation data into individual trials
trialAnnotations = cell(1,nTrials);
currFrame = 1;
for iTrial = 1:nTrials
    lastFrame = currFrame+frameCounts(iTrial)-1;
    
    trialAnnotations{iTrial} = annotationTable(currFrame:lastFrame,:);
    
    % Add column of times to table
    if goodTrials(iTrial)
        frameTime = frameTimes';
        trialAnnotations{iTrial} = [trialAnnotations{iTrial}, table(frameTime)];
    else 
        % Replace times with zeros for trials without enough frames
        frameTime = zeros(1, frameCounts(iTrial))';
        trialAnnotations{iTrial} = [trialAnnotations{iTrial}, table(frameTime)];
    end
    currFrame = currFrame + frameCounts(iTrial);    
end

% Make sure output data file doesn't already exist for this session
saveFilePath = fullfile(saveDir, ['sid_', num2str(sid), '_BehavioralAnnotations.mat']);
assert(exist(saveFilePath, 'file') == 0, 'Error: a file with this name already exists in the save directory')

% Save processed annotation data along with behavior labels
save(saveFilePath, 'trialAnnotations', 'behaviorLabels', 'goodTrials');

end%function