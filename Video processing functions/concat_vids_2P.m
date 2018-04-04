function concat_vids_2P(sid, vidDir, frameRate)
%============================================================================================================================
% CONCATENATE ALL COMBINED PLOTTING VIDEOS FOR THE EXPERIMENT
% Concatenates the behavioral videos for each trial of a 2P experiment, and returns a message string indicating whether
% the operation was a success (and if not, which trial it failed on). The new video will be saved in the same folder as the
% source videos.
%
% Inputs:
%       sid = the session ID of the videos you want to process
%       vidDir = the directory containing the videos you want to combine
%               e.g. 'U:\2P Behavior Video\2017_07_30\_Movies'
%============================================================================================================================

FRAME_RATE = frameRate;
vidFiles = dir(fullfile(vidDir, ['sid_', num2str(sid), '*tid*.mp4']));
vidNames = sort({vidFiles.name});

% Figure out the actual trial numberes from filenames
numLocs = cellfun(@(x) (regexp(x, 'tid_') + 4), vidNames);
trialNums = cell2mat(cellfun(@(x,y) str2num(x(y:y+2)), vidNames, num2cell(numLocs), 'UniformOutput', 0));

% Create save directory if it doesn't already exist
if ~isdir(vidDir);
   mkdir(vidDir) 
end

% Make sure concatenated video file doesn't already exist
assert(exist(fullfile(vidDir, ['sid_', num2str(sid),'_AllTrials.mp4']), 'file')==0, 'Error: concatenated video file already exists in this directory');

% Create vidWriter
myVidWriter = VideoWriter(fullfile(vidDir, ['sid_', num2str(sid),'_AllTrials.mp4']), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)
frameCount = 0;

disp('Concatenating videos...')
% try
vidCount = 0;
for iTrial = trialNums
    
    vidCount = vidCount + 1;
    
    % Pad the trial number if necessary to ensure correct filename sorting
    if iTrial < 10
        padStr = '00';
    elseif iTrial < 100
        padStr = '0';
    else
        padStr = '';
    end
    
    trialStr = ['sid_', num2str(sid), '_tid_', padStr, num2str(iTrial)];
    disp(trialStr)
    
    if ~isempty(dir(fullfile(vidDir, [trialStr, '.mp4']))) % Check to make sure is some video for this trial
        
        % Load movie for the current trial
        myMovie = {};
        myVid = VideoReader(fullfile(vidDir, vidNames{vidCount}));
        while hasFrame(myVid)
            currFrame = readFrame(myVid);
            myMovie(end+1) = {uint8(currFrame)};
        end
        
        % Add frames to movie
        for iFrame = 1:length(myMovie)
            writeVideo(myVidWriter, myMovie{iFrame});
            frameCount = frameCount + 1;
        end
    end%if
end%iTrial


close(myVidWriter)
clear('myMovie')

% Save frame count log
save(fullfile(vidDir, ['sid_', num2str(sid) '_AllTrials_frameCountLog.mat']), 'frameCount')

disp(['Videos concatenated successfully! Total frames = ', num2str(frameCount)])
% catch
%     disp(['Error - video making failed on trial #', num2str(iTrial)])
% end%try
end%function