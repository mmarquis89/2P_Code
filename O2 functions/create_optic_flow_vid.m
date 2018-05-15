function create_optic_flow_vid(parentDir, inputVid, roiDataFile, frameCountFile, varargin)
%=======================================================================================================
% CREATE A MOVIE WITH BEHAVIOR VID COMBINED WITH OPTIC FLOW DATA
%
% Uses previously defined behavior vid ROI to create a movie that combines the behavior video with a
% sliding plot of optic flow around the fly to make behavioral annotation in ANVIL easier. Requires 
% a .txt file with comma-separated frame count data for all the trials in <frame count>,<sid_X_tid_XXX>
% format.
%
% INPUTS:
%       parentDir   = the directory containing the optic flow video%ARTICLE_TITLE%
%
%       inputVid    = the name of the input video file to use (minus the file extension)
%
%       roiDataFile     = 
% 
%       frameCountFile = 
%
%       FRAME_RATE  = the frame rate of the behavior video
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'FrameRate' = (default: 25) the frame rate that the behavior video was recorded at 
% 
%       'OutputDir' = (default: parentDir) the full path to the directory to save the output video in
%
%       'OpticFlowFile' = (default: []) the full path to a .mat file containing optic flow data
%
%       'OutputFileName' = (default: [inputVid, '_With_Optic_Flow']) name for the output video file
%
%========================================================================================================

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'FrameRate', 25);
addParameter(p, 'OutputDir', parentDir);
addParameter(p, 'OpticFlowFile', []);
addParameter(p, 'OutputFileName', [inputVid, '_With_Optic_Flow']);
parse(p, varargin{:});
frameRate = p.Results.FrameRate;
outputDir = p.Results.OutputDir;
flowDataFile = p.Results.OpticFlowFile;
outputFileName = p.Results.OutputFileName;


logFile = fopen('logFile.txt', 'w');


% Load ROI data
load(roiDataFile);

% Load frame count log
frameCountFile = fopen(fullfile(parentDir, frameCountFile), 'r');
frameCounts = [];
vidCount = 0;
currLine = fgetl(frameCountFile);
while ischar(currLine)
    tid = str2double(regexp(currLine, '(?<=tid_)...', 'match'));
    sid = str2double(regexp(currLine, '(?<=sid_).', 'match'));
    frameCounts(tid + 1) = str2double(regexp(currLine, '.*(?=,)', 'match'));
    currLine = fgetl(frameCountFile);
    vidCount = vidCount + 1;
end
framesPerTrial = mode(frameCounts);
fclose(frameCountFile);

% Calculate optic flow for each movie frame
vidFile = fullfile(parentDir, [inputVid, '.avi']);
myVid = VideoReader(vidFile);
frameCount = 0;
opticFlow = opticalFlowFarneback;
while hasFrame(myVid)
    
    frameCount = frameCount + 1;
    currFrame = readFrame(myVid);
    currFrame = currFrame(:,:,1);
    
    % Calculate optic flow within each ROI
    currFrameFlowData = estimateFlow(opticFlow, currFrame);
    currFlowFly = currFrameFlowData.Magnitude;
    currFlowFly(~roiData) = nan;
    meanFlowMag(frameCount) = nanmean(nanmean(currFlowFly));
    
end
nFrames = frameCount;

% Calculate frame times
frameTimes = (1:1:nFrames) ./ frameRate;

% Normalize flow magnitudes
flyFlow = meanFlowMag;
flyFlowNorm = flyFlow ./ max(flyFlow(2:end)); % First frame is artificially high so don't use that

% Save optic flow data
save(fullfile(outputDir, ['sid_', num2str(sid), '_optic_flow_data.mat']), 'flyFlowNorm', 'nFrames', 'frameTimes');

fprintf(logFile, 'optic flow calculation complete\n');

% Recreate vidReader
myVid = VideoReader(vidFile);

% Create vidWriter
myVidWriter = VideoWriter(fullfile(outputDir, outputFileName), 'Motion JPEG AVI');
myVidWriter.FrameRate = frameRate;
open(myVidWriter)

% Plot and write each frame
for iFrame = 1:nFrames
    
    disp(num2str(iFrame))
    if ~mod(iFrame, framesPerTrial)
        fprintf(logFile, num2str(iFrame/framesPerTrial));
        disp(['Plotted trial #' num2str(iFrame/framesPerTrial)])
    end
    currFrame = uint8(readFrame(myVid));
    ySize = size(currFrame, 1);
    xSize = size(currFrame, 2);
    
    % Create figure
    screenSize = [1 1 1824 1026];
    h = figure(10);clf
    h.OuterPosition = [50 50 (xSize) screenSize(4)-50];
    
    % Movie frame plot
    ax = axes('Units', 'Normalized', 'Position', [0 0 1 0.7]);
    imshow(currFrame, []);
    axis normal;
    ax.Units = 'Pixels';
    minPos = min(ax.Position(3:4));
    ax.Position(3:4) = [minPos, minPos];
    h.Position(3) = ax.Position(3);
    axis off
    set(gca, 'xticklabel', []);    
    
    % Calculate optic flow xLims
    currFrameTime = iFrame * (1/frameRate);
    preFrameTime = 4;
    postFrameTime = 16;
    xL = [currFrameTime - preFrameTime, currFrameTime + postFrameTime];
    if xL(1) <= 0
        xL(1) = 0;
        xL(2) = preFrameTime + postFrameTime;
    elseif xL(2) > frameTimes(end)
        xL(1) = frameTimes(end) - (preFrameTime + postFrameTime);
        xL(2) =  frameTimes(end);
    end
    
    % Optic flow plot
    trialBoundTimes = [];
    runningCount = 1;
    for iTrial = 1:(length(frameCounts))
        if runningCount > length(frameTimes) % This prevents an error if the first one or more trials has zero frames
           runningCount = length(frameTimes); 
        end
        trialBoundTimes(iTrial) = frameTimes(runningCount);
        runningCount = runningCount + frameCounts(iTrial);
    end
    trialBoundTimes(end+1) = frameTimes(end);
    
    axes('Units', 'Pixels', 'Position', [0 ax.Position(4) ax.Position(3) (h.Position(4) - ax.Position(4))]);
    hold on
    plot(frameTimes(2:end), smooth(flyFlowNorm(2:end), 5), 'color', 'm');    % Plot fly movmement ROI flow
    plot([currFrameTime, currFrameTime], ylim(), 'LineWidth', 2, 'color', 'r')
    for iTrial = 1:length(trialBoundTimes)
       plot([trialBoundTimes(iTrial), trialBoundTimes(iTrial)], ylim(), 'LineWidth', 2, 'color', 'k') 
    end
    set(gca, 'xticklabel', []);
    xlim(xL);
    ylabel('Optic flow (au)');
    lgd = legend('Fly movmement');
    lgd.LineWidth = 1;
    
    drawnow()
    
    % Write frame to video
    writeFrame = getframe(h);
    disp(size(writeFrame.cdata));
    if iFrame > 1
        % Something's messed up with the first frame's size...
        writeVideo(myVidWriter, writeFrame);   
    end
    if iFrame == 2
        % ...so I'm replacing it with another copy of the second frame
        writeVideo(myVidWriter, writeFrame);
    end
end

close(myVidWriter)

end