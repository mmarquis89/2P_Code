function make_optic_flow_vid(sid, parentDir, FRAME_RATE, roiDataFile)
%=======================================================================================================
% CREATE A MOVIE WITH BEHAVIOR VID COMBINED WITH OPTIC FLOW DATA
%
% Uses previously defined behavior vid ROI to create a movie that combines the behavior video with a
% sliding plot of optic flow around the fly to make behavioral annotation in ANVIL easier. This function
% relies on the video being named "sid_X_AllTrials.mp4".
%
% INPUTS:
%       sid = the session ID of the video you want to process
%
%       parentDir = the directory containing the optic flow video
%
%       FRAME_RATE = the frame rate of the behavior video
%
%       vidFile = <OPTIONAL> the path to the .mat file containing the ROI data mask. Can pass [] to have
%                 the function prompt the user to select a file.
%========================================================================================================

% Load ROI data file
if isempty(roiDataFile)
    roiDataFile = uigetfile([parentDir, '\*.mat'], 'Select an ROI data file');
end
load(fullfile(parentDir, roiDataFile));

% Load frame count log
individualVidFrameCounts = load(fullfile(parentDir, ['sid_', num2str(sid), '_frameCountLog.mat']));
frameCounts = [individualVidFrameCounts.frameCounts.nFrames];
framesPerTrial = mode(frameCounts);

% Calculate optic flow for each movie frame (unless file already exists)

    
disp('Calculating optic flow for behavior video...')
vidFile = fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials.mp4']);
 if isempty(dir(fullfile(parentDir, ['sid_', num2str(sid), '_optic_flow_data.mat'])))
     
    myVid = VideoReader(vidFile);
    frameCount = 0;
    opticFlow = opticalFlowFarneback;
    while hasFrame(myVid)
        
        frameCount = frameCount + 1;
        if ~mod(frameCount, 100)
            disp(['Reading frame ', num2str(frameCount), '...'])
        end
        
        currFrame = readFrame(myVid);
        currFrame = currFrame(:,:,1);
        
        % Calculate optic flow within each ROI
        currFrameFlowData = estimateFlow(opticFlow, currFrame);
        currFlowFly = currFrameFlowData.Magnitude;
        currFlowFly(~roiData) = nan;
        meanFlowMag(frameCount) = nanmean(nanmean(currFlowFly));
        
    end
    nFrames = frameCount;
    disp('Optic flow calculation complete')
    
    
    % Calculate frame times
    frameTimes = (1:1:nFrames) ./ FRAME_RATE;
    
    % Normalize flow magnitudes
    flyFlow = meanFlowMag;
    flyFlowNorm = flyFlow ./ max(flyFlow(2:end)); % First frame is artificially high so don't use that
    
    % Save optic flow data
    savefast(fullfile(parentDir, ['sid_', num2str(sid), '_optic_flow_data.mat']), 'flyFlowNorm', 'nFrames', 'frameTimes');
 
 else
     % Just load the existing data if it exists
    load(fullfile(parentDir, ['sid_', num2str(sid), '_optic_flow_data.mat']));

% TEMP
nFrames = 69899;
frameTimes = (1:1:nFrames) ./ FRAME_RATE;

 end
 
 
 
 % Recreate vidReader
myVid = VideoReader(vidFile);

% Create vidWriter
myVidWriter = VideoWriter(fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials_With_Optic_Flow.mp4']), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)

% Plot and write each frame
for iFrame = 1:nFrames
    
    if ~mod(iFrame, framesPerTrial)
        disp(['Plotting trial #' num2str(iFrame/framesPerTrial)])
    end
    currFrame = uint8(readFrame(myVid));
    ySize = size(currFrame, 1);
    xSize = size(currFrame, 2);
    
    % Create figure
    screenSize = get( groot, 'Screensize' );
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
    currFrameTime = iFrame * (1/FRAME_RATE);
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
    for iTrial = 1:(length(frameCounts)-1)
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
%     h.Position;
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);   
    
end

close(myVidWriter)

end