function make_optic_flow_vid(sid, parentDir, FRAME_RATE, vidFile)







% Load ROI data file
if isempty(vidFile)
    roiDataFile = uigetfile([parentDir, '\*.mat'], 'Select an ROI data file');
else
    roiDataFile = vidFile;
end
load(fullfile(parentDir, roiDataFile));

% Load frame count log
individualVidFrameCounts = load(fullfile(parentDir, ['sid_', num2str(sid), '_frameCountLog.mat']));
frameCounts = [individualVidFrameCounts.frameCounts.nFrames];

% % Load and concatenate optic flow data
% flowDataFiles = uigetfile([parentDir, '\*.mat'], 'Select optic flow data files for all trials', 'MultiSelect', 'on');
% meanFlowMags = [];
% frameCounts = [];
% for iTrial = 1:length(flowDataFiles)
%     currFlowData = load(fullfile(parentDir, flowDataFiles{iTrial}));
%     frameCounts(iTrial) = size(currFlowData.rawFlowMag, 3);
%     
%     % Get mean flow mag within each ROI for the current trial
%     flowMag = zeros(size(currData.rawFlowMag, 3), 2);
%     for iFrame = 1:size(rawFlowCat, 3)
%         currMags = currData.rawFlowMag(:,:,iFrame);
%         
%         currMagsFly = currMags;
%         currMagsFly(~roiData(:,:,1)) = nan;
%         flowMag(iFrame, 1) = nanmean(nanmean(currMagsFly));
%         
%         currMagsWasher = currMags;
%         currMagsWasher(~roiData(:,:,2)) = nan;
%         flowMag(iFrame, 2) = nanmean(nanmean(currMagsWasher));
%         disp(num2str(iFrame))
%     end
%     meanFlowMags = cat(1, meanFlowMags, flowMag); 
% end


% Calculate optic flow for each movie frame
disp('Calculating optic flow for behavior video...')
vidFile = fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials.mp4']);
% nFrames = count_frames(vidFile);
myVid = VideoReader(vidFile); %'sid_', num2str(sid), '_AllTrials
frameCount = 0;
opticFlow = opticalFlowFarneback;
while hasFrame(myVid)
    
    frameCount = frameCount + 1;
    if ~mod(frameCount, 100)
        disp(['Reading frame ', num2str(frameCount), '...'])
    end
    
    currFrame = readFrame(myVid);
    currFrame = currFrame(:,:,1);
%     if frameCount == 1
%         frameSize = size(currFrame);
%         myMovie = uint8(zeros([frameSize, nFrames]));
%     end
%     myMovie(:,:,frameCount) = uint8(currFrame);
    
    % Calculate optic flow within each ROI
    currFrameFlowData = estimateFlow(opticFlow, currFrame);
    currFlowFly = currFrameFlowData.Magnitude;
    currFlowFly(~roiData(:,:,1)) = nan;
    meanFlowMag(frameCount, 1) = nanmean(nanmean(currFlowFly));
    currFlowWasher = currFrameFlowData.Magnitude;
    currFlowWasher(~roiData(:,:,2)) = nan;
    meanFlowMag(frameCount, 2) = nanmean(nanmean(currFlowWasher));
   
end
nFrames = frameCount;
disp('Optic flow calculation complete')


% Calculate frame times
frameTimes = (1:1:nFrames) ./ FRAME_RATE;

% Normalize flow magnitudes
flyFlow = meanFlowMag(:,1);
flyFlowNorm = flyFlow ./ max(flyFlow(2:end)); % First frame is artificially high so don't use that
washerFlow = meanFlowMag(:,2);
washerFlowNorm = washerFlow ./ max(washerFlow(2:end));

% Recreate vidReader
myVid = VideoReader(vidFile); %'sid_', num2str(sid), '_AllTrials

% Create vidWriter
myVidWriter = VideoWriter(fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials_With_Optic_Flow.mp4']), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)

% Plot and write each frame
for iFrame = 1:nFrames
    
    currFrame = uint8(readFrame(myVid));
    ySize = size(currFrame, 1);
    xSize = size(currFrame, 2);
    
    % Create figure
    screenSize = get( groot, 'Screensize' );
    h = figure(10);clf
    h.OuterPosition = [50 50 (xSize) screenSize(4)-50];
    
    % Movie frame plot
    ax = axes('Units', 'Normalized', 'Position', [0 0.3 1 0.7]);
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
    for iTrial = 1:length(frameCounts)
        trialBoundTimes(iTrial) = frameTimes(runningCount);
        runningCount = runningCount + frameCounts(iTrial);
    end
    
    axes('Units', 'Pixels', 'Position', [0 0 ax.Position(3) ax.Position(2)]);
    hold on
    plot(frameTimes(2:end), smooth(flyFlowNorm(2:end), 5), 'color', 'm');    % Plot fly movmement ROI flow
    plot(frameTimes(2:end), smooth(washerFlowNorm(2:end), 5), 'color', 'b'); % Plot washer ROI flow
    plot([currFrameTime, currFrameTime], ylim(), 'LineWidth', 2, 'color', 'r')
    for iTrial = 1:length(trialBoundTimes)
       plot([trialBoundTimes(iTrial), trialBoundTimes(iTrial)], ylim(), 'LineWidth', 2, 'color', 'k') 
    end
    set(gca, 'xticklabel', []);
    xlim(xL);
    ylabel('Optic flow (au)');
    lgd = legend('Fly movmement', 'Washer movement');
    lgd.LineWidth = 1;
    
    drawnow()
    % Write frame to video
    h.Position
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);   
    
end

close(myVidWriter)

end