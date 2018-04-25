%% LOAD IMAGING METADATA

% Load .mat file containing trial data
[myData, m] = load_imaging_metadata();

%% BASIC DATA LOADING/PROCESSING

parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_04_14_exp_2\_Movies\FicTracData';
FRAME_RATE = 25;
trialDuration = 20;
nFrames = FRAME_RATE * trialDuration;

xTick = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;

% Load FicTrac data from each trial into a single array
dataFiles = dir(fullfile(parentDir, '*tid*.dat'));
allData = []; selectData = []; droppedFrames = []; csvData = []; resets = [];
for iTrial = 1:numel(dataFiles)
    currData = csvread(fullfile(parentDir, dataFiles(iTrial).name));
    
    % Deal with dropped frames by repeating rows
    currDroppedFrames = [];
    for iFrame = 1:nFrames
        if size(currData, 1) >= iFrame
            if currData(iFrame, 1) ~= iFrame
                currData = [currData(1:iFrame-1, :); currData(iFrame-1, :); currData(iFrame:end, :)];
                currData(iFrame, [1 23]) = [iFrame, currData(iFrame, 23) + 1]; % Update frame and sequence nums
                currDroppedFrames(end + 1) = iFrame; % --> [trial, framelist]
            end
        else
            currData(iFrame, :) = currData(end, :);
            currDroppedFrames(end + 1) = iFrame; % --> [trial, framelist]
        end
    end
    
    
    currDataCSV = [ones(size(currData, 1), 1) * iTrial, currData];
    
    % Add to full data arrays
    csvData = [csvData; currDataCSV];
    allData(:,:,iTrial) = currData; % --> [frame, var, trial]
    selectData(:,:,iTrial) = currData(:, [1 15 16 17 19 23]); % --> Columns: [Frame Count, xPos, yPos, HD, Speed, Seq num]
    droppedFrames{iTrial} = currDroppedFrames;
    resets(iTrial) = sum(currData(:,23) == 1) - 1;
end

% Save .csv file of concatenated data for all trials
if ~exist(fullfile(parentDir, 'allTrials.csv'), 'file')
    csvwrite(fullfile(parentDir, 'allTrials.csv'), csvData);
end

% LP filter for velocity data
rate = 2 * (6/25);
[kb, ka] = butter(2,rate);

% Split into individual components

frameCounter = squeeze(allData(:, 1, :));
dRotCam = squeeze(allData(:, [2 3 4], :));          % [X, Y, Z]
dRotError = squeeze(allData(:, 5, :));
dRotLab = squeeze(allData(:, [6 7 8], :));          % [X, Y, Z]
absOrientCam = squeeze(allData(:, [9 10 11], :));   % [X, Y, Z]
absOrientLab = squeeze(allData(:, [12 13 14], :));  % [X, Y, Z]
intXY = squeeze(allData(:, [15 16], :));            % [X, Y]
intHD = squeeze(allData(:, 17, :));
moveDirLab = squeeze(allData(:, 18, :));
moveSpeed = squeeze(allData(:, 19, :));
intForwardMove = squeeze(allData(:, 20, :));
intSideMove = squeeze(allData(:, 21, :));
timestamp = squeeze(allData(:, 22, :));
seqNum = squeeze(allData(:, 23, :)); 


%% % Overlay 2D movement data for all trials
currData = selectData(:, [1 6 2 3 5 4],:); % re-order columns to be [Frame Count, Seq num, xPos, yPos, Speed, HD]

figure(13);clf;hold on; ax = gca();
intXY = currData(:, [3 4], :); % --> [frame, var, trial]
HD = squeeze(currData(:, 6, :));        % --> [frame, trial]
nFrames = size(intXY, 1);
cm = jet(nFrames);
cm2 = parula(nFrames);
startTime = 14;
frameRate = 25;
startFrame = startTime * frameRate;

% Movement map
maxes = 0;
for iTrial = 1:size(intXY, 3)%[16 80 28]
    
    currXY = intXY(:, :, iTrial);
    smoothX = smooth(currXY(:, 1), 7);
    smoothY = smooth(currXY(:, 2), 7);
    
    meanHD = mod(smooth(unwrap(HD(:, iTrial)), 11), (2*pi));
    theta = -meanHD(startFrame) + pi/2;
    
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    smoothXYRot = R * [smoothX'; smoothY'];
    
    xData = smoothXYRot(1, startFrame:end);
    yData = smoothXYRot(2, startFrame:end);
    
    xData = xData - xData(1);
    yData = yData - yData(1);
    
    maxes(end + 1) = max([xData, yData]);
    
    cData = permute(repmat(cm(startFrame:end, :), 1, 1, 2), [1 3 2]);
    surface('XData', [xData' xData'], ...
        'YData', [yData' yData'], ...
        'ZData', [zeros(numel(xData), 2)], ...
        'CData', cData, ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'marker', 'none');
end
axis equal
% currXY = intXY(startFrame:end, :, :);
lims = 1.1 * max(abs(maxes));
xlim([-lims lims])
ylim([-lims lims])
legend({'2D movement (mm)'}, 'FontSize', 11)




%% Separate data according to trial type

% Separate trials

yL = [0 0.05];
s = myData.stimSepTrials;
trialGroups = [[s.OdorA + 2 * s.OdorB] .* myData.goodTrials];
legendStr = {'EtOH', 'ACV', 'No stim'};%

trialGroups(1:floor(nTrials/3)) = 0;
trialGroups(floor(nTrials/3):floor(2*nTrials/3)) = 0;
% trialGroups(floor(2*nTrials/3):end) = 0;
% 
% trialGroups = ones(1, nTrials) .* myData.goodTrials;
% trialGroups(floor(nTrials/3):end) = 2;
% trialGroups(floor(2*nTrials/3):end) = 3;
% % trialGroups(~(s.OdorB)) = 0;
% legendStr = {'Early trials', 'Mid trials', 'Late trials'};%

% postOdorA = [0, s.OdorA(1:end-1)];
% postOdorB = [0, s.OdorB(1:end-1)];
% postNoStim = [0, s.NoStim(1:end-1)];
% trialGroups = [(postOdorA & s.NoStim) + 2 * (postOdorB & s.NoStim) + 3 * (postNoStim & s.NoStim)] .* myData.goodTrials;
% legendStr = {'EtOH', 'ACV', 'No stim'};%

% trialGroups = ones(1, nTrials) .* myData.goodTrials;
% trialGroups(floor(nTrials/3):end) = 2;
% trialGroups(floor(2*nTrials/3):end) = 3;
% trialGroups(~(postNoStim)) = 0;
% legendStr = {'Early trials', 'Mid trials', 'Late trials'};%



% Overlay traces separated by trialGroups

currData = selectData(:, [1 6 2 3 5 4],:); % re-order columns to be [Frame Count, Seq num, xPos, yPos, Speed, HD]

figure(13);clf;hold on; ax = gca();
intXY = currData(:, [3 4], :); % --> [frame, var, trial]
HD = squeeze(currData(:, 6, :));        % --> [frame, trial]
nFrames = size(intXY, 1);
nTrials = size(intXY, 3);
startTime = 11;
frameRate = 25;
startFrame = startTime * frameRate;

% cm = [];
% cm{1} = hot(nFrames);
% cm{2} = cool(nFrames);
% cm{3} = winter(nFrames);
cm = [rgb('blue'); rgb('red'); rgb('green')];


% Movement map
maxes = 0; legendPlots = [0 0 0]; legendObj = [];
for iTrial = 1:nTrials%[16 80 28]
    
    if trialGroups(iTrial)
        currXY = intXY(:, :, iTrial);
        smoothX = smooth(currXY(:, 1), 7);
        smoothY = smooth(currXY(:, 2), 7);
        
        meanHD = mod(smooth(unwrap(HD(:, iTrial)), 11), (2*pi));
        theta = -meanHD(startFrame) + pi/2;
        
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        smoothXYRot = R * [smoothX'; smoothY'];
        
        xData = smoothXYRot(1, startFrame:end);
        yData = smoothXYRot(2, startFrame:end);
        
        xData = xData - xData(1);
        yData = yData - yData(1);
        
        maxes(end + 1) = max([xData, yData]);
        
        plt = plot(xData, yData, 'Color', cm(trialGroups(iTrial), :));
%         cData = permute(repmat(cm{trialGroups(iTrial)}(startFrame:end, :), 1, 1, 2), [1 3 2]);
%         surface('XData', [xData' xData'], ...
%             'YData', [yData' yData'], ...
%             'ZData', [zeros(numel(xData), 2)], ...
%             'CData', cData, ...
%             'FaceColor', 'none', 'EdgeColor', 'interp', 'marker', 'none');

        % Save one plot line from each group to use in legend
        if ~legendPlots(trialGroups(iTrial))
            legendObj(trialGroups(iTrial)) = plt;
            legendPlots(trialGroups(iTrial)) = 1;
        end
    end
end%iTrial

axis equal
% currXY = intXY(startFrame:end, :, :);
lims = 1 * max(abs(maxes));
xlim([-lims lims])
ylim([-lims lims])
legend(legendObj, legendStr, 'FontSize', 12)
% legend({'2D movement (mm)', 'b', 'c'}, 'FontSize', 11)


% Plot mean locomotion speed at each time throughout trial

speedData = squeeze(selectData(:, 5, :)); % --> [frame, trial]
cm = [rgb('blue'); rgb('red'); rgb('green')];
dHD = abs(diff(unwrap(intHD, [], 1), 1));

% LP filter for velocity data
rate = 2 * (12/25);
[kb, ka] = butter(2,rate);

figure(14);clf;hold on; ax = gca;
for iGroup = 1:numel(unique(trialGroups(trialGroups > 0)))
    plot(filtfilt(kb, ka, smooth(mean(speedData(:, trialGroups == iGroup), 2), 11)), 'Color', cm(iGroup, :));
end
legend(legendStr)
ax.XTick = xTick;
ax.XTickLabels = xTickLabels;
shadeFrames = [10 12] * FRAME_RATE;
plot_stim_shading(shadeFrames, 'Axes', ax);
title('Mean XY velocity')
xlabel('Time (sec)')
ylim(yL)

% LP filter for velocity data
rate = 2 * (12/25);
[kb, ka] = butter(2,rate);

figure(15);clf;hold on; ax = gca;
for iGroup = 1:numel(unique(trialGroups(trialGroups > 0)))
    plot(filtfilt(kb, ka, movmean(mean(dHD(:, trialGroups == iGroup), 2), 11)), 'Color', cm(iGroup, :));
end
legend(legendStr)
ax.XTick = xTick;
ax.XTickLabels = xTickLabels;
shadeFrames = [10 12] * FRAME_RATE;
plot_stim_shading(shadeFrames, 'Axes', ax);
title('Mean rate of heading change')
xlabel('Time (sec)')
ylim(yL)

%% VARIABLE PLOTTING SCRIPT

currTrial = 28;
currData = selectData(:, [1 6 2 3 5 4],currTrial); % re-order columns to be [Frame Count, Seq num, xPos, yPos, Speed, HD]
currAllData = allData(:, :, currTrial); % --> [frame, variable]

% Split into individual components
frameCounter = currAllData(:, 1);
dRotCam = currAllData(:, [2 3 4]);          % [X, Y, Z]
dRotError = currAllData(:, 5);
dRotLab = currAllData(:, [6 7 8]);          % [X, Y, Z]
absOrientCam = currAllData(:, [9 10 11]);   % [X, Y, Z]
absOrientLab = currAllData(:, [12 13 14]);  % [X, Y, Z]
intXY = currAllData(:, [15 16]);            % [X, Y]
intHD = currAllData(:, 17);
moveDirLab = currAllData(:, 18);
moveSpeed = currAllData(:, 19);
intForwardMove = currAllData(:, 20);
intSideMove = currAllData(:, 21);
timestamp = currAllData(:, 22);
seqNum = currAllData(:, 23);

% figure(1);clf;hold on
% plot(filtfilt(kb, ka, dRotCam))
% legend Right-1-X  Forward-2-Z Down-3-Y
% title dRotCam
%
% figure(3);clf;hold on
% plot(filtfilt(kb, ka, dRotLab))
% legend X Y Z
% title dRotLab
%
% figure(4);clf;hold on
% plot(absOrientCam(:, [1 3 2]))
% legend Right-1-X  Down-3-Y Forward-2-Z
% title absOrientCam
%
% figure(5);clf;hold on
% plot(absOrientLab)
% legend X Y Z
% title absOrientLab

figure(6);clf;hold on; ax = gca();
plot(intXY)
legend X Y x y
title intXY
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;

figure(7);clf;hold on; ax = gca();
uwHD = unwrap(intHD);
smHD = smooth(uwHD, 1);
plot(mod(smHD, (2*pi)))
legend HD hd
title intHD
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;
%
figure(9);clf;hold on
plot(filtfilt(kb, ka, smooth(moveSpeed, 7)))
title moveSpeed

% figure(10);clf;hold on
% plot([intForwardMove, intSideMove, intHD])
% legend Forward Side HD
% title intForward+SideMove

% Movement map
figure(11); clf; hold on;
nSegs = 20;
vectorLen = floor(size(intXY, 1) / nSegs);
cm = jet(nSegs);
for iSeg = 1:nSegs
    if iSeg == 1
        pad = 1;
    else
        pad = 0;
    end
    currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
    plot(smooth(intXY(currSegInd, 1),3), smooth(intXY(currSegInd, 2),3), 'Color', cm(iSeg, :))
end
lims = 1.1 * max(abs(intXY(:)));
xlim([-lims lims])
ylim([-lims lims])

% Anvil annotations
figure(12); clf; hold on
plot(behaviorAnnotArr(currTrial, :), '-*', 'Color', 'k')
cm = [rgb('Navy'); rgb('Cyan'); rgb('maroon'); rgb('Gold')];
for iPlot = 1:4
    currData = behaviorAnnotArr(currTrial, :);
    currData(currData ~= iPlot - 1) = nan;
    plot(currData, '*', 'color', cm(iPlot, :))
end
ylim([-1 4])
ax = gca;
ax.YTick = [0 1 2 3];
ax.YTickLabel = {'Quiescence', 'Locomotion', 'Grooming', 'IsolatedMovement'};

%% CREATE FICTRAC + BEHAVIOR VIDS
vidDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_04_05_exp_2\_Movies'

vidFiles = dir(fullfile(vidDir, ['sid*tid*.mp4']));

for iTrial = 1:size(selectData, 3)
    
    vidName = vidFiles(iTrial).name;
    disp(vidName)
    vidFile = fullfile(vidDir, vidName);
    currData = selectData(:, [1 6 2 3 5 4], iTrial); % --> [Frame Count, Seq num, xPos, yPos, Speed, HD]
    
    % Convert xy data from radians to mm
    r = 4.5;
    mmData = [currData(:,1:2), currData(:, 3:5) * r, currData(:, 6)];
    
    % Smooth velocity data
    smoothVelData = [];
    for iAxis = 1:size(mmData, 2)
        smoothWin = 3;
        smoothVelData(:, iAxis) = smooth(mmData(:, iAxis), smoothWin);  % --> [sample, axis, trial]
    end
    % LP filter velocity data
    rate = 2 * (11/25);
    [kb, ka] = butter(2,rate);
    velData = filtfilt(kb, ka, smoothVelData);
    
    % Heading direction data
    HD = currData(:,end);
    
    % Create VidReader and VidWriter
    myVid = VideoReader(vidFile);
    myVidWriter = VideoWriter(fullfile(parentDir, ['Combined_', vidName]), 'MPEG-4');
    myVidWriter.FrameRate = FRAME_RATE;
    open(myVidWriter)
    
    h = figure(10);
    for iFrame = 1:nFrames
        
        currFrame = uint8(readFrame(myVid));
        %
        clf
        ySize = size(currFrame, 1);
        xSize = size(currFrame, 2);
        
        % Create figure
        screenSize = get( groot, 'Screensize' );
        h.Color = [1 1 1];
        h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];
        xyRatio = h.Position(3) / h.Position(4);
        
        % Create axes
        clf
        M = 0.006;
        P = 0.00;
        axVid = subaxis(3,6,[1 2 7 8], 'S', 0, 'M', M, 'PB', 0.05);
        axMove = subaxis(3,6,[3 4 9 10], 'S', 0, 'M', M, 'PB', 0.06, 'PR', 0.01);
        axHD = subaxis(3,6,[5 6 11 12], 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.01);
        axDisp = subaxis(3,6,[13:15], 'S', 0, 'M', M, 'PB', 0.05, 'PR', 0.01, 'PL', 0.008);
        axVel = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.02);
        
        % Movie frame plot
        imshow(currFrame, [], 'Parent', axVid);
        axis off
        axVid.XTickLabel = [];
        
        % Movement map plot
        axes(axMove);
        hold on
        nSegs = trialDuration;
        vectorLen = floor(size(mmData, 1) / nSegs);
        cm = jet(nSegs);
        for iSeg = 1:nSegs
            if iSeg == 1
                pad = 1;
            else
                pad = 0;
            end
            currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
            plot(axMove, mmData(currSegInd, 3), mmData(currSegInd, 4), 'Color', cm(iSeg, :))
        end
        axis equal; %axis square
        lims = 1.1 * max(abs([ylim(axMove), xlim(axMove)]));
        if lims < 1
            lims = 1;
        end
        xlim([-lims lims])
        ylim([-lims lims])
        
        legend({'2D movement (mm)'}, 'FontSize', 11)
        x = mmData(iFrame, 3);
        y = mmData(iFrame, 4);
        [arrowVec(1), arrowVec(2)] = pol2cart(HD(iFrame)+(pi * 1.5), 0.5);
        %         [testVec(1), testVec(2)] = pol2cart(test(iFrame) + (pi * 1.5), 0.5);
        arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2], 'FaceColor', 'green');
        %         arrow([x - testVec(1)/2, y-testVec(2)/2], [x + testVec(1)/2, y + testVec(2)/2], 'FaceColor', 'red');
        
        % Head direction plot
        axes(axHD)
        plt = polarplot([HD(iFrame), HD(iFrame)], [0 1], '-', 'LineWidth', 5, 'Color', 'k');
        axHD = plt.Parent;
        axHD.FontSize = 12;
        axHD.ThetaZeroLocation = 'bottom';
        axHD.RTickLabel = [];
        axHD.RTick = [];
        axes(axHD);
        hold on
        polarplot(linspace(0, pi * 2, 100), ones(100,1), 'linewidth', 0.5, 'Color', 'k');
        
        % Displacement plot
        axes(axDisp);
        hold on
        plot(mmData(:, [3 4]));
        axDisp.XTick = xTick;
        axDisp.XTickLabel = xTickLabels;
        legend({'X disp (mm)', 'Y disp (mm)'}, 'FontSize', 11, 'Location', 'NorthWest')
        xlabel('Time (sec)')
        plot(iFrame, mmData(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        plot(iFrame, mmData(iFrame, 4), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        
        % Velocity plot
        axes(axVel)
        plot(velData(:,5));
        axVel.XTick = xTick;
        axVel.XTickLabel = xTickLabels;
        legend({'Speed (mm/sec)'}, 'FontSize', 11)
        xlabel('Time (sec)')
        hold on
        plot(iFrame, velData(iFrame, 5), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        
        % Write frame to video(s)
        writeFrame = getframe(h);
        writeVideo(myVidWriter, writeFrame);
    end%iFrame
    
    close(myVidWriter)
end

%% CREATE VID WITH TRIAL-BY-TRIAL FICTRAC PLOTS

% Create VidWriter
myVidWriter = VideoWriter(fullfile(parentDir, 'FicTrac_Trial_Plots.avi'));
myVidWriter.FrameRate = 1;
open(myVidWriter)

for iTrial = 1:size(selectData, 3)
    
    disp(num2str(iTrial))
    currData = selectData(:, [1 6 2 3 5 4], iTrial); % --> [Frame Count, Seq num, xPos, yPos, Speed, HD]
    
    % Convert xy data from radians to mm
    r = 4.5;
    mmData = [currData(:,1:2), currData(:, 3:5) * r, currData(:, 6)];
    
    % Smooth velocity data
    smoothVelData = [];
    for iAxis = 1:size(mmData, 2)
        smoothWin = 3;
        smoothVelData(:, iAxis) = smooth(mmData(:, iAxis), smoothWin);  % --> [sample, axis, trial]
    end
    % LP filter velocity data
    rate = 2 * (11/25);
    [kb, ka] = butter(2,rate);
    velData = filtfilt(kb, ka, smoothVelData);
    
    % Heading direction data
    HD = currData(:,end);
    
    % Create figure
    h = figure(100); clf
    screenSize = get( groot, 'Screensize' );
    h.Color = [1 1 1];
    h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];
    
    % Create axes
    M = 0.01;
    P = 0.00;
    axMove = subaxis(3,6,[1:3 7:9], 'S', 0, 'M', M, 'PB', 0.05);
    axAnnot = subaxis(3,6,[4:6], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axVel = subaxis(3,6,[10:12], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axDisp = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axHD = subaxis(3,6,[13:15], 'S', 0, 'M', M, 'PB', 0.05, 'PR', 0.02);
    
    % Movement map plot
    axes(axMove);
    hold on
    nSegs = trialDuration;
    vectorLen = floor(size(mmData, 1) / nSegs);
    cm = jet(nSegs);
    for iSeg = 1:nSegs
        if iSeg == 1
            pad = 1;
        else
            pad = 0;
        end
        currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
        plot(axMove, smooth(mmData(currSegInd, 3), 3), smooth(mmData(currSegInd, 4), 3), 'Color', cm(iSeg, :))
    end
    axis equal;
    lims = 1.1 * (max(abs([ylim(axMove), xlim(axMove)])));
    if lims < 1
        lims = 1;
    end
    xlim([-lims lims])
    ylim([-lims lims])
    legend({'2D movement (mm)'}, 'FontSize', 11)
    x = mmData(iFrame, 3);
    y = mmData(iFrame, 4);
    
    % Head direction plot
    axes(axHD)
    plot(HD);
    legend({'Heading (rad)'}, 'FontSize', 11, 'Location', 'NorthWest');
    axHD.XTick = xTick;
    axHD.XTickLabel = xTickLabels;
    xlabel('Time (sec)')
    
    % Displacement plot
    axes(axDisp);
    hold on
    plot(mmData(:, [3 4]));
    axDisp.XTick = xTick;
    axDisp.XTickLabel = xTickLabels;
    legend({'X disp (mm)', 'Y disp (mm)'}, 'FontSize', 11, 'Location', 'NorthWest')
    xlabel('Time (sec)')
    
    % Velocity plot
    axes(axVel)
    plot(velData(:,5));
    axVel.XTick = xTick;
    axVel.XTickLabel = xTickLabels;
    legend({'Speed (mm/sec)'}, 'FontSize', 11)
    hold on
    
    % Anvil annotations
    axes(axAnnot)
    hold on
    plot(behaviorAnnotArr(iTrial, :), '-*', 'Color', 'k')
    cm = [rgb('Navy'); rgb('Cyan'); rgb('maroon'); rgb('Gold')];
    for iPlot = 1:4
        currData = behaviorAnnotArr(iTrial, :);
        currData(currData ~= iPlot - 1) = nan;
        plot(currData, '*', 'color', cm(iPlot, :))
    end
    ylim([-1 4])
    axAnnot.YTick = [0 1 2 3];
    axAnnot.YTickLabel = {'Quiescence', 'Locomotion', 'Grooming', 'IsolatedMovement'};
    axAnnot.XTick = xTick;
    axAnnot.XTickLabel = xTickLabels;
    
    % Write frame to video(s)
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
end

close(myVidWriter)