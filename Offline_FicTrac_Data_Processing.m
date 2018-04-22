
%% BASIC DATA LOADING/PROCESSING

parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_04_05_exp_2\_Movies\FicTracData';
FRAME_RATE = 25;
trialDuration = 20;
nFrames = FRAME_RATE * trialDuration;

xTick = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;

% Load FicTrac data from each trial into a single array
dataFiles = dir(fullfile(parentDir, '*tid*.dat'));
allData = []; selectData = []; droppedFrames = []; csvData = [];
for iTrial = 1:numel(dataFiles)
    currData = csvread(fullfile(parentDir, dataFiles(iTrial).name));
    
    % Deal with dropped frames by repeating rows
    currDroppedFrames = [];
    for iFrame = 1:nFrames
        if size(currData, 1) >= iFrame
            if currData(iFrame, 1) ~= iFrame
                currData = [currData(1:iFrame, :); currData(iFrame, :); currData(iFrame + 1:end, :)];
                currData(iFrame, [1 23]) = iFrame; % Update frame and sequence nums
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
end

% Save .csv file of concatenated data for all trials
if ~exist(fullfile(parentDir, 'allTrials.csv'), 'file')
    csvwrite(fullfile(parentDir, 'allTrials.csv'), csvData);
end

% LP filter velocity data
rate = 2 * (11/25);
[kb, ka] = butter(2,rate);

%% VARIABLE PLOTTING SCRIPT

currTrial = 33;
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
plot(intHD)
legend HD hd
title intHD
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;
% 
% figure(9);clf;hold on
% plot(filtfilt(kb, ka, moveSpeed))
% title moveSpeed

figure(10);clf;hold on
plot([intForwardMove, intSideMove, intHD])
legend Forward Side HD
title intForward+SideMove

% MOVEMENT MAP
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
        lims = max(abs([ylim(axMove), xlim(axMove)]));
        xlim([-lims lims])
        ylim([-lims lims])
        axis equal; %axis square
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
    axMove = subaxis(3,6,[1:3 7:9 13:15], 'S', 0, 'M', M, 'PB', 0.03);
    axHD = subaxis(3,6,[4:6], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axDisp = subaxis(3,6,[10:12], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axVel = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    
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
    xlabel('Time (sec)')
    hold on

    % Write frame to video(s)
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
    
end

close(myVidWriter)


%% OLD PLOTTING CODE
% Get yaw velocity
dHD = [0; diff(unwrap(currData(:,6)))];

% Convert xy data from radians to mm and add yaw vel
r = 4.5;
mmData = [currData(:,1:2), currData(:, 3:5) * r, currData(:, 6), dHD]; % --> columns: [Frame Count, Seq num, xPos, yPos, Speed, HD, HD vel]           

% Smooth  data
smoothVelData = mmData(:, [5, 7]);
for iAxis = 1:size(smoothVelData, 2)
    smoothWin = 3;
    smoothVelData(:, iAxis) = smooth(smoothVelData(:, iAxis), smoothWin);  % --> [sample, axis, trial]
end

% LP filter velocity data
rate = 2 * (5/25);
[kb, ka] = butter(2,rate);
velData = filtfilt(kb, ka, smoothVelData); % --> columns: [xyVel, hdVel]

% DISPLACEMENT
figure(1); clf; hold on
plot(mmData(:, [3 4 6]));
plot(unwrap(mmData(:, 6)));
ax = gca;
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;
legend({'X', 'Y', 'HD'})
title('mmData')

% VELOCITY
figure(2); clf; hold on
plot(velData(:,1));
ax = gca;
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;
legend({'Speed (mm/sec)'})
title('velData')

% MOVEMENT MAP
figure(3); clf; hold on;
nSegs = 20;
vectorLen = floor(size(mmData, 1) / nSegs);
cm = jet(nSegs);
for iSeg = 1:nSegs
    if iSeg == 1
        pad = 1;
    else
        pad = 0;
    end
    currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
    plot(smooth(mmData(currSegInd, 3)), smooth(mmData(currSegInd, 4)), 'Color', cm(iSeg, :))   
end
lims = max(max(abs(mmData(:, [3 4]))));
xlim([-lims lims])
ylim([-lims lims])

% HD VELOCITY
figure(4); clf; hold on
plot(rad2deg(velData(:,2)));
ax = gca;
ax.XTick = xTick;
ax.XTickLabel = xTickLabels;
legend({'Yaw Speed (deg/frame)'})
title('hdvelData')

%% YEF DATA PROCESSING

sR = 20000;

xRad = xRaw * (2*pi/10);
yRad = yRaw * (2*pi/10);
hdRad = HDraw * (2*pi/2);
figure(1);clf;hold on;
plot(xRad); plot(yRad);

uwX = unwrap(xRad);
uwY = unwrap(yRad);
uwHD = unwrap(hdRad);
zeroX = uwX - uwX(1);
zeroY = uwY - uwY(1);
zeroHD = uwHD - uwHD(1);

smX = smooth(zeroX, 10);
smY = smooth(zeroY, 10);
smHD = smooth(zeroHD, 10);

dsX = downsample(smX, 400);
dsY = downsample(smY, 400);
dsHD = downsample(smHD, 400);


figure(2);clf;hold on
plot(xCut); plot(yCut);

xCut = dsX(1:10000) * 9;
yCut = dsY(1:10000) * 9;
hdCut = dsHD(1:10000);

xyMove = [xCut, yCut];
xyMoveDiff = [0,0 ; diff(xyMove, 1)];

figure(1);clf;hold on
currPoint = [0 0];
cm = jet(length(hdCut));
xyPos = [0 0];
for iTheta = 1:length(hdCut)
    disp(num2str(iTheta))
%    iTheta = 55;
   theta = dsHD(iTheta);
   R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
%    rotMat = R;
   rotXY = R * [xyMoveDiff(iTheta, 1); xyMoveDiff(iTheta, 2)];
   
   oldPoint = currPoint;
   currPoint = currPoint + rotXY';
   xyPos(iTheta, :) = currPoint;
   
   plot([oldPoint(1), currPoint(1)], [oldPoint(2), currPoint(2)], 'color', cm(iTheta,:))
end
   
yL = ylim;
xL = xlim;

lims = max([xL, yL]);
ylim([-lims, lims])
xlim([-lims, lims])


myVidWriter = VideoWriter(fullfile(parentDir, ['testVid.mp4']), 'MPEG-4');
myVidWriter.FrameRate = 50;
open(myVidWriter)
nSegs = 20;
h = figure(10);
h.Position = [50 50 800 800];
for iFrame = 1:length(hdCut)
    clf
    
    % Movement map plot
    hold on
    nSegs = 20;
    vectorLen = floor(length(hdCut) / nSegs);
    cm = jet(nSegs);
    for iSeg = 1:nSegs
        if iSeg == 1
            pad = 1;
        else
            pad = 0;
        end
        currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
        plot(xyPos(currSegInd, 1), xyPos(currSegInd, 2), 'Color', cm(iSeg, :))
    end
%     lims = max(abs([ylim, xlim]));
%     xlim([-lims lims])
%     ylim([-lims lims])

    xlim([-20 5]);
    ylim([-5 20]);
    axis equal;
    legend({'2D movement (mm)'}, 'FontSize', 11)
    x = xyPos(iFrame, 1);
    y = xyPos(iFrame, 2);
    [arrowVec(1), arrowVec(2)] = pol2cart(hdCut(iFrame)+(pi * 1.5), 0.1);
    
    arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2]);
    
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
end
close(myVidWriter)
