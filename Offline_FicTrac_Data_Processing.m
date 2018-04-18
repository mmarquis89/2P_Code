
parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_04_14_exp_1\_Movies\FicTracData';
FRAME_RATE = 25;
trialDuration = 20;
nFrames = FRAME_RATE * trialDuration;

xTick = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;

% Load FicTrac data from each trial into a single array
dataFiles = dir(fullfile(parentDir, '*sid*.dat'));
allData = []; selectData = []; droppedFrames = [];
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
    
    % Add to full data arrays
    allData(:,:,iTrial) = currData; % --> [frame, var, trial]
    selectData(:,:,iTrial) = currData(:, [1 15 16 17 19 23]); % --> Columns: [Frame Count, xPos, yPos, HD, Speed, Seq num]
    droppedFrames{iTrial} = currDroppedFrames;
end

%%
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


figure(1);clf;hold on
plot(filtfilt(kb, ka, dRotCam))
legend Right-1-X  Forward-2-Z Down-3-Y 
title dRotCam

figure(3);clf;hold on
plot(filtfilt(kb, ka, dRotLab))
legend X Y Z
title dRotLab

figure(4);clf;hold on
plot(absOrientCam(:, [1 3 2]))
legend Right-1-X  Down-3-Y Forward-2-Z
title absOrientCam

figure(5);clf;hold on
plot(absOrientLab)
legend X Y Z
title absOrientLab

figure(6);clf;hold on
plot(intXY)
legend X Y
title intXY

figure(7);clf;hold on
plot(intHD)
title intHD

figure(9);clf;hold on
plot(filtfilt(kb, ka, moveSpeed))
title moveSpeed

figure(10);clf;hold on
plot([intForwardMove, intSideMove])
legend Forward Side
title intForward+SideMove







%%
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
lims = 30;
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