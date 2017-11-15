%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA

  
    % Load .mat file containing trial metadata
    myData = load_imaging_data();

for iFold = 1 
    expDate = myData.expDate;
    nPlanes = myData.nPlanes;
    nVolumes = myData.nVolumes;
    refImg = myData.refImg;
    nFrames = myData.nFrames;
    nTrials = myData.nTrials;
    stimTypes = myData.stimTypes;
    trialDuration = myData.trialDuration;
    volumeRate = myData.volumeRate;
    volFrames = myData.volFrames;
    
    myData.stimDuration = [trialDuration(1), trialDuration(2)]; stimDuration = myData.stimDuration; % [stim start time, stim length] in seconds
    myData.stimStart = myData.stimDuration(1); stimStart = myData.stimStart;
    myData.stimEnd = sum(myData.stimDuration); stimEnd = myData.stimEnd;
    
    % Create hardcoded parameters
    myData.ROIdata = [];
    myData.MAX_INTENSITY = 2000; MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots
    myData.FRAME_RATE = 25; FRAME_RATE = 25; % This is the frame rate of the behavior video, not the GCaMP imaging
 
end%iFold
   
%% LOAD ROI DATA

parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sid, '\Analysis'];
fileName = 'ROI_Data.mat';

load(fullfile(parentDir, fileName));
myData.ROIdata = ROIdata;
 

%% PLOT AND SAVE VISUALIZATION OF BEHAVIOR DATA ANNOTATIONS

saveFigs = 0;
trialGroups = []; %myData.stimSepTrials.windTrials + 2 * (~myData.stimSepTrials.windTrials);

for iFold = 1
    
% Create array of annotation data (row = trial, col = frame)
annotationArr = [];
annotTrials = 1:myData.nTrials;
for iTrial = annotTrials(myData.goodTrials) % Skip any trials with dropped frames 
   annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.actionNums; %--> [trial, frame]
end

%----------Plot 2D data with seconds on X-axis----------
titleStr = [regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Behavior Summary']; % regex to add escape characters
[~, ~, f] = plot_behavior_summary_2D(myData, annotationArr, [], titleStr, trialGroups);

%----- Plot 1D trial-averaged movement data -----

% Average across trials for wind and control trials
annotArrLog = annotationArr ~= 0;
annotArrSum_wind = sum(annotArrLog(myData.stimSepTrials.windTrials,:), 1);
annotArrSum_noWind = sum(annotArrLog(~myData.stimSepTrials.windTrials,:), 1);

% Plot stimulus trial data
h = figure(2); clf;
h.Position = [100 50 1600 950];
subplot(2,1,1); ax1 = gca();
[~, ~, ~] = plot_behavior_summary_1D(myData, annotArrSum_wind, ax1, 'Wind Trials');

% Add shading during stimulus presentation
yL = ylim();
rectPos = [stimStart*FRAME_RATE, yL(1), (stimEnd-stimStart)*FRAME_RATE, diff(yL)]; % [x y width height]
rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
ylim(yL);

% Plot control trial data
subplot(2,1,2); ax2 = gca();
[~, ~, ~] = plot_behavior_summary_1D(myData, annotArrSum_noWind, ax2, 'Control Trials');

suptitle([regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Summed Movement Frames']) % regex to add escape characters

% ----------Save data as a .fig file and a .png file for each figure----------
if saveFigs
    % Create analysis directory if necessary
    saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sid, '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [expDate, '_Behavior_Summary.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [expDate, '_Behavior_Summary.png']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [expDate, '_Summed_Movement_Frames.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [expDate, '_Summed_Movement_Frames.png']), 'file') ~= 0
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        savefig(f, fullfile(saveDir, [expDate, '_Behavior_Summary']));
        export_fig(fullfile(saveDir, [expDate, '_Behavior_Summary']), '-png', f);
        savefig(h, fullfile(saveDir, [expDate, '_Summed_Movement_Frames']));
        export_fig(fullfile(saveDir, [expDate, '_Summed_Movement_Frames']), '-png', h);
    end
end%if
end%iFold

%% PLOT AND SAVE SUMMARY OF BALL STOPPING ANNOTATIONS

saveFig = 1;

for iFold = 1
    
% Create array of annotation data (row = trial, col = frame)
annotationArr = [];
annotTrials = 1:myData.nTrials;
for iTrial = annotTrials(myData.goodTrials) % Skip any trials with dropped frames 
   annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.ballStopNums; %--> [trial, frame]
end

%----------Plot 2D data with seconds on X-axis----------
f = figure(5);clf
f.Position = [200 100 1120 840];
ax = gca;
titleStr = [regexprep(expDate, '_(?<num>..)', '\\_$<num>'), ' ball-stopping Summary']; % regex to add escape characters
[~, ~, ~] = plot_behavior_summary_2D(myData, annotationArr, ax, titleStr, []);


% ----------Save data as a .fig file and a .png file for each figure----------
if saveFig
    % Create analysis directory if necessary
    saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sid, '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [expDate, '_Ball_Stopping_Summary.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [expDate, '_Ball_Stopping_Summary.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        savefig(f, fullfile(saveDir, [expDate, '_Ball_Stopping_Summary']));
        export_fig(fullfile(saveDir, [expDate, '_Ball_Stopping_Summary']), '-png', f);

    end
end%if

end%iFold

%% =================================================================================================
%            WIND STIMULUS ANALYSES                                   
%%%=================================================================================================
for iFold = 1
    %% CALCULATE MEAN dF/F FOR WIND STIM ONSET AND OFFSET

    combineStimTypes = 0;
    baselineLength = stimDuration(2);
    respLength = stimDuration(2);

    for iFold = 1;
        % Divide data into different stim types
        stimTypeData = sep_stim_types(myData, combineStimTypes); % --> [stimType]{x, y, plane, volume, trial}

        % Figure out which volumes are needed for each period
        onsetBaselineVols = ceil((stimStart - baselineLength) * volumeRate):floor(stimStart * volumeRate);
        offsetBaselineVols = ceil((stimEnd - baselineLength) * volumeRate):floor(stimEnd * volumeRate);
        onsetRespVols = ceil(stimStart * volumeRate):floor(stimEnd * volumeRate);
        offsetRespVols = ceil(stimEnd * volumeRate):floor((stimEnd + respLength) * volumeRate);

        % Initialize arrays
        onsetBaseline = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(onsetBaselineVols), nTrials, length(stimTypeData)]);
        offsetBaseline = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(offsetBaselineVols), nTrials, length(stimTypeData)]);
        onsetResp = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(onsetRespVols), nTrials, length(stimTypeData)]);
        offsetResp = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(offsetRespVols), nTrials, length(stimTypeData)]);
        onsetDff = zeros(size(onsetBaseline(:,:,:,:,1,:)));
        offsetDff = zeros(size(offsetBaseline(:,:,:,:,1,:)));
        onsetDffAvg = zeros(size(onsetBaseline(:,:,:,1,1,:)));
        offsetDffAvg = zeros(size(offsetBaseline(:,:,:,1,1,:)));

        for iStim = 1:length(stimTypeData)

            % Extract baseline and response period volumes from the each stim type
            disp('Extracting volumes...')
            onsetBaseline(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, onsetBaselineVols,:);                           % --> [x, y, plane, volume, trial, stimType]
            offsetBaseline(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, offsetBaselineVols,:);                         % --> [x, y, plane, volume, trial, stimType]
            onsetResp(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, onsetRespVols,:);                                   % --> [x, y, plane, volume, trial, stimType]
            offsetResp(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, offsetRespVols,:);                                 % --> [x, y, plane, volume, trial, stimType]

            % Calculate dF/F for trial-averaged data
            disp('Calculating trial-averaged dF/F...')
            onsetDff(:,:,:,:,iStim) = calc_dFF(onsetResp(:,:,:,:,:,iStim), onsetBaseline(:,:,:,:,:,iStim), 5);          % --> [x, y, plane, volume, stimType]
            offsetDff(:,:,:,:,iStim) = calc_dFF(offsetResp(:,:,:,:,:,iStim), offsetBaseline(:,:,:,:,:,iStim), 5);       % --> [x, y, plane, volume, stimType]

            % Avgerage dF/F data across volumes as well
            disp('Calculating trial- and volume-averaged dF/F...')
            onsetDffAvg(:,:,:,iStim) = calc_dFF(onsetResp(:,:,:,:,:,iStim), onsetBaseline(:,:,:,:,:,iStim), [4 5]);     % --> [x, y, plane, stimType]
            offsetDffAvg(:,:,:,iStim) = calc_dFF(offsetResp(:,:,:,:,:,iStim), offsetBaseline(:,:,:,:,:,iStim), [4 5]);  % --> [x, y, plane, stimType]
            disp('Calculations complete');
        end
    end%iFold

                    %% PLOT WIND STIM ONSET RESPONSE HEATMAPS FOR EACH PLANE

    smoothingSigma = [0.5];
    scalingFactor = [0.6];
    makeVid = 0;
    saveDir = [];
    fileName = 'Wind_Onset_Plane_Heatmaps';   

    % Calculate absolute max dF/F value across all planes and stim types
    range = calc_range(onsetDffAvg, scalingFactor);

    % Plot figures
    [~, ~] = plot_heatmaps(onsetDffAvg, myData, range, stimTypes, [], [], smoothingSigma, makeVid, saveDir, fileName);

                    %% PLOT WIND STIM OFFSET RESPONSE HEATMAPS FOR EACH PLANE

    smoothingSigma = [0.5];
    scalingFactor = [1];
    makeVid = 0;
    saveDir = [];
    fileName = 'Wind_Offset_Plane_Heatmaps';       

    % Calculate absolute max dF/F value across all planes and stim types
    range = calc_range(offsetDffAvg, scalingFactor);

    % Plot figures
    [~, ~] = plot_heatmaps(offsetDffAvg, myData, range, stimTypes, [], [], smoothingSigma, makeVid, saveDir, fileName);

                    %% CREATE VIDEO OF MEAN dF/F THROUGHOUT WIND RESPONSES

    fileName = ['Wind_Stim_Responses'];
    sid = 1;
    smoothingSigma = 0.5;
    scalingFactor = 0.5;
    offsetAlign = 0;

    % Make sure wind stim types were combined for dF/F calculation
    if combineStimTypes

        % Select correct alignment
        if offsetAlign
            windDffVols = offsetDff(:,:,:,:, 1);
            baselineVols = offsetBaselineVols;
            respVols = offsetRespVols;
            titleStr = 'offset';
        else
            windDffVols = onsetDff(:,:,:,:, 1)
            baselineVols = onsetBaselineVols;
            respVols = onsetRespVols;
            titleStr = 'onset';
        end


        % Calculate absolute max dF/F value across all planes and action states
        range = calc_range(windDffVols,scalingFactor);

        % Calculate volume times in seconds relative to wind onset
        baselineVolTimes = -(1:length(baselineVols)) / myData.volumeRate;
        stimVolTimes = 1:length(respVols) / myData.volumeRate;
        relTimes = [baselineVolTimes(end:-1:1), stimVolTimes];

        % Create cell array with titles for each frame
        titleStrings = [];
        for iVol = 1:size(windDffVols, 4)
            if iVol <= size(baselineVols, 4)
                titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before wind ', titleStr];
            else
                titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after wind ', titleStr];
            end
        end

        % Create video
        savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
        make_heatmap_vid(windDffVols, myData, range, fileName, titleStrings, savePath, [], [], [], smoothingSigma);
    else
        disp('Error: wind stim types must be combined to create this video')
    end
end%iFold
    
%% =================================================================================================
%           PURE BEHAVIOR ANALYSES                                   
%%%=================================================================================================
for iFold = 1
    %% CALCULATE MEAN dF/F ACROSS BEHAVIORAL STATES

    locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
    actionLabel = 2;

for iFoldIn = 1
    % Identify behavioral state during each volume
    actionVols = zeros(myData.nTrials, myData.nVolumes); stoppedVols = actionVols;
    for iTrial = 1:myData.nTrials
        if myData.goodTrials(iTrial)

            % Pull out action numbers for each volume
            currActions = myData.trialAnnotations{iTrial}.actionNums;
            volActions = currActions(volFrames);

            % Identify volume actions
            actionVols(iTrial, :) = (volActions == actionLabel);    %--> [trial, vol]
            stoppedVols(iTrial, :) = (volActions == noActionLabel); %--> [trial, vol]
        else
            % So data from invalid trials won't ever be matched to an action state
            actionVols(iTrial, :) = 0;
            stoppedVols(iTrial, :) = 0;
        end
    end

    % Calculate average values for each plane across behavioral states
    meanActionVols = [];
    meanStoppedVols = [];
    test = [];
    imgData = myData.wholeSession; %--> [x, y, plane, volume, trial]
    for iTrial = 1:myData.nTrials
       currImgData = imgData(:,:,:,:,iTrial);
       currActionVols = logical(actionVols(iTrial,:));
       currStoppedVols = logical(stoppedVols(iTrial,:));
    test = 1;
       % Exclude volumes that occurred during or just after the wind stim
       if myData.stimSepTrials.windTrials(iTrial)
           stimStartVol = ceil(stimStart * volumeRate);
           stimEndVol = floor(stimEnd * volumeRate);
           excludeVols = stimStartVol:(stimEndVol + ceil(volumeRate));
           currImgData(:,:,:,excludeVols) = [];
           currActionVols(excludeVols) = [];
           currStoppedVols(excludeVols) = [];
       end
    test = 2;
       % Pull out running volumes, if any exist, from the current trial
       if sum(currActionVols) > 0
           meanActionVols(:,:,:,end+1) = mean(currImgData(:,:,:,currActionVols),4);    %--> [x, y, plane, trial]
       end 
    test = 3;
       % Pull out stopping volumes, if any exist, from the current trial
       if sum(currStoppedVols) > 0
           meanStoppedVols(:,:,:,end+1) = mean(currImgData(:,:,:,currStoppedVols),4);  %--> [x, y, plane, trial]
       end

    end
    actionMean = mean(meanActionVols, 4);   %--> [x, y, plane]
    stoppedMean = mean(meanStoppedVols, 4); %--> [x, y, plane] 

    % Get dF/F values for action relative to quiescence
    actionDff = (actionMean - stoppedMean) ./ stoppedMean; % --> [x, y, plane]
    
end%iFoldIn
                    %% PLOT BEHAVIORAL STATE HEATMAPS FOR EACH PLANE

    smoothingSigma = [0.5];    
    scalingFactor = [0.5];
    makeVid = 1;
    saveDir = [];
    fileName = 'Locomotion_Plane_Heatmaps';
    titleStr = {'dF/F - Locomotion vs. Quiescence'}

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(actionDff, scalingFactor);

    % Plot figures
    [~, ~] = plot_heatmaps(actionDff, myData, range, titleStr, [], [], smoothingSigma, makeVid, saveDir, fileName);


    %% CALCULATE MEAN dF/F AROUND LOCOMOTION ONSET

    %---------- Identify behavioral state during each volume ----------
    runVols = zeros(myData.nTrials, myData.nVolumes); 
    stoppedVols = runVols; legMoveVols = runVols; volActions = runVols;
    onsetVols = [];
    for iTrial = 1:myData.nTrials

        if myData.goodTrials(iTrial)

            % Pull out action numbers for each volume
            currActions = myData.trialAnnotations{iTrial}.actionNums;
            volActions(iTrial,:) = currActions(volFrames);

            % Identify volume actions
            locomotionLabel = 2; noActionLabel = 0; isoMovementLabel = 4;
            runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel);       % [trial, vol]
            stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel);     % [trial, vol]
            legMoveVols(iTrial, :) = (volActions(iTrial,:) == isoMovementLabel);  % [trial, vol]

            % Find onsets of running bouts >respLen volumes in duration and preceded by >baseLen volumes of quiescence
            baseLen = 12;
            respLen = 24;
            actionPattern = [zeros(1,baseLen), (ones(1,respLen) * locomotionLabel)];
            patternLen = length(actionPattern);
            currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);

            % Discard any bouts that occurred during or just after the wind stim
            if myData.stimSepTrials.windTrials(iTrial)
                stimStartVol = ceil(stimStart * volumeRate);
                stimEndVol = floor(stimEnd * volumeRate);
                currTrialOnsets(currTrialOnsets > (stimStartVol - patternLen) & currTrialOnsets < (stimEndVol + ceil(volumeRate))) = [];
            end

            onsetVols{iTrial} = currTrialOnsets;

        else
            % So data from invalid trials won't ever be matched to an action state
            runVols(iTrial, :) = 0;
            stoppedVols(iTrial, :) = 0;
            legMoveVols(iTrial, :) = 0;
            volActions(iTrial, :) = -1000;
            onsetVols{iTrial} = [];
        end
    end

    %---------- Get imaging data for running onsets ----------
    onsetData = [];
    for iTrial = 1:myData.nTrials
        % myData.wholeSession = [x, y, plane, volume, trial]
        if ~isempty(onsetVols{iTrial})
            onsets = onsetVols{iTrial};
            for iOnset = 1:length(onsets)
                volIdx = onsets(iOnset):onsets(iOnset) + patternLen-1;
                onsetData(:,:,:,:,end+1) =  myData.wholeSession(:,:,:, volIdx, iTrial); % --> [x, y, plane, onsetVolume, onsetNum]
            end
        end
    end

    % Calculate dF/F before and after movment onset using pre-movement period as baseline
    onsetBaselines = onsetData(:,:,:, 1:baseLen,:);                            % --> [x, y, plane, volume, onsetNum]
    onsetBaselineMean = mean(mean(onsetBaselines, 5), 4);                      % --> [x, y, plane]
    onsetBaselineMeanRep = repmat(onsetBaselineMean, 1, 1, 1, patternLen);     % --> [x, y, plane, volume]
    onsetMean = mean(onsetData, 5);                                            % --> [x, y, plane, volume]
    onsetDffVols = (onsetMean - onsetBaselineMeanRep) ./ onsetBaselineMeanRep; % --> [x, y, plane, volume]
    onsetMeanDff = mean(onsetDffVols, 4);                                      % --> [x, y, plane]

                    %% PLOT MOVEMENT ONSET HEATMAPS FOR EACH PLANE

    smoothingSigma = [0.5];    
    scalingFactor = [0.5];
    makeVid = 1;
    saveDir = [];
    fileName = 'Behavior_Onset_Plane_Heatmaps';

    % Calculate dF/F range
    range = calc_range(onsetMeanDff, scalingFactor);

    % Plot figures
    [~, ~] = plot_heatmaps(onsetMeanDff, myData, range, {'dF/F - Locomotion onset'}, [], [], smoothingSigma, makeVid, saveDir, fileName);


                    %% CREATE VIDEO OF MEAN dF/F FOR EACH PLANE THROUGHOUT MOVEMENT ONSET

    smoothingSigma = [0.5];    
    rangeScalar = 0.4;

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(onsetDffVols, rangeScalar);


    % Calculate volume times in seconds relative to movement onset
    baselineVolTimes = -(1:size(onsetBaselines, 4))/myData.volumeRate;
    stimVolTimes = ((1:(size(onsetData, 4)-size(onsetBaselines,4)))/myData.volumeRate);
    relTimes = [baselineVolTimes(end:-1:1), stimVolTimes];

    % Generate titles for each volume
    titleStrings = {};
    for iVol = 1:size(onsetDffVols, 4)
        if iVol <= size(onsetBaselines, 4);
            titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before movement onset'];
        else
            titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after movement onset'];
        end
    end

    % Create video
    savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
    fileName = 'Move_Onset_Responses';
    make_heatmap_vid(onsetDffVols, myData, range, fileName, titleStrings, savePath, [], [], [], smoothingSigma);


    %% CREATE COMBINED PLOTTING VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

    baselineLength = 12; % Number of baseline volumes to analyze/plot for each movement bout

for iFoldIn = 1
    % Identify starting and ending indices of all movement bouts
    volActionCell = num2cell(volActions, 2);
    for iTrial = 1:myData.nTrials
        volActionCell{iTrial} = num2str(volActionCell{iTrial}, '%u'); % Providing FormatSpec so it doesn't add spaces
    end
    baselineStr = num2str(zeros(1,baselineLength), '%u');
    [startInd, endInd] = regexp(volActionCell, [baselineStr, '2+']);

    % Create array with information needed to identify bout data
    allBouts = [];
    for iTrial = 1:myData.nTrials
        if ~isempty(startInd{iTrial})
            for iBout = 1:length(startInd{iTrial})
              allBouts(end+1,:) = [iTrial, startInd{iTrial}(iBout), endInd{iTrial}(iBout)]; % --> [trialNum, startInd, endInd]
            end
        end
    end

    % Create save directory
    savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', myData.expDate, '\sid_0\Analysis'];
    fileName = 'Combined_Plotting_All_Movement_Bouts';
    if ~isdir(savePath)
        mkdir(savePath);
    end

    % Warn user and offer to cancel save if this will overwrite an existing file
    overwrite = 1;
    if exist(fullfile(savePath, [fileName, '.avi']), 'file') ~= 0
        dlgAns = questdlg('Creating this video will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end

    % Create video writer
    myVid = VideoWriter(fullfile(savePath, fileName));
    myVid.FrameRate = 25;
    open(myVid)

    for iBout = 1:size(allBouts,1)

        disp(['Plotting bout ', num2str(iBout), ' of ' num2str(size(allBouts, 1))]);

        % Calculate dF/F for the current bout ****Is this an acceptable way to calculate dF/F?*****
        boutData = double(squeeze(myData.wholeSession(:,:,:, allBouts(iBout,2):allBouts(iBout,3),allBouts(iBout,1))));   % --> [x, y, plane, volume]
        boutBaseline = mean(boutData(:,:,:,1:baselineLength),4);                                                         % --> [x, y, plane]
        boutBaselineRep = repmat(boutBaseline,1,1,1,length(allBouts(iBout,2):allBouts(iBout,3)));                        % --> [x, y, plane, volume]
        boutDff = (boutData - boutBaselineRep) ./ boutBaselineRep;                                                       % --> [x, y, plane, volume]

        boutDff(isinf(boutDff)) = 0; % To eliminate inf values from dividiing by zero above...baseline shouldn't be zero in valid data anyways
        boutDff(isnan(boutDff)) = 0;

        % Average each volume with its neighbors to improve SNR
        boutDff = movmean(boutDff, 3, 3);

        % Calculate dF/F value range
        range = calc_range(boutDff,[0.8]);

        % Load behavior vid for the current trial
        vidDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', myData.expDate, '_Movies');
        trialStr = ['sid_', num2str(myData.sid), '_tid_', sprintf('%03.f',(allBouts(iBout, 1)))];
        myMovie = [];
        flyVid = VideoReader(fullfile(vidDir, [trialStr, '.mp4']));
        while hasFrame(flyVid)
            currFrame = readFrame(flyVid);
            myMovie(:,:,end+1) = rgb2gray(currFrame);
        end
        myMovie = uint8(myMovie(:,:,2:end)); % Adds a black first frame for some reason, so drop that

        % Calculate frame times in seconds relative to behavior bout onset
        boutFrames = find(volFrames == allBouts(iBout,2),1):find(volFrames == allBouts(iBout,3),1);
        baselineLengthFrames = floor(find(volFrames == baselineLength, 1, 'last'));
        baselineFrameTimes = frameTimes(baselineLengthFrames:-1:1);
        boutFrameTimes = frameTimes(1:(numel(boutFrames) - baselineLengthFrames));
        relFrameTimes = [baselineFrameTimes; boutFrameTimes];
        for iFrame = 1:numel(boutFrames)
            if iFrame <= baselineLengthFrames
                relTimeStr = 'before';
            else
                relTimeStr = 'after';
            end
            titleStrings{iFrame} = ['Time = ', sprintf('%05.2f', relFrameTimes(iFrame)), ' sec ', relTimeStr, ' locomotion onset'];
        end


        for iFrame = 1:numel(boutFrames)

            % Create fig
            f = figure(1); clf
            f.Position = [50 45, 1620, 855];

            for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral

                % Plot dF/F for each plane
                ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);

                if iPlane == 1
                    imshow(myMovie(:,:,boutFrames(iFrame)));
                    axis image; axis off
                else
                    imagesc(boutDff(:, :, iPlane, volFrames(iFrame)))
                    caxis(range)
                    colormap(ax, bluewhitered) %bluewhitered % 'parula'
                    axis equal; axis off
                end

                % Label positions
                if iPlane == myData.nPlanes
                    title('Ventral')
                elseif iPlane == 2
                    title('Dorsal')
                end
            end

            suptitle(['Trial #', num2str(allBouts(iBout,1)), ', Movement bout #', num2str(iBout), ' of ', ...
                num2str(size(allBouts,1)), '  -  ', titleStrings{iFrame}]);

            % Write frame to video
            writeFrame = getframe(f);
            writeVideo(myVid, writeFrame);
        end
    end%for
    close(myVid)
end%iFoldIn
end%iFold

%% =================================================================================================
%           BALL STOPPING ANALYSES                                   
%%%=================================================================================================
for iFold = 1
    %% CALCULATE MEAN dF/F AROUND BALL STOPPING

    % Extract ball-stopping event timing for each trial
    ballStopVols = [];
    for iTrial = 1:nTrials
        currAnnot = num2str(myData.trialAnnotations{iTrial}.ballStopNums');
        currAnnot = currAnnot(~isspace(currAnnot));
        moveStart = regexp(currAnnot, '[02]1') + 1;
        ballStop = regexp(currAnnot, '[10]2') + 1;
        ballRelease = regexp(currAnnot, '2[10]') + 1;
        moveEnd = regexp(currAnnot, '10000') + 1;
        ballStopFrames = [moveStart(1), ballStop, ballRelease, moveEnd];                   % --> [startMoveUp, ballStop, ballRelease, endMoveDown]
        for iTime = 1:length(ballStopFrames) 
           [~, ballStopVols(iTrial, iTime)] = min(abs(volFrames - ballStopFrames(iTime))); % --> [startMoveUp, ballStop, ballRelease, endMoveDown]
        end
    end

    % Find minimum duration of each annotation period
    moveUpDur = min(ballStopVols(:,2) - ballStopVols(:,1));
    ballStopDur = min(ballStopVols(:,3) - ballStopVols(:,2));
    moveDownDur = min(ballStopVols(:,4) - ballStopVols(:,3));

    % Extract imaging data for each phase
    wsSize = size(myData.wholeSession);
    moveUpBaselineData = zeros([wsSize(1:3), moveUpDur, nTrials]);
    ballStopBaselineData = zeros([wsSize(1:3), ballStopDur, nTrials]);
    moveDownBaselineData = zeros([wsSize(1:3), moveDownDur, nTrials]);
    moveUpData = moveUpBaselineData;
    ballStopData = ballStopBaselineData;
    moveDownData = moveDownBaselineData;

    for iTrial = 1:nTrials
        moveUpBaselineData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, (ballStopVols(iTrial, 1) - moveUpDur):(ballStopVols(iTrial, 1) - 1), iTrial);
        moveUpData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, ballStopVols(iTrial, 1):(ballStopVols(iTrial, 1) + moveUpDur - 1), iTrial);
        ballStopBaselineData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, (ballStopVols(iTrial, 2) - ballStopDur):(ballStopVols(iTrial, 2) - 1), iTrial);
        ballStopData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, ballStopVols(iTrial, 2):(ballStopVols(iTrial, 2) + ballStopDur - 1), iTrial);
        moveDownBaselineData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, (ballStopVols(iTrial, 3) - moveDownDur):(ballStopVols(iTrial, 3) - 1), iTrial);
        moveDownData(:,:,:,:,iTrial) = myData.wholeSession(:,:,:, ballStopVols(iTrial, 3):(ballStopVols(iTrial, 3) + moveDownDur - 1), iTrial);
    end

    % 8:13  15:20 22:26  29:33  35:39
    % [1 4 14 21 27 34 50] 

    earlyTrials =   1:(nTrials/2);
    lateTrials =   ((nTrials/2)+1):nTrials;

    % Calc dF/F for trial-averaged data
    moveUpDff = calc_dFF(moveUpData, moveUpBaselineData, 5);                    % --> [x, y, plane, volume]
    ballStopDff = calc_dFF(ballStopData, ballStopBaselineData, 5);              % --> [x, y, plane, volume]
    moveDownDff = calc_dFF(moveDownData, moveDownBaselineData, 5);              % --> [x, y, plane, volume]
    moveDownRelDff = calc_dFF(moveDownData, ballStopBaselineData, 5);           % --> [x, y, plane, volume]
    ballStopEarlyTrialsDff = calc_dFF(ballStopData(:,:,:,:, earlyTrials), ballStopBaselineData(:,:,:,:, earlyTrials), 5);   % --> [x, y, plane, volume]
    moveDownEarlyTrialsDff = calc_dFF(moveDownData(:,:,:,:, earlyTrials), moveDownBaselineData(:,:,:,:, earlyTrials), 5);   % --> [x, y, plane, volume]
    ballStopLateTrialsDff = calc_dFF(ballStopData(:,:,:,:, lateTrials), ballStopBaselineData(:,:,:,:, lateTrials), 5);      % --> [x, y, plane, volume]
    moveDownLateTrialsDff = calc_dFF(moveDownData(:,:,:,:, lateTrials), moveDownBaselineData(:,:,:,:, lateTrials), 5);      % --> [x, y, plane, volume] 

    % Average the dF/F data over volumes as well
    moveUpDffAvg = calc_dFF(moveUpData, moveUpBaselineData, [4 5]);             % --> [x, y, plane]
    ballStopDffAvg = calc_dFF(ballStopData, ballStopBaselineData, [4 5]);       % --> [x, y, plane]
    moveDownDffAvg = calc_dFF(moveDownData, moveDownBaselineData, [4 5]);       % --> [x, y, plane]
    moveDownRelDffAvg = calc_dFF(moveDownData, ballStopBaselineData, [4 5]);    % --> [x, y, plane]
    ballStopEarlyTrialsDffAvg = calc_dFF(ballStopData(:,:,:,:, earlyTrials), ballStopBaselineData(:,:,:,:, earlyTrials), [4 5]);    % --> [x, y, plane]
    moveDownEarlyTrialsDffAvg = calc_dFF(moveDownData(:,:,:,:, earlyTrials), moveDownBaselineData(:,:,:,:, earlyTrials), [4 5]);    % --> [x, y, plane]
    ballStopLateTrialsDffAvg = calc_dFF(ballStopData(:,:,:,:, lateTrials), ballStopBaselineData(:,:,:,:, lateTrials), [4 5]);       % --> [x, y, plane]
    moveDownLateTrialsDffAvg = calc_dFF(moveDownData(:,:,:,:, lateTrials), moveDownBaselineData(:,:,:,:, lateTrials), [4 5]);       % --> [x, y, plane]

    
                    %% PLOT HEATMAPS ALIGNED TO BALL STOPPING

    currPhases = [6 8 7 9];
    smoothingSigma = [0.5];    
    scalingFactor = [0.05];
    makeVid = 1;
    saveDir = [];
    fileName = 'ballStoppingHeatmaps_EarlyLate';

for iFoldIn = 1
    phasesDff = zeros([size(moveUpDffAvg), 3]);
    phasesDff(:,:,:,1) = moveUpDffAvg;                      % --> [x, y, plane, trialPhase]
    phasesDff(:,:,:,2) = ballStopDffAvg;                    % --> [x, y, plane, trialPhase]
    phasesDff(:,:,:,3) = moveDownDffAvg;                    % --> [x, y, plane, trialPhase]
    phasesDff(:,:,:,4) = moveDownRelDffAvg;                 % --> [x, y, plane, trialPhase]
    phasesDff(:,:,:,5) = moveDownRelDffAvg - ballStopDffAvg;
    phasesDff(:,:,:,6) = ballStopEarlyTrialsDffAvg;
    phasesDff(:,:,:,7) = moveDownEarlyTrialsDffAvg;
    phasesDff(:,:,:,8) = ballStopLateTrialsDffAvg;
    phasesDff(:,:,:,9) = moveDownLateTrialsDffAvg;


    titleStrings = [];
    titleStrings{1} = {'Actuator move up onset'};
    titleStrings{2} = {'Ball stopping'};
    titleStrings{3} = {'Ball release'};
    titleStrings{4} = {'Ball release (pre-stop baseline)'};
    titleStrings{5} = {'Release (pre-stop baseline) - stopping'};
    titleStrings{6} = {'Ball stopping - early trials'};
    titleStrings{7} = {'Ball release - early trials'};
    titleStrings{8} = {'Ball stopping - late trials'};
    titleStrings{9} = {'Ball release - late trials'};

    % Calculate dF/F range
    range = calc_range(phasesDff(:,:,:, currPhases), scalingFactor);

    % Plot figures
    [~, ~] = plot_heatmaps(phasesDff(:,:,:, currPhases), myData, range, titleStrings(currPhases), [], [], smoothingSigma, makeVid, saveDir, fileName);

end%iFoldIn
    %% SEPARATE TRIALS BASED ON FLY'S BEHAVIOR AROUND BALL STOPPING

    %---------- Identify behavioral state during each volume ----------

    preStopTime = 1;
    postStopTime = 2;
    preReleaseTime = 1;
    postReleaseTime = 2;

for iFoldIn = 1
    volActions = zeros(myData.nTrials, myData.nVolumes); % --> [trial, volume]
    ballStopActions = volActions;
    preStopMove = zeros(myData.nTrials, 1);
    postStopMove = preStopMove; preReleaseMove = preStopMove; postReleaseMove = preStopMove;
    for iTrial = 1:myData.nTrials

        if myData.goodTrials(iTrial)

            % Pull out action numbers and ball stopping data for each volume
            currActions = myData.trialAnnotations{iTrial}.actionNums;
            currBallStop = myData.trialAnnotations{iTrial}.ballStopNums;
            volActions(iTrial,:) = currActions(volFrames);
            ballStopActions(iTrial, :) = currBallStop(volFrames);
            preStopVols = floor(ballStopVols(iTrial, 2) - (preStopTime * volumeRate):(ballStopVols(iTrial, 2) - 1));
            postStopVols = ballStopVols(iTrial, 2):floor(ballStopVols(iTrial, 2) + (postStopTime * volumeRate));
            preReleaseVols = floor(ballStopVols(iTrial, 3) - (preReleaseTime * volumeRate):(ballStopVols(iTrial, 3) - 1));
            postReleaseVols = ballStopVols(iTrial, 3):floor(ballStopVols(iTrial, 3) + (postReleaseTime * volumeRate));       

            % Determine if fly was active during each period
            if sum(volActions(iTrial, preStopVols)) > 0
               preStopMove(iTrial) = 1; 
            end
            if sum(volActions(iTrial, postStopVols)) > 0
               postStopMove(iTrial) = 1; 
            end
            if sum(volActions(iTrial, preReleaseVols)) > 0
               preReleaseMove(iTrial) = 1; 
            end
            if sum(volActions(iTrial, postReleaseVols)) > 0
               postReleaseMove(iTrial) = 1; 
            end        
        else
            % So data from invalid trials won't ever be matched to an action state
            volActions(iTrial, :) = -1000;
        end        
    end

    % Separate trials
    tc = [preStopMove, postStopMove, preReleaseMove, postReleaseMove]; 
    ballStopNoMove          = tc(:,1) == 0 & tc(:,2) == 0; % --------------------------%
    ballStopStartMove       = tc(:,1) == 0 & tc(:,2) == 1; % dF/F aligned to ball stop %         
    ballStopEndMove         = tc(:,1) == 1 & tc(:,2) == 0; %                           %
    ballStopContMove        = tc(:,1) == 1 & tc(:,2) == 1; %---------------------------%
    ballReleaseNoMove       =                               tc(:,3) == 0 & tc(:,4) == 0;  %------------------------------% 
    ballReleaseStartMove    =                               tc(:,3) == 0 & tc(:,4) == 1;  % dF/F aligned to ball release %            
    ballReleaseEndMove      =                               tc(:,3) == 1 & tc(:,4) == 0;  %                              %
    ballReleaseContMove     =                               tc(:,3) == 1 & tc(:,4) == 1;  %------------------------------%            
    trialConds = [ballStopNoMove,ballStopStartMove,ballStopEndMove,ballStopContMove,ballReleaseNoMove,ballReleaseStartMove,ballReleaseEndMove,ballReleaseContMove];
    sum(trialConds)
    trialCondNames = {'ballStopNoMove','ballStopStartMove','ballStopEndMove', 'ballStopContMove','ballReleaseNoMove','ballReleaseStartMove','ballReleaseEndMove','ballReleaseContMove'}

    %%% CALCULATE MEAN dF/F FOR EACH TRIAL CONDITION
    ws = size(myData.wholeSession);
    wholeTrialAvg = zeros([ws(1:4), length(trialCondNames)]);
    preStopAvg = zeros([ws(1:3), length(trialCondNames)]);
    postStopAvg = preStopAvg; preReleaseAvg = preStopAvg; postReleaseAvg = preStopAvg; 
    for iCond = 1:length(trialCondNames)

        %---------- Get trial averaged baseline and stimulus data ----------                        % myData.wholeSession = [x, y, plane, volume, trial]                                                                                                          
        disp(['Trial condition ' num2str(iCond)]);
        wholeTrialAvg(:,:,:,:,iCond) = mean(myData.wholeSession(:,:,:,:,trialConds(:,iCond)), 5);   % --> [x, y, plane, volume, trialCondition]
        preStopAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,preStopVols, iCond), 4);                 % --> [x, y, plane, trialCondition]          
        postStopAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,postStopVols, iCond), 4);               % --> [x, y, plane, trialCondition]
        preReleaseAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,preReleaseVols, iCond), 4);           % --> [x, y, plane, trialCondition]
        postReleaseAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,postReleaseVols, iCond), 4);         % --> [x, y, plane, trialCondition]

    end%for

    % Calculate dF/F values
    dffAvgBallStop = (postStopAvg - preStopAvg) ./ preStopAvg;                                      % --> [x, y, plane, trialCondition]
    dffAvgBallRelease = (postReleaseAvg - preReleaseAvg) ./ preReleaseAvg;                          % --> [x, y, plane, trialCondition]

    % Want the baseline/response aligned to ball stop for the first 4 trial conditions, ball release for the rest
    dffAvg = dffAvgBallStop;
    dffAvg(:,:,:, 5:8) = dffAvgBallRelease(:,:,:, 5:8);                                             % --> [x, y, plane, trialCondition]

    % Eliminate any inf values caused by dividing by zero above
    dffAvg(isinf(dffAvg)) = 0;                                                                      % --> [x, y, plane, trialCondition]


    ballStopMoveDiff = dffAvgBallStop(:,:,:,2) - dffAvgBallStop(:,:,:,1);                           % --> [x, y, plane] 
    ballReleaseMoveDiff = dffAvgBallStop(:,:,:,6) - dffAvgBallStop(:,:,:,5);                        

end%iFoldIn
                    %% PLOT BALL STOP/RELEASE HEATMAPS FOR SOME TRIAL CONDITIONS

    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [1 2 5 6];
    smoothingSigma = [0.5];    
    scalingFactor = [0.25];
    makeVid = 1;
    saveDir = [];
    fileName = 'Ball_Stop_Behavior_Interaction_Heatmaps';

    dffConds = reshape(dffAvg(:,:,:,currConds), [1 numel(dffAvg(:,:,:,currConds))]);
    range = calc_range(dffConds, scalingFactor);

    % Plot figures
    dffCurrConds = dffAvg(:,:,:,currConds); %cat(4, ballStopMoveDiff, ballReleaseMoveDiff);
    [f, plotAxes] = plot_heatmaps(dffCurrConds, myData, range, trialCondNames(currConds), [], [], smoothingSigma, makeVid, saveDir, fileName); %trialCondNames(currConds)

end%iFold

%% =================================================================================================
%           OTHER ANALYSES                                   
%%%=================================================================================================
for iFold = 1

%% CREATE VIDEO OF MEAN dF/F THROUGHOUT ENTIRE TRIAL

%++++++++++ Calculate whole-trial dF/F using overall mean as baseline ++++++++++

allWindTrials = myData.wholeSession(:,:,:,:,myData.stimSepTrials.windTrials);   % --> [x, y, plane, volume, trial]   

% Calculate dF/F using whole trial average as baseline
baseline = mean(mean(allWindTrials, 5), 4);                           % --> [x, y, plane]
baselineMeanRep = repmat(baseline, 1, 1, 1, size(allWindTrials, 4));  % --> [x, y, plane, volume]
trialAvg = mean(allWindTrials, 5);                                    % --> [x, y, plane, volume]
wholeTrialDff = (trialAvg - baselineMeanRep) ./ baselineMeanRep;      % --> [x, y, plane, volume]

% Calculate absolute max dF/F value across all planes and action states
rangeScalar = .3;
range = calc_range(wholeTrialDff, rangeScalar);

smoothingSigma = [0.5];

%----------Create video of dF/F for each plane throughout trial----------

% Add title above all subplots
titleStrings = [];
for iVol = 1:size(wholeTrialDff, 4)
    if iVol <= stimStart * volumeRate
        titleStr = 'Before wind onset';
    elseif iVol <= stimEnd * volumeRate
        titleStr = 'During wind stimulus';
    else
        titleStr = 'After wind stimulus';
    end
    titleStrings{iVol} = [titleStr, '  -  time = ', num2str(round(iVol./volumeRate, 1))];
end

% Create video
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', myData.expDate, '\sid_0\Analysis'];
fileName = 'Full_Trial_dFF';
make_heatmap_vid(wholeTrialDff, myData, range, fileName, titleStrings, savePath, [], [], [], smoothingSigma);

%----------Create video of mean raw fluorescence signal throughout trial----------

% Calculate range
rangeScalar = 0.25;
range = calc_range(trialAvg, rangeScalar);

% Add title above all subplots
titleStrings = [];
for iVol = 1:size(wholeTrialDff, 4)
    if iVol <= stimStart * volumeRate
        titleStr = 'Before wind onset';
    elseif iVol <= stimEnd * volumeRate
        titleStr = 'During wind stimulus';
    else
        titleStr = 'After wind stimulus';
    end
    titleStrings{iVol} = [titleStr, '  -  time = ', num2str(round(iVol./volumeRate, 1))];
end

% Update file name and create video
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', myData.expDate, '\sid_0\Analysis'];
fileName = 'Full_Trial_RawF';
make_heatmap_vid(trialAvg, myData, range, fileName, titleStrings, savePath, [], [], [], smoothingSigma);


%% SEPARATE TRIALS BASED ON WHETHER FLY WAS MOVING DURING WIND STIM

%---------- Identify behavioral state during each volume ----------
volActions = zeros(myData.nTrials, myData.nVolumes);
preStimMove = zeros(myData.nTrials, 1); stimMove = preStimMove; postStimMove = preStimMove;
for iTrial = 1:myData.nTrials
    
    if myData.goodTrials(iTrial)
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions(iTrial,:) = currActions(volFrames);
        
        preStimTime = 1;
        stimTime = 1;
        postStimTime = 1;
        preStimVols = floor((stimStart-preStimTime)*volumeRate):floor(stimStart*volumeRate);
        stimVols = floor(stimStart*volumeRate):floor((stimStart+stimTime)*volumeRate);
        postStimVols = ceil(stimEnd*volumeRate):floor((stimEnd+postStimTime)*volumeRate);
        
        % Determine if fly was active just before the wind onset
        if sum(volActions(iTrial, preStimVols)) > 0
            preStimMove(iTrial) = 1;
        end
        
        % Determine if fly was moving in the first few seconds after wind onset
        if sum(volActions(iTrial, stimVols)) > 0
            stimMove(iTrial) = 1;
        end
        
        % Determine if fly was moving in the first few seconds after wind offset
        if sum(volActions(iTrial, postStimVols)) > 0
            postStimMove(iTrial) = 1;
        end
        
    else
        % So data from invalid trials won't ever be matched to an action state
        volActions(iTrial, :) = -1000;
    end        
end

% Separate trials 
tc = [preStimMove, stimMove, postStimMove, myData.stimSepTrials.windTrials'];
windOnsetNoMove         = tc(:,1) == 0 & tc(:,2) == 0                & tc(:,4) == 1;
windOffsetNoMove        = tc(:,1) == 0 & tc(:,2) == 0 & tc(:,3) == 0 & tc(:,4) == 1; % dF/F calculation should be centered on wind offset for this condition
windOnsetStartMove      = tc(:,1) == 0 & tc(:,2) == 1                & tc(:,4) == 1;
windOffsetStartMove     = tc(:,1) == 0 & tc(:,2) == 0 & tc(:,3) == 1 & tc(:,4) == 1; % dF/F calculation should be centered on wind offset for this condition
windContMove            = tc(:,1) == 1 & tc(:,2) == 1                & tc(:,4) == 1;

noWindNoMove            = tc(:,1) == 0 & tc(:,2) == 0                & tc(:,4) == 0;
noWindStartMove         = tc(:,1) == 0 & tc(:,2) == 1                & tc(:,4) == 0;
noWindContMove          = tc(:,1) == 1 & tc(:,2) == 1                & tc(:,4) == 0;

trialConds = [windOnsetNoMove,windOffsetNoMove,windOnsetStartMove,windOffsetStartMove,windContMove,noWindNoMove,noWindStartMove,noWindContMove];
sum(trialConds)
trialCondNames = {'windOnsetNoMove','windOffsetNoMove','windOnsetStartMove', 'windOffsetStartMove','windContMove','noWindNoMove','noWindStartMove','noWindContMove'}

%%% CALCULATE MEAN dF/F FOR EACH TRIAL CONDITION

wholeTrialAvg = []; baselineAvg = []; dffAvg = []; stimAvg = [];
for iCond = 1:length(trialCondNames)
    
    %---------- Get trial averaged baseline and stimulus data ----------                                                % myData.wholeSession = [x, y, plane, volume, trial]                                                                                                          
    wholeTrialAvg(:,:,:,:,iCond) = mean(myData.wholeSession(:,:,:,:,trialConds(:,iCond)), 5);                                                 % --> [x, y, plane, volume, trialCondition]
        
    % Pre-stim
    if floor((stimStart-preStimTime)*volumeRate) > 0 
        baselineAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:, preStimVols, iCond), 4);                                   % --> [x, y, plane, StimType]
    else
        % if stimDuration > preStimDuration, start baseline one second after beginning of trial
        baselineAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),iCond), 4);  	% --> [x, y, plane, trialCondition]
    end
    
    % During stim
    stimAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:, stimVols, iCond), 4);                                              % --> [x, y, plane, trialCondition]
    
    % Post-stim 
    postStimAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:, postStimVols, iCond), 4);                                      % --> [x, y, plane, trialCondition]
    
end%for

% Calculate dF/F values
dffAvgOnset = (stimAvg - baselineAvg) ./ baselineAvg;                                                                   % --> [x, y, plane, trialCondition]
dffAvgOffset = (postStimAvg - stimAvg) ./ stimAvg;                                                                      % --> [x, y, plane, trialCondition]

% Want the baseline/response aligned to stim offset for two trial conditions, stim onset for the rest
dffAvg = dffAvgOnset;
dffAvg(:,:,:,[2 4]) = dffAvgOffset(:,:,:,[2 4]);                                                                        % --> [x, y, plane, trialCondition]

% Eliminate inf values from dividing by zero above...baseline shouldn't be zero in valid data anyways
dffAvg(isinf(dffAvg)) = 0;

                %% PLOT WIND STIM ONSET RESPONSE HEATMAPS FOR SOME TRIAL CONDITIONS
    
% Calculate absolute max dF/F value across all planes and stim types
currConds = [1 2 4];
smoothingSigma = [0.5];    
scalingFactor = [0.1];
makeVid = 0;
saveDir = [];
fileName = 'Wind_Behavior_Interaction_Heatmaps';

dffConds = reshape(dffAvg(:,:,:,currConds), [1 numel(dffAvg(:,:,:,currConds))]);
range = calc_range(dffConds, scalingFactor);

% Plot figures
dffCurrConds = dffAvg(:,:,:,currConds);
[f, plotAxes] = plot_heatmaps(dffCurrConds, myData, range, trialCondNames(currConds), [], [], smoothingSigma, makeVid, saveDir, fileName);


%% MAKE VIDEOS OF WIND AND MOVEMENT RESPONSES FOR JUST A SINGLE PLANE


baselineDur = 2;
responseDur = 4;
planeNum = 6;

% WIND RESPONSE DATA
stimSepTrials = [];

% Separate out center wind trials
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
end 
windTrials = myData.wholeSession(:,:,:,:,logical(stimSepTrials.CenterWind + stimSepTrials.RightWind));                                                   % --> [x, y, plane, volume, trial]

% Calculate dF/F before and after wind onset
stimLength = stimEnd - stimStart;
stimLengthVols = floor(responseDur * volumeRate);
baselineVols = windTrials(:,:,:,floor((stimStart-baselineDur)*volumeRate):floor(stimStart*volumeRate),:);               % --> [x, y, plane, volume, trial]

stimVols = windTrials(:,:,:,ceil((stimStart*volumeRate)-size(baselineVols, 4)):ceil((stimStart*volumeRate)+stimLengthVols),:); % --> [x, y, plane, volume, trial]  
stimVolsMean = mean(stimVols, 5);                                                                                                    % --> [x, y, plane, volume]
nVols = size(stimVols, 4);
baselineMean = mean(mean(baselineVols, 5), 4);                                                                                       % --> [x, y, plane]
baselineMeanRep = repmat(baselineMean, 1, 1, 1, nVols);                                                                              % --> [x, y, plane, volume]
windDffVols = (stimVolsMean - baselineMeanRep) ./ baselineMeanRep;                                                                   % --> [x, y, plane, volume]

% MOVEMENT RESPONSE DATA
%---------- Identify behavioral state during each volume ----------
runVols = zeros(myData.nTrials, myData.nVolumes); 
stoppedVols = runVols; legMoveVols = runVols; volActions = runVols;
onsetVols = [];
for iTrial = 1:myData.nTrials
    
    if myData.goodTrials(iTrial)
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions(iTrial,:) = currActions(volFrames);
        
        % Identify volume actions
        locomotionLabel = 2; noActionLabel = 0; isoMovementLabel = 4;
        runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel);   % [trial, vol]
        stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel); % [trial, vol]
        legMoveVols(iTrial, :) = (volActions(iTrial,:) == isoMovementLabel);  % [trial, vol]
        
        % Find onsets of running bouts
        baseLen = size(baselineVols, 4);
        respLen = size(stimVols,4) - size(baselineVols, 4);
        actionPattern = [zeros(1,baseLen), ones(1,respLen) * locomotionLabel];   % [0 0 0 0 0 0 0 2 2 2 2 2 2 2];
        patternLen = length(actionPattern);
        currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);
        
        % Discard any bouts that occurred during or just after the wind stim
        if myData.stimSepTrials.windTrials(iTrial)
            stimStartVol = ceil(stimStart * volumeRate);
            stimEndVol = floor(stimEnd * volumeRate);
            currTrialOnsets(currTrialOnsets > (stimStartVol - patternLen) & currTrialOnsets < (stimEndVol + ceil(volumeRate))) = [];
        end
        
        onsetVols{iTrial} = currTrialOnsets;

    else
        % So data from invalid trials won't ever be matched to an action state
        runVols(iTrial, :) = 0;
        stoppedVols(iTrial, :) = 0;
        legMoveVols(iTrial, :) = 0;
        volActions(iTrial, :) = -1000;
        onsetVols{iTrial} = [];
    end
end

%---------- Get imaging data for running onsets ----------
onsetData = [];
for iTrial = 1:myData.nTrials
    % myData.wholeSession = [x, y, plane, volume, trial]
    if ~isempty(onsetVols{iTrial})
        onsets = onsetVols{iTrial};
        for iOnset = 1:length(onsets)
            volIdx = onsets(iOnset):onsets(iOnset) + patternLen-1;
            onsetData(:,:,:,:,end+1) =  myData.wholeSession(:,:,:, volIdx, iTrial); % --> [x, y, plane, onsetVolume, onsetNum]
        end
    end
end

% Calculate dF/F before and after movment onset using pre-movement period as baseline
onsetBaselines = onsetData(:,:,:, 1:baseLen,:);                            % --> [x, y, plane, volume, onsetNum]
onsetBaselineMean = mean(mean(onsetBaselines, 5), 4);                      % --> [x, y, plane]
onsetBaselineMeanRep = repmat(onsetBaselineMean, 1, 1, 1, patternLen);     % --> [x, y, plane, volume]
onsetMean = mean(onsetData, 5);                                            % --> [x, y, plane, volume]
onsetDffVols = (onsetMean - onsetBaselineMeanRep) ./ onsetBaselineMeanRep; % --> [x, y, plane, volume]
onsetMeanDff = mean(onsetDffVols, 4);                                      % --> [x, y, plane]

                %% MAKE COMBINED VIDEO
    % Calculate dF/F ranges
    windRange = calc_range(windDffVols, []);
    moveRange = calc_range(onsetDffVols,[]);

    % Create save directory if necessary
    savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
    fileName = ['CombinedResponses_Plane_', num2str(planeNum)];
    if ~isdir(savePath)
        mkdir(savePath);
    end

    % Warn user and offer to cancel save if this will overwrite an existing file
     overwrite = 1;
    if exist(fullfile(savePath, [fileName, '.avi']), 'file') ~= 0
        dlgAns = questdlg('Creating this video will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end

    if overwrite

        % Calculate volume times in seconds relative to wind onset
        baselineVolTimes = -(1:size(baselineVols, 4))/myData.volumeRate;
        stimVolTimes = ((1:(size(stimVols, 4)-size(baselineVols,4)))/myData.volumeRate);
        relTimes = [baselineVolTimes(end:-1:1), stimVolTimes];

        % Create video writer
        myVid = VideoWriter(fullfile(savePath, fileName));
        myVid.FrameRate = 7;
        open(myVid)

        for iVol = 1:nVols

            % Create figure
            f = figure(1); clf
            f.Position = [50 45, 1620, 950];

            % Plot reference image for the current plane
            ax1 = subaxis(2,2,1, 'Spacing', 0.01, 'MB', 0.025);
            imshow(myData.refImg{planeNum}, [0 myData.MAX_INTENSITY])

            % Plot wind dF/F
            ax2 = subaxis(2, 2, 3, 'Spacing', 0.01, 'MB', 0.025);
            imagesc(windDffVols(:, :, planeNum, iVol))
            caxis(windRange)
            colormap(ax2, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off

            textAx = axes();
            textAx.Position = [0.14 0.42 0.2 0.1];
            textAx.Visible = 'Off';
            txt = text();
            txt.String = 'Wind Onset';
            txt.Units = 'Normalized';
            txt.Position = [0.6 0.12];
            txt.FontSize = 16;

            % Plot movement dF/F
            ax3 = subaxis(2, 2, 4, 'Spacing', 0.01, 'MB', 0.025);
            imagesc(onsetDffVols(:, :, planeNum, iVol))
            caxis(moveRange)
            colormap(ax3, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off

            textAx2 = axes();
            textAx2.Position = [0.565 0.40 0.2 0.1];
            textAx2.Visible = 'Off';
            txt2 = text();
            txt2.String = 'Movement Onset';
            txt2.Units = 'Normalized';
            txt2.Position = [0.45 0.3];
            txt2.FontSize = 16;

            % Display time relative to wind/movment onset
            titleTextAx = axes();
            titleTextAx.Position = [0.38 0.40 0.2 0.1];
            titleTextAx.Visible = 'Off';
            titleText = text();
            titleText.String = ['Time = ', sprintf('%04.2f', relTimes(iVol))];
            titleText.Units = 'Normalized';
            titleText.Position = [0.45 0.3];
            titleText.FontSize = 16;

            % Write frame to video
            writeFrame = getframe(f);
            writeVideo(myVid, writeFrame);
        end
        close(myVid)
    end

%% PCA

% Need rows = volumes and columns = pixels (linearized)

windTrials = myData.wholeSession(:,:,:,:,~logical(myData.stimSepTrials.windTrials));   % --> [x, y, plane, volume, trial]


%% PCA PLOTTING


% Pull out data for one plane
planeNum = 5;
planeData = squeeze(windTrials(:,:,planeNum,:,:)); % --> [x, y, volume, trial]
[w,x,y,z] = size(planeData);
planeDataRS = reshape(planeData, [w, x, y*z]); % --> [x, y, volume]


pcaData = mean(squeeze(myData.wholeSession(:,:,planeNum,:,:)),4); % --> [x, y, volume]
[n,m,d] = size(pcaData);
data2D = double(reshape(pcaData, [(n*m), d])); % --> [pixel, volume]

tData2D = data2D'; % --> [volume, pixel]

[coeff, score, latent, ~, explained] = pca(tData2D);

figure(1);clf;plot(explained(1:10))

coeffReshaped = reshape(coeff, [n, m, d-1]);

figure(2); clf;
subplot(2,2,1)
imshow(myData.refImg{planeNum},[0 myData.MAX_INTENSITY])
colormap(gca, 'gray')
% colormap('parula')
for iPlot = 2:4
    subplot(2, 2, iPlot); imagesc(coeffReshaped(:,:,iPlot-1));
    colormap(gca, 'bluewhitered')
end

figure(3); clf;
subplot(2,2,1)
for iPlot = 1:4
    subplot(2, 2, iPlot); imagesc(coeffReshaped(:,:,iPlot+3));
    colormap(gca, 'bluewhitered')
end

end%iFold

%%
%%% ICA
% 
% icaData = tData2D;
% nIC = 3;
% [icasig, A, W] = fastica(icaData, 'numOfIC', nIC);
% 
% output = W * icaData;
% % icasig = output;
% 
% icasigReshaped = reshape(icasig', [n, m, nIC]);
% 
% % imagesc(icasigReshaped(:,:,6))
% 
% figure(4); clf;
% subplot(2,2,1)
% imshow(mean(pcaData,3),[0 myData.MAX_INTENSITY])
% colormap(gca, 'gray')
% colormap('parula')
% for iPlot = 2:4
%     subplot(2, 2, iPlot); imagesc(icasigReshaped(:,:,iPlot-1));
% end
% 
% figure(5); clf;
% subplot(2,2,1)
% % colormap(gca, 'gray')
% for iPlot = 1:4
%    
%     subplot(2, 2, iPlot); imagesc(icasigReshaped(:,:,iPlot+3));
% %      colormap(gca, 'gray')
% end
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% CALCULATE MEAN dF/F FOR WIND STIM ONSET AND OFFSET
% stimSepTrials = []; wholeTrialAvg = []; baselineAvg = []; baselineF = []; dff = []; dffRaw = []; stimAvg = []; postStimAvg = [];
% for iStim = 1:length(stimTypes)
%     
%     % Separate trials by stimulus type
%     stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));
%     
%     %---------- Get trial averaged baseline and stimulus data ----------                                                                                              % myData.wholeSession = [x, y, plane, volume, trial]                                                                                                          
%     wholeTrialAvg(:,:,:,:,iStim) = mean(myData.wholeSession(:,:,:,:,stimSepTrials.(stimTypes{iStim})), 5);                                                            % --> [x, y, plane, volume, StimType]
%     
%     % Pre-stim
%     if floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate) > 0 
%         baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate):floor(stimStart*volumeRate),iStim), 4); % --> [x, y, plane, StimType]
%     else
%         % if stimDuration > preStimDuration, start baseline one second after beginning of trial
%         baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),iStim), 4);                                                 % --> [x, y, plane, StimType]
%     end
%     
%     % Stim period
%     stimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,ceil(stimStart*volumeRate):floor(stimEnd*volumeRate),iStim), 4);                                                  % --> [x, y, plane, StimType]
%     
%     % Post stim
%     if floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate) < size(wholeTrialAvg, 4)
%         postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate),iStim), 4);   % --> [x, y, plane, StimType]
%     else
%         % if stimEnd + stimLength > fullTrialDuration, end post-stim period at end of trial
%         postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):size(wholeTrialAvg, 4)-floor(volumeRate),iStim), 4);                            % --> [x, y, plane, StimType]
%     end
%     
% end%for
% 
% % Calculate dF/F values
% dffAvg = (stimAvg - baselineAvg) ./ baselineAvg; % --> [x, y, plane, StimType]
% dffAvgPost = (postStimAvg - stimAvg) ./ stimAvg; % --> [x, y, plane, StimType]
% 
% %% CALCULATE MEAN dF/F AROUND WIND RESPONSES
% stimSepTrials = [];
% 
% % Separate out  wind trials
% for iStim = 1:length(stimTypes)
%     stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
% end 
% windTrials = myData.wholeSession(:,:,:,:,logical(myData.stimSepTrials.windTrials));  % --> [x, y, plane, volume, trial]
% 
% % Calculate dF/F before and after wind onset using an equal period before onset as baseline
% stimLength = stimEnd - stimStart;
% stimLengthVols = floor(stimLength * volumeRate);
% 
% % Pre-stim
% if floor(stimStart*volumeRate) - stimLengthVols > 0
%     baselineVols = windTrials(:,:,:,floor(stimStart*volumeRate) - stimLengthVols:floor(stimStart*volumeRate),:);                                   % --> [x, y, plane, volume, trial]
% else
%     % If stimDuration > preStimDuration, start baseline one second after beginning of trial
%     baselineVols = windTrials(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),:);                                                              % --> [x, y, plane, volume, trial]
% end
% 
% % Stim period + 3 seconds after the stim offset
% stimVols = windTrials(:,:,:,ceil((stimStart*volumeRate)-size(baselineVols, 4)):floor((stimStart*volumeRate) + (3*volumeRate) + stimLengthVols),:); % --> [x, y, plane, volume, trial]  
% 
% stimVolsMean = mean(stimVols, 5);                                                                                                                  % --> [x, y, plane, volume]
% baselineMean = mean(mean(baselineVols, 5), 4);                                                                                                     % --> [x, y, plane]
% baselineMeanRep = repmat(baselineMean, 1, 1, 1, size(stimVols, 4));                                                                                % --> [x, y, plane, volume]
% windDffVols = (stimVolsMean - baselineMeanRep) ./ baselineMeanRep;                                                                                 % --> [x, y, plane, volume]
% 
% %% CREATE VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS
% 
% baselineLength = 12;
% 
% % Identify starting and ending indices of all movement bouts
% volActionCell = num2cell(volActions, 2);
% for iTrial = 1:myData.nTrials
%     volActionCell{iTrial} = num2str(volActionCell{iTrial}, '%u'); % Providing FormatSpec so it doesn't add spaces
% end
% baselineStr = num2str(zeros(1,baselineLength), '%u');
% [startInd, endInd] = regexp(volActionCell, [baselineStr, '2+']);
% 
% % Create array with information needed to identify bout data
% allBouts = [];
% for iTrial = 1:myData.nTrials
%     if ~isempty(startInd{iTrial})
%         for iBout = 1:length(startInd{iTrial})
%           allBouts(end+1,:) = [iTrial, startInd{iTrial}(iBout), endInd{iTrial}(iBout)]; % --> [trialNum, startInd, endInd]
%         end
%     end
% end
% 
% % Create save directory and open video writer
% savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
% fileName = 'Movement_Bouts_All_Trials';
% if ~isdir(savePath)
%     mkdir(savePath);
% end
% 
% % Warn user and offer to cancel save if this will overwrite an existing file
% overwrite = 1;
% if exist(fullfile(savePath, [fileName, '.avi']), 'file') ~= 0
%     dlgAns = questdlg('Creating this video will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
%     if strcmp(dlgAns, 'No')
%         overwrite = 0;
%         disp('Saving cancelled')
%     end
% end
% 
% % Create video writer
% myVid = VideoWriter(fullfile(savePath, fileName));
% myVid.FrameRate = 10;
% open(myVid)
% 
% for iBout = 1:size(allBouts,1)
%     
%     disp(['Plotting bout ', num2str(iBout), ' of ' num2str(size(allBouts, 1))]);
%     
%     % Calculate dF/F for the current bout ****Is this an acceptable way to calculate dF/F?*****
%     boutData = double(squeeze(myData.wholeSession(:,:,:, allBouts(iBout,2):allBouts(iBout,3),allBouts(iBout,1))));   % --> [x, y, plane, volume]
%     boutBaseline = mean(boutData(:,:,:,1:baselineLength),4);                                                         % --> [x, y, plane]
%     boutBaselineRep = repmat(boutBaseline,1,1,1,length(allBouts(iBout,2):allBouts(iBout,3)));                        % --> [x, y, plane, volume]
%     boutDff = (boutData - boutBaselineRep) ./ boutBaselineRep;                                                       % --> [x, y, plane, volume]
%     
%     boutDff(isinf(boutDff)) = 0; % To eliminate inf values from dividiing by zero above...baseline shouldn't be zero in valid data anyways
%     boutDff(isnan(boutDff)) = 0;
%     
%     % Avergage each volume with its neighbors to improve SNR
%     boutDff = movmean(boutDff, 3, 3);
%     
%     % Calculate dF/F value range
%     range = calc_range(boutDff,[]);
%     
%     for iVol = 1:size(boutDff, 4)
%         
%         % Create fig
%         f = figure(1); clf
%         f.Position = [50 45, 1620, 855];
%         
%         for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral
%             
%             % Plot dF/F for each plane
%             ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
%             imagesc(boutDff(:, :, iPlane, iVol))
%             caxis(range)
%             colormap(ax, bluewhitered) %bluewhitered % 'parula'
%             axis equal
%             axis off
%             
%             % Label positions
%             if iPlane == myData.nPlanes
%                 title('Ventral')
%             elseif iPlane == 1
%                 title('Dorsal')
%             end
%         end
%         
%         % Add title above all subplots
%         if iVol <= floor(baselineLength)
%             titleStr = 'Before movement onset';
%         else
%             titleStr = 'After movement onset';
%         end
%         
%         suptitle(['Trial #', num2str(allBouts(iBout,1)), ', Movement bout #', num2str(iBout), ' of ', ...
%             num2str(size(allBouts,1)), '  -  ', titleStr]);
%         
%         % Write frame to video
%         writeFrame = getframe(f);
%         writeVideo(myVid, writeFrame);
%     end
% end%if
% close(myVid)

