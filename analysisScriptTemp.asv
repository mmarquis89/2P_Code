
%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA

% Prompt user for data file
[dataFile, pathName, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
if dataFile == 0
    disp('Initialization cancelled')
    imgData = []; % Skip loading if user clicked "Cancel"
else
    disp(['Loading ' dataFile, '...'])
    imgData = load([pathName, dataFile]);
    disp([dataFile, ' loaded'])
    
    % Prompt user for behavioral annotation data file
    [annotDataFile, pathName, ~] = uigetfile('*.mat', 'Select a behavioral annotation data file if desired', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
    if annotDataFile == 0
        disp('No behavioral annotation data selected')
        annotData = [];
    else
        disp(['Loading ' annotDataFile, '...'])
        annotData = load([pathName, annotDataFile]);
        disp([annotDataFile, ' loaded'])
    end
    myData = setstructfields(imgData, annotData); % Combine into one structure
    
    % Load .mat file containing trial metadata
    
    
    
    % Process raw data structure
    if ~isfield(myData, 'wholeSession') % myData.wholeSession = [x, y, plane, volume, trial]
        myData.wholeSession = myData.regProduct;
    end
    myData.nTrials = size(myData.wholeSession, 5);
    singleTrial = squeeze(myData.wholeSession(:,:,:,:,1)); % --> [x, y, plane, volume]
    myData.nPlanes = size(singleTrial, 3);
    myData.nVolumes = size(singleTrial, 4);
    
    % Create reference image for each plane
    myData.refImg = [];
    for iPlane = 1:myData.nPlanes
        myData.refImg{iPlane} = squeeze(mean(mean(myData.wholeSession(:,:,iPlane,:,:),4),5)); % --> [x, y]
    end
    
    % Extract session number
    origFileName = myData.origFileNames{1};
    sidLoc = strfind(origFileName, 'sid_');
    myData.sessionNum = origFileName(sidLoc+4);
    
    % For compatibility with early experiments
    if ~isfield(myData, 'expDate')
        cellDate = inputdlg('Please enter experiment date in YYYY_MM_DD format');
        myData.expDate = cellDate{1};
    end 

    % Add hardcoded parameters
    myData.volumeRate = 6.5; volumeRate = myData.volumeRate;
    myData.trialDuration = 15; trialDuration = myData.trialDuration;
    if ~isempty(annotData)
        myData.nFrames = max(cellfun(@height, myData.trialAnnotations)); nFrames = myData.nFrames;
    end
    myData.stimDuration = [4 trialDuration-7]; % [start time, length] in seconds
    myData.maxIntensity = 2000;
    stimStart = myData.stimDuration(1);
    stimEnd = sum(myData.stimDuration);
    stimTypes = sort(unique(myData.trialType));
    
end%if

%% PLOT AND SAVE VISUALIZATION OF BEHAVIOR DATA ANNOTATIONS

expDate = myData.expDate;
frameRate = 25;

% Separate out wind trials
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));
end
stimSepTrials.windTrials = logical(stimSepTrials.CenterWind + stimSepTrials.RightWind);

% Create array of annotation data (row = trial, col = frame)
annotationArr = [];
annotTrials = 1:myData.nTrials;
for iTrial = annotTrials(annotData.goodTrials) % Skip any trials with dropped frames     ~logical(cellfun(@(x) strcmp(x, 'OdorNoWind'), myData.trialType)))
   annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.actionNums; 
end

% Average across trials for wind and control trials
annotArrLog = annotationArr ~= 0;
annotArrMean_wind = sum(annotArrLog(stimSepTrials.windTrials,:), 1);
annotArrMean_noWind = sum(annotArrLog(~stimSepTrials.windTrials,:), 1);

%----------Plot 2D data with seconds on X-axis----------

% Create figure and plot data
f = figure(1); clf; ax = axes();
f.Position = [200 100 1120 840]; 
imagesc([annotationArr(stimSepTrials.windTrials,:) ; annotationArr(~stimSepTrials.windTrials,:)*1.25]);

% Format figure
f.Color = [1 1 1];
ax.FontSize = 12;
ax.XLabel.String = 'Time (sec)';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Trial number';
ax.YLabel.FontSize = 14;
ax.Title.String = [regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Behavior Summary  (top: wind trials,  bottom: control trials)']; % regex to add escape characters
ax.XTick = [0:(1/trialDuration):1]*nFrames;
ax.XTickLabel = [0:(1/trialDuration):1]*trialDuration;

%----------Plot total number of movement frames across trials for wind stim and control trials ----------

%----- Create figure and plot wind trial data -----
h = figure(2); clf;
h.Position = [100 50 1600 950];
subplot(2,1,1); ax2 = gca();
plt = plot(smooth(annotArrMean_wind,3));

% Format figure
h.Color = [1 1 1];
plt.LineWidth = 1;
ax2.FontSize = 12;
ax2.XLabel.String = 'Time (sec)';
ax2.XLabel.FontSize = 16;
ax2.YLabel.String = 'Total movement frames';
ax2.YLabel.FontSize = 14;
ax2.Title.String = 'Wind Trials';
ax2.XTick = [0:(1/trialDuration):1]*nFrames;
ax2.XTickLabel = [0:(1/trialDuration):1]*trialDuration;
ax2.XLim = [0 nFrames];

% Add shading during stimulus presentation
yL = ylim();
rectPos = [stimStart*frameRate, yL(1), (stimEnd-stimStart)*frameRate, diff(yL)]; % [x y width height]
rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
ylim(yL);

% ----- Plot control trial data-----
subplot(2,1,2); ax2 = gca();
plt = plot(smooth(annotArrMean_noWind,3));

% Format figure
h.Color = [1 1 1];
plt.LineWidth = 1;
ax2.FontSize = 12;
ax2.XLabel.String = 'Time (sec)';
ax2.XLabel.FontSize = 16;
ax2.YLabel.String = 'Total movement frames';
ax2.YLabel.FontSize = 14;
ax2.Title.String = 'Control Trials';
ax2.XTick = [0:(1/trialDuration):1]*nFrames;
ax2.XTickLabel = [0:(1/trialDuration):1]*trialDuration;
ax2.XLim = [0 nFrames];
suptitle([regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Summed Movement Frames']) % regex to add escape characters



% %----------Plot onset time of each movement bout----------
% 
% runVols = zeros(myData.nTrials, myData.nVolumes); 
% stoppedVols = runVols; legMoveVols = runVols; volActions = runVols;
% onsetVols = [];
% for iTrial = 1:myData.nTrials
%     
%     if myData.goodTrials(iTrial)
%         
%         % Match frame times to volumes
%         volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
%         frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
%         for iVol = 1:myData.nVolumes
%             [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
%             volFrames = volFrames';
%         end
%         
%         % Pull out action numbers for each volume
%         currActions = myData.trialAnnotations{iTrial}.actionNums;
%         volActions(iTrial,:) = currActions(volFrames);
%         
%         % Identify volume actions
%         locomotionLabel = 2; noActionLabel = 0; legMoveLabel = 4;
%         runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel); % [trial, vol]
%         stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel); % [trial, vol]
%         legMoveVols(iTrial, :) = (volActions(iTrial,:) == legMoveLabel); % [trial, vol]
%         
%         % Find all onsets of running bouts 
%         actionPattern = [0 2];
%         patternLen = length(actionPattern); % Make sure this is an even number
%         currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);
%         onsetVols{iTrial} = currTrialOnsets;
% 
%     else
%         % So data from invalid trials won't ever be matched to an action state
%         runVols(iTrial, :) = 0;
%         stoppedVols(iTrial, :) = 0;
%         legMoveVols(iTrial, :) = 0;
%         volActions(iTrial, :) = -1000;
%         onsetVols{iTrial} = [];
%     end
% end
% 
% onsetTimes = [onsetVols{:}];
% figure, clf; plot(onsetTimes, 1:numel(onsetTimes), '.');%, '.')
% set(gca, 'XTick', [0:(1/trialDuration):1]*myData.nVolumes);
% set(gca, 'XTickLabel', [0:(1/trialDuration):1]*trialDuration);



% ----------Save data as a .fig file and a .png file for each figure----------

% Create analysis directory if necessary
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sessionNum, '\Analysis'];
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

%% CALCULATE MEAN dF/F FOR WIND STIM ONSET AND OFFSET
stimSepTrials = []; wholeTrialAvg = []; baselineAvg = []; baselineF = []; dff = []; dffRaw = []; stimAvg = [];
for iStim = 1:length(stimTypes)
    
    % Separate trials by stimulus type
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));
    
    %---------- Get trial averaged baseline and stimulus data ----------                                                                                              % myData.wholeSession = [x, y, plane, volume, trial]                                                                                                          
    wholeTrialAvg(:,:,:,:,iStim) = mean(myData.wholeSession(:,:,:,:,stimSepTrials.(stimTypes{iStim})), 5);                                                            % --> [x, y, plane, volume, StimType]
    
    % Pre-stim
    if floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate) > 0 
        baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate):floor(stimStart*volumeRate),iStim), 4); % --> [x, y, plane, StimType]
    else
        % if stimDuration > preStimDuration, start baseline one second after beginning of trial
        baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),iStim), 4);                                                 % --> [x, y, plane, StimType]
    end
    
    % Stim period
    stimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,ceil(stimStart*volumeRate):floor(stimEnd*volumeRate),iStim), 4);                                                  % --> [x, y, plane, StimType]
    
    % Post stim
    if floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate) < size(wholeTrialAvg, 4)
        postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate),iStim), 4);   % --> [x, y, plane, StimType]
    else
        % if stimEnd + stimLength > fullTrialDuration, end post-stim period at end of trial
        postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):size(wholeTrialAvg, 4)-floor(volumeRate),iStim), 4);                            % --> [x, y, plane, StimType]
    end
    
end%for

% Calculate dF/F values
dffAvg = (stimAvg - baselineAvg) ./ baselineAvg; % --> [x, y, plane, StimType]
dffAvgPost = (postStimAvg - stimAvg) ./ stimAvg; % --> [x, y, plane, StimType]

    %% PLOT WIND STIM ONSET RESPONSE HEATMAPS FOR EACH PLANE
    
% Calculate absolute max dF/F value across all planes and stim types
range = calc_range(dffAvg, []);

% Plot figures
[f, plotAxes] = plot_heatmaps(dffAvg, myData, range, stimTypes, [], []);


    %% PLOT WIND STIM OFFSET RESPONSE HEATMAPS FOR EACH PLANE
    
% Calculate absolute max dF/F value across all planes and stim types
range = calc_range(dffAvgPost,[]);
    
% Plot figures
[f, plotAxes] = plot_heatmaps(dffAvgPost, myData, range, stimTypes, [], []);


%% CALCULATE MEAN dF/F AROUND WIND RESPONSES

stimSepTrials = [];

% Separate out center wind trials
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
end 
centerWindTrials = myData.wholeSession(:,:,:,:,logical(stimSepTrials.CenterWind));                                                   % --> [x, y, plane, volume, trial]

% Calculate dF/F before and after wind onset using an equal period before onset as baseline
stimLength = stimEnd - stimStart;
stimLengthVols = floor(stimLength * volumeRate);
if floor(stimStart*volumeRate) - stimLengthVols > 0
    baselineVols = centerWindTrials(:,:,:,floor(stimStart*volumeRate) - stimLengthVols:floor(stimStart*volumeRate),:);               % --> [x, y, plane, volume, trial]
else
    % If stimDuration > preStimDuration, start baseline one second after beginning of trial
    baselineVols = centerWindTrials(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),:);                                          % --> [x, y, plane, volume, trial]
end


stimVols = centerWindTrials(:,:,:,ceil((stimStart*volumeRate)-size(baselineVols, 4)):ceil((stimStart*volumeRate)+2*volumeRate +stimLengthVols),:); % --> [x, y, plane, volume, trial]  
stimVolsMean = mean(stimVols, 5);                                                                                                    % --> [x, y, plane, volume]
nVols = size(stimVols, 4);
baselineMean = mean(mean(baselineVols, 5), 4);                                                                                       % --> [x, y, plane]
baselineMeanRep = repmat(baselineMean, 1, 1, 1, nVols);                                                                              % --> [x, y, plane, volume]
windDffVols = (stimVolsMean - baselineMeanRep) ./ baselineMeanRep;                                                                   % --> [x, y, plane, volume]

% Calculate absolute max dF/F value across all planes and action states
range = calc_range(windDffVols,[]);

    %% CREATE VIDEO OF MEAN dF/F THROUGHOUT WIND RESPONSES

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
fileName = 'Wind_Onset_Responses2';
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
        
        % Create fig
        f = figure(iPlane); clf
        f.Position = [50 45, 1620, 855];
        
        %     % Plot reference image for the current plane
        %     ax1 = subplot(3,5,1);
        %     imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
        
        for iPlane = myData.nPlanes:-1:1 % Plot planes from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
            imagesc(windDffVols(:, :, iPlane, iVol))
            caxis(range)
            colormap(ax, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off
            
            % Label positions
            if iPlane == myData.nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
        end
        
        % Add title above all subplots
        if iVol <= size(baselineVols, 4)
            suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before wind onset'])
        else
            suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after wind onset'])
        end
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
    close(myVid)
end


%% CALCULATE MEAN dF/F ACROSS BEHAVIORAL STATES

% Identify behavioral state during each volume
runVols = zeros(myData.nTrials, myData.nVolumes); stoppedVols = runVols;
for iTrial = 1:myData.nTrials
    if myData.goodTrials(iTrial)
        
        % Match frame times to volumes
        volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
        frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
        for iVol = 1:myData.nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
            volFrames = volFrames';
        end
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions = currActions(volFrames);
        
        % Identify volume actions
        locomotionLabel = 2; noActionLabel = 0;
        runVols(iTrial, :) = (volActions == locomotionLabel); %--> [trial, vol]
        stoppedVols(iTrial, :) = (volActions == noActionLabel); %--> [trial, vol]
    else
        % So data from invalid trials won't ever be matched to an action state
        runVols(iTrial, :) = 0;
        stoppedVols(iTrial, :) = 0;
    end
end

% Calculate average values for each plane across behavioral states
meanRunVols = [];
meanStoppedVols = [];
imgData = myData.wholeSession;
for iTrial = 1:myData.nTrials
   currRunVols = logical(runVols(iTrial,:));
   currStoppedvols = logical(stoppedVols(iTrial,:));
   
   % Pull out running volumes, if any exist, from the current trial
   if sum(currRunVols) > 0
       meanRunVols(:,:,:,end+1) = mean(imgData(:,:,:,currRunVols,iTrial),4); %--> [x, y, plane, trial]
   end 
   
   % Pull out stopping volumes, if any exist, from the current trial
   if sum(currStoppedvols) > 0
       meanStoppedVols(:,:,:,end+1) = mean(imgData(:,:,:,currStoppedvols,iTrial),4); %--> [x, y, plane, trial]
   end  
end
runMean = mean(meanRunVols, 4); %--> [x, y, plane] 
stoppedMean = mean(meanStoppedVols, 4); %--> [x, y, plane] 

% Get dF/F values for running relative to quiescence
runDff = (runMean - stoppedMean) ./ stoppedMean; % --> [x, y, plane]

    %% PLOT BEHAVIORAL STATE HEATMAPS FOR EACH PLANE
    
% Calculate absolute max dF/F value across all planes and action states
range = calc_range(runDff,[]);

% Plot figures
[f, plotAxes] = plot_heatmaps(runDff, myData, range, {'dF/F - Locomotion vs. Quiescent'}, [], []);


%% CALCULATE MEAN dF/F AROUND LOCOMOTION ONSET

%---------- Identify behavioral state during each volume ----------
runVols = zeros(myData.nTrials, myData.nVolumes); 
stoppedVols = runVols; legMoveVols = runVols; volActions = runVols;
onsetVols = [];
for iTrial = 1:myData.nTrials
    
    if myData.goodTrials(iTrial)
        % Match frame times to volumes
        volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
        frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
        for iVol = 1:myData.nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
            volFrames = volFrames';
        end
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions(iTrial,:) = currActions(volFrames);
        
        % Identify volume actions
        locomotionLabel = 2; noActionLabel = 0; legMoveLabel = 4;
        runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel);   % [trial, vol]
        stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel); % [trial, vol]
        legMoveVols(iTrial, :) = (volActions(iTrial,:) == legMoveLabel);  % [trial, vol]
        
        % Find onsets of running bouts > 1 sec in duration and preceded by > 1 sec of quiescence
        baseLen = 12;
        respLen = 12;
        actionPattern = [zeros(1,baseLen), ones(1,respLen)*2];   % [0 0 0 0 0 0 0 2 2 2 2 2 2 2];
        patternLen = length(actionPattern);
        currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);
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
onsetBaselines = onsetData(:,:,:, 1:baseLen,:);                            % --> [x, y, plane, onsetVolume, onsetNum]
onsetBaselineMean = mean(mean(onsetBaselines, 5), 4);                      % --> [x, y, plane]
onsetBaselineMeanRep = repmat(onsetBaselineMean, 1, 1, 1, patternLen);     % --> [x, y, plane, onsetVolume]
onsetMean = mean(onsetData, 5);                                            % --> [x, y, plane, onsetVolume]
onsetDffVols = (onsetMean - onsetBaselineMeanRep) ./ onsetBaselineMeanRep; % --> [x, y, plane, onsetVolume]
onsetMeanDff = mean(onsetDffVols, 4);                                      % --> [x, y, plane]

    %% PLOT MOVEMENT ONSET HEATMAPS FOR EACH PLANE

% Calculate dF/F range
rangeScalar = 0.75;
range = calc_range(onsetMeanDff, rangeScalar);

% Plot figures
[f, plotAxes] = plot_heatmaps(onsetMeanDff, myData, range, {'dF/F - Locomotion onset'}, [], []);


    %% CREATE VIDEO OF MEAN dF/F FOR EACH PLANE THROUGHOUT MOVEMENT ONSET

% Calculate absolute max dF/F value across all planes and action states
rangeScalar = 0.4;
range = calc_range(onsetDffVols, rangeScalar);

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
fileName = 'Move_Onset_Responses_1_2';
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
    baselineVolTimes = -(1:size(onsetBaselines, 4))/myData.volumeRate;
    stimVolTimes = ((1:(size(onsetData, 4)-size(onsetBaselines,4)))/myData.volumeRate);
    relTimes = [baselineVolTimes(end:-1:1), stimVolTimes];
    
    % Create video writer
    myVid = VideoWriter(fullfile(savePath, fileName));
    myVid.FrameRate = 5;
    open(myVid)
    
    for iVol = 1:patternLen
        
        % Create fig
        f = figure(iPlane); clf
        f.Position = [50 45, 1620, 855];
        
        %     % Plot reference image for the current plane
        %     ax1 = subplot(3,5,1);
        %     imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
        
        for iPlane = myData.nPlanes:-1:1 % Plot planes from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
            imagesc(onsetDffVols(:, :, iPlane, iVol))
            caxis(range)
            colormap(ax, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off
            
            % Label positions
            if iPlane == myData.nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
        end
        
        % Add title above all subplots
        if iVol <= baseLen
            suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before movement onset'])
        else
            suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after movement onset'])
        end
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
    close(myVid)
end

    %% CREATE VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

baselineLength = 12;

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

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
fileName = 'Movement_Bouts_All_Trials';
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
myVid.FrameRate = 10;
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
    
    % Avergage each volume with its neighbors to improve SNR
    boutDff = movmean(boutDff, 3, 3);
    
    % Calculate dF/F value range
    range = calc_range(boutDff,[]);
    
    for iVol = 1:size(boutDff, 4)
        
        % Create fig
        f = figure(1); clf
        f.Position = [50 45, 1620, 855];
        
        for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
            imagesc(boutDff(:, :, iPlane, iVol))
            caxis(range)
            colormap(ax, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off
            
            % Label positions
            if iPlane == myData.nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
        end
        
        % Add title above all subplots
        if iVol <= floor(baselineLength)
            titleStr = 'Before movement onset';
        else
            titleStr = 'After movement onset';
        end
        
        suptitle(['Trial #', num2str(allBouts(iBout,1)), ', Movement bout #', num2str(iBout), ' of ', ...
            num2str(size(allBouts,1)), '  -  ', titleStr]);
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
end%if
close(myVid)


    %% CREATE COMBINED PLOTTING VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

baselineLength = 12; % Number of baseline volumes to analyze/plot for each movement bout

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
    range = calc_range(boutDff,[]);

    % Load behavior vid for the current trial
    vidDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', myData.expDate, '_Movies');
    trialStr = ['sid_', num2str(myData.sessionNum), '_tid_', sprintf('%03.f',(allBouts(iBout, 1)))];
    myMovie = [];
    flyVid = VideoReader(fullfile(vidDir, [trialStr, '.mp4']));
    while hasFrame(flyVid)
        currFrame = readFrame(flyVid);
        myMovie(:,:,end+1) = rgb2gray(currFrame);
    end
    myMovie = uint8(myMovie(:,:,2:end)); % Adds a black first frame for some reason, so drop that
    
    % Match frame times to volumes
    volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
    frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
    for iFrame = 1:size(myMovie, 3)
       [~, volFrames(iFrame)] = min(abs(volTimes - frameTimes(iFrame))); 
    end
    
    boutFrames = find(volFrames == allBouts(iBout,2),1):find(volFrames == allBouts(iBout,3),1);
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
        
        % Add title above all subplots
        if iFrame <= floor(find(volFrames == baselineLength, 1, 'last'))
            titleStr = 'before movement onset';
        else
            titleStr = 'after movement onset';
        end
        suptitle(['Trial #', num2str(allBouts(iBout,1)), ', Movement bout #', num2str(iBout), ' of ', ...
            num2str(size(allBouts,1)), '  -  Time = ', sprintf('%05.2f', frameTimes(boutFrames(iFrame))), ...
            ' sec ', titleStr]);
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
end%if
close(myVid)


%% CREATE VIDEO OF MEAN dF/F THROUGHOUT ENTIRE TRIAL

%----------Calculate whole-trial dF/F using overall mean as baseline----------

% Separate out wind trials
stimSepTrials = [];
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
end 
centerWindTrials = myData.wholeSession(:,:,:,:,stimSepTrials.CenterWind);                                   % --> [x, y, plane, volume, trial]   
allWindTrials = myData.wholeSession(:,:,:,:,logical(stimSepTrials.CenterWind + stimSepTrials.RightWind));   % --> [x, y, plane, volume, trial]   
nVols = size(allWindTrials, 4);

% Calculate dF/F using whole trial average as baseline
baseline = mean(mean(allWindTrials, 5), 4);                 % --> [x, y, plane]
baselineMeanRep = repmat(baseline, 1, 1, 1, nVols);         % --> [x, y, plane, volume]
trialAvg = mean(allWindTrials, 5);                          % --> [x, y, plane, volume]
trialDff = (trialAvg - baselineMeanRep) ./ baselineMeanRep; % --> [x, y, plane, volume]

% Calculate absolute max dF/F value across all planes and action states
rangeScalar = 0.25;
range = calc_range(trialDff, rangeScalar);

%----------Create video of dF/F for each plane throughout trial----------

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', myData.expDate, '\sid_0\Analysis'];
fileName = 'Full_Trial_dFF';
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
    % Create video writer
    myVid = VideoWriter(fullfile(savePath, fileName));
    myVid.FrameRate = 7;
    open(myVid)
    
    for iVol = 1:nVols
        
        % Create fig
        f = figure(iPlane); clf
        f.Position = [50 45, 1620, 855];
        
        for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
            imagesc(trialDff(:, :, iPlane, iVol))
            caxis(range)
            colormap(ax, bluewhitered) %bluewhitered % 'parula'
            axis equal
            axis off
            
            % Label positions
            if iPlane == myData.nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
        end
        
        % Add title above all subplots
        if iVol <= stimStart * volumeRate
            titleStr = 'Before wind onset';
        elseif iVol <= stimEnd * volumeRate
            titleStr = 'During wind stimulus';
        else
            titleStr = 'After wind stimulus';
        end
        suptitle([titleStr, '  -  time = ', num2str(round(iVol./volumeRate, 1))]);
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end%for
    close(myVid)
end%if
%----------Create video of mean raw fluorescence signal throughout trial----------

% Calculate range
rangeScalar = 0.25;
range = calc_range(trialAvg, rangeScalar);

% Update file name
fileName = 'Full_Trial_RawF';

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
    % Create video writer
    myVid = VideoWriter(fullfile(savePath, fileName));
    myVid.FrameRate = 7;
    open(myVid)
    
    for iVol = 1:nVols
        
        % Create fig
        f = figure(iPlane); clf
        f.Position = [50 45, 1620, 855];
        
        for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);
            imagesc(trialAvg(:, :, iPlane, iVol))
            caxis(range)
            colormap(ax, 'bluewhitered') %bluewhitered % 'parula'
            axis equal
            axis off
            
            % Label positions
            if iPlane == myData.nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
        end
        
        % Add title above all subplots
        if iVol <= stimStart * volumeRate
            titleStr = 'Before wind onset';
        elseif iVol <= stimEnd * volumeRate
            titleStr = 'During wind stimulus';
        else
            titleStr = 'After wind stimulus';
        end
        suptitle([titleStr, '  -  time = ', num2str(round(iVol./volumeRate, 1))]);
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end%for
    close(myVid)
end%if


%% SEPARATE TRIALS BASED ON WHETHER FLY WAS MOVING DURING WIND STIM

% Separate out wind trials
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));
end
stimSepTrials.windTrials = stimSepTrials.CenterWind + stimSepTrials.RightWind;

%---------- Identify behavioral state during each volume ----------
volActions = zeros(myData.nTrials, myData.nVolumes);
preStimMove = zeros(myData.nTrials, 1); stimMove = preStimMove;
volFrames = zeros(1,myData.nVolumes);
for iTrial = 1:myData.nTrials
    
    if myData.goodTrials(iTrial)
        
        % Match frame times to volumes
        volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
        frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
        for iVol = 1:myData.nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
        end
        volFrames = volFrames';
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions(iTrial,:) = currActions(volFrames);
        
        preStimTime = 2;
        stimTime = 3;
        preStimVols = floor((stimStart-preStimTime)*volumeRate):floor(stimStart*volumeRate);
        stimVols = floor(stimStart*volumeRate):floor((stimStart+stimTime)*volumeRate);
        
        % Determine if fly was active just before the wind onset
        if sum(volActions(iTrial, preStimVols)) > 0
            preStimMove(iTrial) = 1;
        end
        
        % Determine if fly was moving in the first few seconds after wind onset
        if sum(volActions(iTrial, stimVols)) > 0
            stimMove(iTrial) = 1;
        end
        
    else
        % So data from invalid trials won't ever be matched to an action state
        volActions(iTrial, :) = -1000;
    end        
end

% Separate trials
trialConditions = [preStimMove, stimMove.*3, stimSepTrials.windTrials'.*5];
windNoMove = (sum(trialConditions,2) == 5);
windStartMove = (sum(trialConditions,2) == 8);
windContMove = (sum(trialConditions,2) == 9);
noWindNoMove = (sum(trialConditions,2) == 0);
noWindStartMove = (sum(trialConditions,2) == 3);
noWindContMove = (sum(trialConditions,2) == 4);
trialConds = [windNoMove,windStartMove,windContMove,noWindNoMove,noWindStartMove,noWindContMove];
trialCondNames = {'windNoMove','windStartMove','windContMove','noWindNoMove','noWindStartMove','noWindContMove'};

%% CALCULATE MEAN dF/F FOR EACH TRIAL CONDITION

wholeTrialAvg = []; baselineAvg = []; dffAvg = []; stimAvg = [];
for iCond = 1:length(trialCondNames)
    
    %---------- Get trial averaged baseline and stimulus data ----------                                                                      % myData.wholeSession = [x, y, plane, volume, trial]                                                                                                          
    wholeTrialAvg(:,:,:,:,iCond) = mean(myData.wholeSession(:,:,:,:,trialConds(:,iCond)), 5);                                                 % --> [x, y, plane, volume, trialCondition]
    
    % Pre-stim
    if floor((stimStart-preStimTime)*volumeRate) > 0 
        baselineAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,floor((stimStart-preStimTime)*volumeRate):floor(stimStart*volumeRate),iCond), 4); % --> [x, y, plane, StimType]
    else
        % if stimDuration > preStimDuration, start baseline one second after beginning of trial
        baselineAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),iCond), 4);                         % --> [x, y, plane, trialCondition]
    end
    
    % Stim period
    stimAvg(:,:,:,iCond) = mean(wholeTrialAvg(:,:,:,ceil(stimStart*volumeRate):floor((stimStart+stimTime)*volumeRate),iCond), 4);             % --> [x, y, plane, trialCondition]
end%for

% Calculate dF/F values
dffAvg = (stimAvg - baselineAvg) ./ baselineAvg; % --> [x, y, plane, StimType]

% Eliminate inf values from dividiing by zero above...baseline shouldn't be zero in valid data anyways
dffAvg(isinf(dffAvg)) = 0;

 %% PLOT WIND STIM ONSET RESPONSE HEATMAPS FOR SOME TRIAL CONDITIONS
    
% Calculate absolute max dF/F value across all planes and stim types
currConds = [2 5 6];
rangeScalar = 0.25;
dffConds = reshape(dffAvg(:,:,:,currConds), [1 numel(dffAvg(:,:,:,currConds))]);
range = calc_range(dffConds, rangeScalar);

% Plot figures
dffCurrConds = dffAvg(:,:,:,currConds);
[f, plotAxes] = plot_heatmaps(dffCurrConds, myData, range, trialCondNames, [], []);



%% MAKE VIDEOS OF WIND AND MOVEMENT RESPONSES FOR JUST A SINGLE PLANE


baselineDur = 2;
responseDur = 4;
planeNum = 4;

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
        % Match frame times to volumes
        volTimes = (1:myData.nVolumes)' ./ myData.volumeRate;
        frameTimes = myData.trialAnnotations{find(myData.goodTrials, 1)}.frameTime;
        for iVol = 1:myData.nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
            volFrames = volFrames';
        end
        
        % Pull out action numbers for each volume
        currActions = myData.trialAnnotations{iTrial}.actionNums;
        volActions(iTrial,:) = currActions(volFrames);
        
        % Identify volume actions
        locomotionLabel = 2; noActionLabel = 0; legMoveLabel = 4;
        runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel);   % [trial, vol]
        stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel); % [trial, vol]
        legMoveVols(iTrial, :) = (volActions(iTrial,:) == legMoveLabel);  % [trial, vol]
        
        % Find onsets of running bouts
        baseLen = size(baselineVols, 4);
        respLen = size(stimVols,4) - size(baselineVols, 4);
        actionPattern = [zeros(1,baseLen), ones(1,respLen)*2];   % [0 0 0 0 0 0 0 2 2 2 2 2 2 2];
        patternLen = length(actionPattern);
        currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);
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
onsetBaselines = onsetData(:,:,:, 1:baseLen,:);                            % --> [x, y, plane, onsetVolume, onsetNum]
onsetBaselineMean = mean(mean(onsetBaselines, 5), 4);                      % --> [x, y, plane]
onsetBaselineMeanRep = repmat(onsetBaselineMean, 1, 1, 1, patternLen);     % --> [x, y, plane, onsetVolume]
onsetMean = mean(onsetData, 5);                                            % --> [x, y, plane, onsetVolume]
onsetDffVols = (onsetMean - onsetBaselineMeanRep) ./ onsetBaselineMeanRep; % --> [x, y, plane, onsetVolume]
onsetMeanDff = mean(onsetDffVols, 4);                                      % --> [x, y, plane]



%% MAKE COMBINED VIDEO
% Calculate dF/F ranges
windRange = calc_range(windDffVols, []);
moveRange = calc_range(onsetDffVols,[]);

% Create save directory and open video writer
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
        imshow(myData.refImg{planeNum}, [0 myData.maxIntensity])
        
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
        
%         % Add title above all subplots
%         if iVol <= size(baselineVols, 4)
%             suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol))])
%         else
%             suptitle(['Time = ', sprintf('%05.2f', relTimes(iVol))])
%         end
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
    close(myVid)
end



