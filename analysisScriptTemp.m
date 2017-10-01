
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

% Create array of annotation data (row = trial, col = frame)
annotationArr = [];
annotTrials = 1:myData.nTrials;
for iTrial = annotTrials(annotData.goodTrials) % Skip any trials with dropped frames     ~logical(cellfun(@(x) strcmp(x, 'OdorNoWind'), myData.trialType)))
   annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.actionNums; 
end


% Average across trials
annotArrLog = annotationArr ~= 0;
annotArrMean = sum(annotArrLog, 1);

%----------Plot 2D data with seconds on X-axis----------

% Create figure and plot data
f = figure(1); clf; ax = axes();
f.Position = [200 100 1120 840]; 
imagesc(annotationArr);

% Format figure
f.Color = [1 1 1];
ax.FontSize = 12;
ax.XLabel.String = 'Time (sec)';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Trial number';
ax.YLabel.FontSize = 14;
ax.Title.String = [regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Behavior Summary']; % regex to add escape characters
ax.XTick = [0:(1/trialDuration):1]*nFrames;
ax.XTickLabel = [0:(1/trialDuration):1]*trialDuration;

%----------Plot total number of movement frames across trials----------

% Create figure and plot data
h = figure(2); clf; ax2 = axes();
h.Position = [100 200 1600 600];
plt = plot(smooth(annotArrMean,3));

% Add shading during stimulus presentation
yL = ylim();
rectPos = [stimStart*frameRate, yL(1), (stimEnd-stimStart)*frameRate, diff(yL)]; % [x y width height]
rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
ylim(yL);

% Format figure
h.Color = [1 1 1];
plt.LineWidth = 1;
ax2.FontSize = 12;
ax2.XLabel.String = 'Time (sec)';
ax2.XLabel.FontSize = 16;
ax2.YLabel.String = 'Total movement frames';
ax2.YLabel.FontSize = 14;
ax2.Title.String = [regexprep(expDate, '_(?<num>..)', '\\_$<num>'), '  Summed Movement Frames']; % regex to add escape characters
ax2.XTick = [0:(1/trialDuration):1]*nFrames;
ax2.XTickLabel = [0:(1/trialDuration):1]*trialDuration;
ax2.XLim = [0 nFrames];

% Add shading during stimulus presentation
yL = ylim();
rectPos = [stimStart*frameRate, yL(1), (stimEnd-stimStart)*frameRate, diff(yL)]; % [x y width height]
rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
ylim(yL);

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
% figure, clf; plot(onsetTimes, 1:96, '.');%, '.')
% set(gca, 'XTick', [0:(1/trialDuration):1]*nVols);
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
minMax = [];
minMax(1) = max(dffAvg(:));
minMax(2) = min(dffAvg(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

for iPlane = myData.nPlanes:-1:1 % Plot windows arranged from dorsal --> ventral
    
    % Create fig
    f = figure(iPlane); clf
    f.Position = [50 45, 1800, 950];
    
    % Plot reference image for the current plane
    ax1 = subplot(221);
    imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
       
    % Plot first stim type
    ax2 = subplot(222);
    imagesc(dffAvg(:,:,iPlane,1)) % [x, y, plane, StimType]
    caxis(range)
    colormap(ax2, bluewhitered) %bluewhitered % 'parula'
    axis equal 
    axis off
    title(stimTypes{1})
    colorbar
    
    % Second stim type
    ax3 = subplot(223);
    imagesc(dffAvg(:,:,iPlane,2))
    caxis(range)
    colormap(ax3, bluewhitered)
    axis equal 
    axis off
    title(stimTypes{2})
    colorbar
    
    % Third stim type
    ax4 = subplot(224);
    imagesc(dffAvg(:,:,iPlane,3))
    caxis(range)
    colormap(ax4, bluewhitered)
    axis equal 
    axis off
    title(stimTypes{3})
    colorbar
    
end

    %% PLOT WIND STIM OFFSET RESPONSE HEATMAPS FOR EACH PLANE
    
% Calculate absolute max dF/F value across all planes and stim types
minMax = [];
minMax(1) = max(dffAvgPost(:));
minMax(2) = min(dffAvgPost(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);
    
for iPlane = myData.nPlanes:-1:1 % Plot windows arranged from dorsal --> ventral
    
    % Create fig
    f = figure(iPlane); clf
    f.Position = [50 45, 1800, 950];
    
    % Plot reference image for the current plane
    ax1 = subplot(221);
    imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
       
    % Plot first stim type
    ax2 = subplot(222);
    imagesc(dffAvg(:,:,iPlane,1)) % [x, y, plane, StimType]
    caxis(range)
    colormap(ax2, bluewhitered) %bluewhitered % 'parula'
    axis equal 
    axis off
    title(stimTypes{1})
    colorbar
    
    % Second stim type
    ax3 = subplot(223);
    imagesc(dffAvg(:,:,iPlane,2))
    caxis(range)
    colormap(ax3, bluewhitered)
    axis equal 
    axis off
    title(stimTypes{2})
    colorbar
    
    % Third stim type
    ax4 = subplot(224);
    imagesc(dffAvg(:,:,iPlane,3))
    caxis(range)
    colormap(ax4, bluewhitered)
    axis equal 
    axis off
    title(stimTypes{3})
    colorbar
    
end


%% CALCULATE MEAN dF/F AROUND WIND RESPONSES

stimSepTrials = [];

% Separate out center wind trials
for iStim = 1:length(stimTypes)
    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
end 
centerWindTrials = myData.wholeSession(:,:,:,:,logical(stimSepTrials.CenterWind)); % --> [x, y, plane, volume, trial]

% Calculate dF/F before and after wind onset using an equal period before onset as baseline
stimLength = stimEnd - stimStart;
stimLengthVols = floor(stimLength * volumeRate);
if floor(stimStart*volumeRate) - stimLengthVols > 0
    baselineVols = centerWindTrials(:,:,:,floor(stimStart*volumeRate) - stimLengthVols:floor(stimStart*volumeRate),:); % --> [x, y, plane, volume, trial]
else
    % If stimDuration > preStimDuration, start baseline one second after beginning of trial
    baselineVols = centerWindTrials(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),:); % --> [x, y, plane, volume, trial]
end


stimVols = centerWindTrials(:,:,:,ceil((stimStart*volumeRate)-size(baselineVols, 4)):ceil((stimStart*volumeRate)+stimLengthVols),:); % --> [x, y, plane, volume, trial]  
stimVolsMean = mean(stimVols, 5);  % --> [x, y, plane, volume]
nVols = size(stimVols, 4);
baselineMean = mean(mean(baselineVols, 5), 4); % --> [x, y, plane]
baselineMeanRep = repmat(baselineMean, 1, 1, 1, nVols);     % --> [x, y, plane, volume]
windDffVols = (stimVolsMean - baselineMeanRep) ./ baselineMeanRep;  % --> [x, y, plane, volume]

% Calculate absolute max dF/F value across all planes and action states
minMax = [];
minMax(1) = max(windDffVols(:));
minMax(2) = min(windDffVols(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

    %% CREATE VIDEO OF MEAN dF/F THROUGHOUT WIND RESPONSES

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
fileName = '_Wind_Onset_Responses_Right+Center';
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
            suptitle('Before wind onset')
        else
            suptitle('After wind onset')
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

% Calculate absolute max dF/F value across all planes and action states
minMax = [];
minMax(1) = max(runDff(:));
minMax(2) = min(runDff(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

    %% PLOT BEHAVIORAL STATE HEATMAPS FOR EACH PLANE

for iPlane = myData.nPlanes:-1:1 % Plot windows arranged from dorsal --> ventral
     % Create fig
    f = figure(iPlane); clf
    f.Position = [50 45, 1800, 950];
    f.PaperOrientation = 'landscape';
    
    % Plot reference image for the current plane
    ax1 = subplot(121);
    imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
    
    % Plot dF/F for behavioral state
    ax2 = subplot(122);
    imagesc(runDff(:,:,iPlane))
    caxis(range)
    colormap(ax2, bluewhitered) %bluewhitered % 'parula'
    axis equal 
    axis off
    title('dF/F - Locomotion vs. Quiescent')
    colorbar

end


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
        runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel); % [trial, vol]
        stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel); % [trial, vol]
        legMoveVols(iTrial, :) = (volActions(iTrial,:) == legMoveLabel); % [trial, vol]
        
        % Find onsets of running bouts > 1 sec in duration and preceded by > 1 sec of quiescence
        baseLen = 6;
        respLen = 6;
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
minMax = [];
minMax(1) = max(onsetMeanDff(:));
minMax(2) = min(onsetMeanDff(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

for iPlane = myData.nPlanes:-1:1 % Plot windows arranged from dorsal --> ventral
    
    % Create fig
    f = figure(iPlane); clf
    f.Position = [50 45, 1800, 950];
    
    % Plot reference image for the current plane
    ax1 = subplot(121);
    imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
    
    % Plot dF/F for behavioral state
    ax2 = subplot(122);
    imagesc(onsetMeanDff(:,:,iPlane))
    caxis(range)
    colormap(ax2, bluewhitered) %bluewhitered % 'parula'
    axis equal
    axis off
    title('dF/F - Locomotion onset')
    colorbar
end

    %% CREATE VIDEO OF MEAN dF/F FOR EACH PLANE THROUGHOUT MOVEMENT ONSET

% Calculate absolute max dF/F value across all planes and action states
minMax = [];
minMax(1) = max(onsetDffVols(:));
minMax(2) = min(onsetDffVols(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

% Create save directory and open video writer
savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0\Analysis'];
fileName = 'Move_Onset_Responses_1_1';
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
            suptitle('Before movement onset')
        else
            suptitle('After movement onset')
        end
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end
    close(myVid)
end

    %% CREATE VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

baselineLength = 10;

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
myVid.FrameRate = 7;
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
    
    % Average every two volumes together to improve SNR
    n = 2;
    avgDff = [];
    count = 1;
    for iPlane = 1:n:size(boutDff,4)-n+1
        avgDff(:,:,:,count) = mean(boutDff(:,:,:,iPlane:iPlane+n-1),4);
        count = count + 1;
    end 
    boutDff = avgDff;
    
    % Calculate dF/F value range
    minMax = []; range = [];
    minMax(1) = max(boutDff(:));
    minMax(2) = min(boutDff(:));
    range(2) = max([abs(min(minMax)) max(minMax)]);
    range(1) = -range(2);
%     range = range * 0.5;
    
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
        if iVol <= floor(baselineLength/2)
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

%% COMBINED PLOTTING VIDEO WITH FLY BEHAVIOR


    %% CREATE VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

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
%     
    boutDff = movmean(boutDff, 3, 3);
%     % Average every two volumes together to improve SNR
%     n = 2;
%     avgDff = [];
%     count = 1;
%     for iPlane = 1:n:size(boutDff,4)-n+1
%         avgDff(:,:,:,count) = mean(boutDff(:,:,:,iPlane:iPlane+n-1),4);
%         count = count + 1;
%     end 
%     boutDff = avgDff;
    
    % Calculate dF/F value range
    minMax = []; range = [];
    minMax(1) = max(boutDff(:));
    minMax(2) = min(boutDff(:));
    range(2) = max([abs(min(minMax)) max(minMax)]);
    range(1) = -range(2);

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
            titleStr = 'Before movement onset';
        else
            titleStr = 'After movement onset';
        end
        suptitle(['Time = ', sprintf('%05.2f', frameTimes(boutFrames(iFrame))), ' sec  -  Trial #', ...
            num2str(allBouts(iBout,1)), ', Movement bout #', num2str(iBout), ' of ', ...
            num2str(size(allBouts,1)), '  -  ', titleStr]);
        
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
centerWindTrials = myData.wholeSession(:,:,:,:,stimSepTrials.CenterWind); % --> [x, y, plane, volume, trial]   
allWindTrials = myData.wholeSession(:,:,:,:,logical(stimSepTrials.CenterWind + stimSepTrials.RightWind + stimSepTrials.LeftWind)); % --> [x, y, plane, volume, trial]   
nVols = size(allWindTrials, 4);

% Calculate dF/F using whole trial average as baseline
baseline = mean(mean(allWindTrials, 5), 4);                 % --> [x, y, plane]
baselineMeanRep = repmat(baseline, 1, 1, 1, nVols);         % --> [x, y, plane, volume]
trialAvg = mean(allWindTrials, 5);                          % --> [x, y, plane, volume]
trialDff = (trialAvg - baselineMeanRep) ./ baselineMeanRep; % --> [x, y, plane, volume]

% Calculate absolute max dF/F value across all planes and action states
minMax = [];
minMax(1) = max(trialDff(:));
minMax(2) = min(trialDff(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);
range = range * 0.25;

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
minMax = [];
minMax(1) = max(trialAvg(:));
minMax(2) = min(trialAvg(:));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);
range = range * 0.25;

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
        stimTime = 5;
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

 %% PLOT WIND STIM ONSET RESPONSE HEATMAPS FOR 3 TRIAL CONDITIONS
    
% Calculate absolute max dF/F value across all planes and stim types
currConds = [1 2 5];
minMax = [];
minMax(1) = max(reshape(dffAvg(:,:,:,currConds), [1 numel(dffAvg(:,:,:,currConds))]));
minMax(2) = min(reshape(dffAvg(:,:,:,currConds), [1 numel(dffAvg(:,:,:,currConds))]));
range(2) = max([abs(min(minMax)) max(minMax)]);
range(1) = -range(2);

for iPlane = myData.nPlanes:-1:1 % Plot windows arranged from dorsal --> ventral
    
    % Create fig
    f = figure(iPlane); clf
    f.Position = [50 45, 1800, 950];
    
    % Plot reference image for the current plane
    ax1 = subplot(221);
    imshow(myData.refImg{iPlane}, [0 myData.maxIntensity])
       
    % Plot first stim type
    ax2 = subplot(222);
    imagesc(dffAvg(:,:,iPlane,currConds(1))) % [x, y, plane, StimType]
    caxis(range)
    colormap(ax2, bluewhitered) %bluewhitered % 'parula'
    axis equal 
    axis off
    title(trialCondNames{currConds(1)})
    colorbar
    
    % Second stim type
    ax3 = subplot(223);
    imagesc(dffAvg(:,:,iPlane,currConds(2)))
    caxis(range)
    colormap(ax3, bluewhitered)
    axis equal 
    axis off
    title(trialCondNames{currConds(2)})
    colorbar
    
    % Third stim type
    ax4 = subplot(224);
    imagesc(dffAvg(:,:,iPlane,currConds(3)))
    caxis(range)
    colormap(ax4, bluewhitered)
    axis equal 
    axis off
    title(trialCondNames{currConds(3)})
    colorbar
    
end







