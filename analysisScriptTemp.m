%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA

  
% Load .mat file containing trial data
myData = load_imaging_data();

badTrials = [];

for iFold = 1
    
    % Throw specifc trials out of the dataset if necessary
    if ~isempty(badTrials)
        myData.goodTrials(badTrials) = 0;
%         rmData = myData;
%         rmData.goodTrials(badTrials) = [];
%         rmData.nTrials = myData.nTrials - numel(badTrials); nTrials = rmData.nTrials;
%         rmData.origFileNames(badTrials) = [];
%         fNames = fieldnames(myData.stimSepTrials);
%         for iField = 1:numel(fNames)
%             rmData.stimSepTrials.(fNames{iField})(badTrials) = [];
%         end
%         rmData.trialAnnotations(badTrials) = [];
%         rmData.trialType(badTrials) = [];
%         rmData.wholeSession(:,:,:,:,badTrials) = [];
%         rmData.annotArr(badTrials,:,:) = [];
%         myData = rmData;
%         clear('rmData');
    end
    
    % Copy variables for convenience
    wholeSession = myData.wholeSession;
    sessionSize = size(wholeSession);
    expDate = myData.expDate;
    sid = myData.sid;
    nPlanes = myData.nPlanes;
    nVolumes = myData.nVolumes;
    refImg = myData.refImg;
    if ~isempty(myData.nFrames)
        nFrames = myData.nFrames;
    else
        nFrames = nVolumes; 
    end
    nTrials = myData.nTrials;
    nGoodTrials = sum(myData.goodTrials);
    stimTypes = myData.stimTypes;
    trialDuration = myData.trialDuration;
    volumeRate = myData.volumeRate;
    volFrames = myData.volFrames;
    goodTrials = myData.goodTrials;
    
    preStimDur = myData.preStimDur;
    postStimDur = myData.postStimDur;
    stimPeriodDur = myData.stimPeriodDur;
    
    ballStop = myData.ballStop;
    nStops = myData.nStops;
    stopDur = myData.stopDur;
    interStopInterval = myData.interStopInterval;
    
    odorStim = myData.odorStim;
    nOdorStims = myData.nOdorStims;
    odorStimDur = myData.odorStimDur;
    interOdorInterval = myData.interOdorInterval;
    odorStartTimes = myData.odorStartTimes;
    odorEndTimes = myData.odorEndTimes;
        
    % Create hardcoded parameters
    myData.ROIdata = [];
    myData.MAX_INTENSITY = 1000; MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots
    myData.FRAME_RATE = 25; FRAME_RATE = 25; % This is the frame rate of the behavior video, not the GCaMP imaging
    if isempty(nFrames)
        nFrames = sum(trialDuration) * FRAME_RATE;
    end
    volTimes = (1:nVolumes)' ./ volumeRate;
    frameTimes = (1:nFrames)' ./ FRAME_RATE;
    
end%iFold

%% VIEW RAW DATA FOR A SINGLE TRIAL AND PLANE

planeNum = 6;
trialNum = 8; % Does not account for any skipped trials

preview_trial_movie(myData.wholeSession, planeNum, trialNum, [], [], []);

%% INITIAL DATA PROCESSING STEPS

skipTrials = []; 
myData.skipTrials = skipTrials;
nSkippedTrials = length(skipTrials); myData.nSkippedTrials = nSkippedTrials;

%============== Create array of annotation/event data ==============================================

%   annotArr: (row = trial, col = frame, Z-dim = event type)
%   Event types are: [odor stim, behavior, ball stopping]

annotationTypes = [];

% ----------------------------------------------------------------------------------------------
% Odor stim
% ----------------------------------------------------------------------------------------------

% Add odor stim frames
odorAnnotArr = zeros(nTrials, nFrames);
if odorStim
    for iStim = 1:nOdorStims
        odorFrames = floor((odorStartTimes(iStim) * FRAME_RATE)):floor((odorEndTimes(iStim) * FRAME_RATE));
        goodOdorTrials = logical(myData.stimSepTrials.odorTrials .* goodTrials);
        odorAnnotArr(goodOdorTrials, odorFrames) = 4; %--> [trial, frame]
    end
end

% All odor events
odorAnnotations = annotationType(myData, odorAnnotArr, skipTrials, 'odor');
odorAnnotations = get_event_vols(odorAnnotations, '04', '40');
annotationTypes{end + 1} = odorAnnotations;

% Odor A events
annotArr_OdorA = odorAnnotArr;
annotArr_OdorA(~myData.stimSepTrials.OdorA, :) = 0;
odorAnnotations_A = annotationType(myData, annotArr_OdorA, skipTrials, 'odor_A');
odorAnnotations_A = get_event_vols(odorAnnotations_A, '04', '40');
annotationTypes{end + 1} = odorAnnotations_A;

% Odor B events
annotArr_OdorB = odorAnnotArr;
annotArr_OdorB(~myData.stimSepTrials.OdorB, :) = 0;
odorAnnotations_B = annotationType(myData, annotArr_OdorB, skipTrials, 'odor_B');
odorAnnotations_B = get_event_vols(odorAnnotations_B, '04', '40');
annotationTypes{end + 1} = odorAnnotations_B;

% ----------------------------------------------------------------------------------------------
% Behavior
% ----------------------------------------------------------------------------------------------

% Add behavior annotations
behaviorAnnotArr = zeros(nTrials, nFrames);
if ~isempty(myData.trialAnnotations)
    annotTrials = 1:nTrials;
    for iTrial = annotTrials(goodTrials) % Skip any trials with dropped frames
        behaviorAnnotArr(iTrial, :) = myData.trialAnnotations{iTrial}.actionNums;           %--> [trial, frame]
    end
end

% Locomotion events
locomotionAnnotations = annotationType(myData, behaviorAnnotArr, skipTrials, 'locomotion');
locomotionAnnotations = get_event_vols(locomotionAnnotations, '[034]2', '2[034]'); 
annotationTypes{end + 1} = locomotionAnnotations;

% Isolated movement events
isoMoveAnnotations = annotationType(myData, behaviorAnnotArr, skipTrials, 'isoMove');
isoMoveAnnotations = get_event_vols(isoMoveAnnotations, '[023]4', '4[023]'); 
annotationTypes{end + 1} = isoMoveAnnotations;

% Grooming events
groomingAnnotations = annotationType(myData, behaviorAnnotArr, skipTrials, 'groom');
groomingAnnotations = get_event_vols(groomingAnnotations, '[024]3', '3[024]');
annotationTypes{end + 1} = groomingAnnotations;

% All behavior events
onsetRegExpStr = '[034]2|[023]4|[024]3';
offsetRegExpStr = '2[034]|4[023]|3[024]';
behaviorAnnotations = annotationType(myData, behaviorAnnotArr, skipTrials, 'move');
behaviorAnnotations = get_event_vols(behaviorAnnotations, onsetRegExpStr, offsetRegExpStr);
annotationTypes{end + 1} = behaviorAnnotations;

% ----------------------------------------------------------------------------------------------
% Ball stopping
% ----------------------------------------------------------------------------------------------

% Add ball stop annotations
bStopAnnotArr = zeros(nTrials, nFrames);
if ~isempty(myData.trialAnnotations)
    if ballStop
        annotTrials = 1:nTrials;
        for iTrial = annotTrials(goodTrials) % Skip any trials with dropped frames
            bStopAnnotArr(iTrial, :) = myData.trialAnnotations{iTrial}.ballStopNums * 2; %--> [trial, frame]
        end
    end
end

% Ball stopping events
bStopAnnotations = annotationType(myData, bStopAnnotArr, skipTrials, 'bStop');
bStopAnnotations = get_event_vols(bStopAnnotations, '04', '40');
annotationTypes{end + 1} = bStopAnnotations;

% ==================================================================================================

% Make list of all annotation type names
for iType = 1:numel(annotationTypes)
    annotationTypeNames{iType} = annotationTypes{iType}.name;
end
annotationTypeSummary = table([1:numel(annotationTypeNames)]', annotationTypeNames', 'VariableNames', {'Index', 'AnnotationType'});


%% =================================================================================================
%           ROI-BASED ANALYSES                                   
%%%=================================================================================================
for iFold = 0
    
%% LOAD ROI DATA

parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(myData.sid), '\Analysis\ROIs'];
fileName = 'ROI_Data.mat';

load(fullfile(parentDir, fileName));
myData.ROIdata = ROIdata;
disp('ROIs loaded')

%% EXTRACT SESSION DATA WITHIN ROIs


nROIs = numel(myData.ROIdata.mask);
ROIDataAvg = [];
for iROI = 1:nROIs  
    disp(['Extracting data for ROI #', num2str(iROI), ' of ', num2str(nROIs), '...'])
    currMask = myData.ROIdata.mask{iROI};
    startInd = strfind(myData.ROIdata.plane{iROI}, '#') + 1;
    currPlane = str2double(myData.ROIdata.plane{iROI}(startInd:end));
    currPlaneData = squeeze(wholeSession(:,:,currPlane,:,:));                                   % --> [y, x, volume, trial]
    currPlaneData(~currMask(:,:,ones(1, nVolumes), ones(1, nTrials))) = nan;                    % --> [ y, x, volume, trial]
    
    currDataLin = reshape(currPlaneData, size(currPlaneData, 1)*size(currPlaneData, 2), ...
                                nVolumes, nTrials);                                             % --> [pixel, volume, trial, ROI]
    ROIDataAvg(:,:,iROI) = squeeze(mean(currDataLin, 1, 'omitnan'));                            % --> [volume, trial, ROI] 
end
disp('ROI extraction complete')
myData.ROIDataAvg = ROIDataAvg;


%% CALCULATE MEAN dF/F WITHIN ROIs THROUGHOUT ENTIRE EXPERIMENT

    
    % Using bottom 5% of entire ROI's mean value throughout each trial as baseline
    ROIDataAvgSorted = sort(ROIDataAvg, 1);                                     % --> [volume, trial, ROI] 
    baselineMean = mean(ROIDataAvgSorted(1:round(nVolumes * 0.05), :, :), 1);   % --> [volume, trial, ROI] 
    baselineMeanRep = baselineMean(ones(1, nVolumes), :, :);                    % --> [volume, trial, ROI] 
    ROIDffAvg = (ROIDataAvg - baselineMeanRep) ./ baselineMeanRep;              % --> [volume, trial, ROI] 
    myData.ROIDffAvg = ROIDffAvg;
    
    
    %% PLOT MEAN dF/F WITHIN ROIs THROUGHOUT EXPERIMENT
    
    % Concatate the trial-by-trial dF/F values into one long vector
    nROIs = size(ROIDffAvg, 3);
    concatROIDff = reshape(ROIDffAvg, nVolumes * nTrials, nROIs);   % --> [volume, ROI]
    
    % Get linear arrays of event annotations
    annotVols = annotationTypes{7}.eventVolsLin;
    locomotionVols = annotationTypes{4}.eventVolsLin;
    isoMoveVols = annotationTypes{5}.eventVolsLin;
    odorVols = annotationTypes{1}.eventVolsLin;
    for iROI = 1:nROIs
        
        f = figure(iROI); clf; hold on
        f.Position = [-1850 200 1800 600];
        title(num2str(iROI))
        
        % Plot dF/F and event annotations
        plot(smooth(concatROIDff(:,iROI), 3), 'b');
%         plot(isoMoveVols, 'c');
%         plot(locomotionVols, 'r');
        plot(odorVols, 'r');
        
        % Add trial delineators and numbers
        allVols = 1:(nTrials * nVolumes);
        yL = ylim();
        for iVol = allVols(~logical(mod(allVols, nVolumes)))
            plot([iVol, iVol], yL,  'Color', 'k')
            text(iVol + 50, 0.8 * yL(2), num2str(round(iVol / nVolumes)+1));
        end
        
        xlim([0, 3 * nVolumes]);
    end
    
 

%%

plotTitlePrefixes = { ...
    ''
};
% saveDir = 0;
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');


filterVecs = [myData.stimSepTrials.OdorA; myData.stimSepTrials.OdorB];

ROIlist = 1:size(ROIDffAvg, 3);
ROIlist = [1 2 3 4 5 6 7 8 9 10];

for iROI = ROIlist
    
    % Separate dF/F data for the current stim type
    currDffAvg = ROIDffAvg(:,logical(myData.stimSepTrials.OdorA), iROI); % --> [volume, trial]
    
    % Create figure
    f = figure(iROI); clf
    f.Position = [-1850 200 1800 600];
    f.Color = [1 1 1];
    
    % Plot reference image and ROI outline
    startInd = strfind(myData.ROIdata.plane{iROI}, '#') + 1;
    currPlane = str2double(myData.ROIdata.plane{iROI}(startInd:end));
    xData = myData.ROIdata.xi{iROI};
    yData = myData.ROIdata.yi{iROI};
    subaxis(3, 6, [1 2 7 8], 'ML', 0.005, 'MR', 0.01, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
    hold on
    imshow(myData.refImg{currPlane}, [0 MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
    title(myData.ROIdata.plane{iROI})
    
    % Odor A plot
    subaxis(3, 6, [3:6])
    ax1 = gca;
    plot_ROI_data(ax1, currDffAvg, myData, annotationTypes{1});
    yL1 = ylim(ax1);
    title('ACV (top), Empty vial (middle), 50 msec empty vial pulses (bottom)');
    
    % Odor B plot
    currDffAvg = ROIDffAvg(:,logical(myData.stimSepTrials.OdorB), iROI); % --> [volume, trial]
    subaxis(3, 6, [9:12])
    ax2 = gca;
    plot_ROI_data(ax2, currDffAvg, myData, annotationTypes{1});
    yL2 = ylim(ax2); 
    
    % Odor B Pulse plot
    currDffAvg = ROIDffAvg(:,logical(myData.stimSepTrials.OdorBPulse), iROI); % --> [volume, trial]
    subaxis(3, 6, [15:18])
    ax3 = gca;
    plot_ROI_data(ax3, currDffAvg, myData, annotationTypes{1});
    yL3 = ylim(ax3); 
    
    % Scale y-axes to match
    yMin = min([yL1, yL2, yL3]);
    yMax = max([yL1, yL2, yL3]);
    ylim(ax1, [yMin yMax]);
    ylim(ax2, [yMin yMax]);
    
     
%     
%     % Odor A plot ------------------------------------------------------------------------
%    
%     hold on
%     ylabel('dF/F');
%     xlabel('Time (s)')
%     xlim([0, max(volTimes)]);
%     
%     currDffAvg = ROIDffAvg(:,logical(myData.stimSepTrials.OdorA), iROI); % --> [volume, trial]
%     odorADff = mean(currDffAvg, 2);
%     stdDev = std(currDffAvg, 0, 2);
%     
%     % Discard any trials that are >5 SDs from mean
%     outliers = zeros(1, size(currDffAvg, 2));
%     for iTrial = 1:size(currDffAvg, 2)
%         if sum(abs(currDffAvg(:, iTrial) - odorADff) > (5 * stdDev))
%             outliers(iTrial) = 1;
%         end
%     end
%     currDffAvg(:, logical(outliers)) = [];
%     odorADff = mean(currDffAvg, 2);
%     stdDev = std(currDffAvg, 0, 2);
%     
%     % Plot individual trials in background
%     for iTrial = 1:size(currDffAvg, 2)
%         currData = currDffAvg(:, iTrial);
%                 plot(volTimes, smooth(currData, 3), 'color', [0 0 1 0.25], 'LineWidth', 0.1)
%     end
%        
%     % Shade one SD above and below mean
%     upper = smooth(odorADff, 3) + stdDev;
%     lower = smooth(odorADff, 3) - stdDev;
%     jbfill(volTimes, upper', lower', 'b', 'b', 1, 0.2);
%     
%     % Plot mean response line
%     plot(volTimes, smooth(odorADff, 3), 'LineWidth', 2, 'Color', 'b');
%     
%     % Plot odor stim timing
%     yL = ylim();
%     for iOnset = 1:numel(odorOnsetTimes)
%         currTime = odorOnsetTimes(iOnset);
%         fill([currTime, currTime, currTime + 1, currTime + 1], [yL(1), yL(2), yL(2), yL(1)], 'r', 'facealpha', 0.20, 'edgealpha', 0);
%     end
%     title('ACV (top), Empty vial (bottom)')
%     
%     % Odor B plot --------------------------------------------------------------------------
%     subaxis(2, 6, [9:12])
%     hold on
%     ylabel('dF/F');
%     xlim([0, max(volTimes)]);
%     
%     currDffAvg = ROIDffAvg(:,logical(myData.stimSepTrials.OdorB), iROI); % --> [volume, trial]
%     stdDev = std(currDffAvg, 0, 2);
%     odorBDff = mean(currDffAvg, 2);
%     
%     % Discard any trials that are >5 SDs from mean
%     outliers = zeros(1, size(currDffAvg, 2));
%     for iTrial = 1:size(currDffAvg, 2)
%         if sum(abs(currDffAvg(:, iTrial)-odorBDff) > (5 * stdDev))
%             outliers(iTrial) = 1;
%         end
%     end
%     currDffAvg(:, logical(outliers)) = [];
%     odorBDff = mean(currDffAvg, 2);
%     stdDev = std(currDffAvg, 0, 2);
%     
%     % Plot individual trials in background
%     for iTrial = 1:size(currDffAvg, 2)
%         currData = currDffAvg(:, iTrial);
%                 plot(volTimes, smooth(currData, 3), 'color', [0 0 1 0.25], 'LineWidth', 0.1)
%     end
%        
%     % Shade one SD above and below mean
%     upper = smooth(odorBDff, 3) + stdDev;
%     lower = smooth(odorBDff, 3) - stdDev;
%     jbfill(volTimes, upper', lower', 'b', 'b', 1, 0.2);
%     
%     % Plot mean response line
%     odorVols = annotationTypes{1}.eventVolsLin;
%     plot(volTimes, smooth(odorBDff, 3), 'LineWidth', 2, 'Color', 'b');
%     plot(volTimes, odorVols, 'r', 'LineWidth', 1);
%     
    
    % Save figure -------------------------------------------------------------------------------
    if saveDir
        fileName = ['ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_OdorResponses'];
        savefig(f, fullfile(saveDir, fileName));
        export_fig(fullfile(saveDir, fileName), '-png', f);
    end
end

end%iFold 


%% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS

saveFig = 0;
plotTypes = [2 1]; % 1 = odor stims, 2 = behavior, 3 = ball stopping
s = myData.stimSepTrials;
trialGroups =  [[s.OdorA + 2 * s.OdorB + 3 * s.OdorBPulse] .* goodTrials];%[]%            %[s.OdorA + 2 * s.OdorB + 3 * s.NoOdor] .* goodTrials; %[s.odorTrials + 2 * (~s.odorTrials)]; % 
plotTitleSuffix = '(ACV vs Empty Vial vs 50 msec air pulses)';%''%
fileNameSuffix = '_OdorAvsOdorBvsAirPulse'; %'_AllTrials';%
for iFold = 1

% Create plot titles
nPlots = length(plotTypes);
titleStrings = [];
plotNames = [];
annotArr = [];
for iPlot = 1:nPlots
    if plotTypes(iPlot) == 1
        % Odor stim
        plotNames{iPlot} = 'Odor Delivery';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' Summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = odorAnnotArr;
    elseif plotTypes(iPlot) == 2
         % Behavior
        plotNames{iPlot} = 'Behavior Annotation';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' Summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = behaviorAnnotArr;
    elseif plotTypes(iPlot) == 3
        % Ball stopping
        plotNames{iPlot} = 'Ball Stopping';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' Summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = bStopAnnotArr;
    end%if
end%for

% Create figure
f = figure(7);clf
if nPlots == 1
    f.Position = [200 45 1020 950];
elseif nPlots == 2
    f.Position = [200 45 900 950];
elseif nPlots == 3
    f.Position = [200 45 700 950];
end
f.Color = [1 1 1];

% Plot and format each figure
for iPlot = 1:nPlots
    
    subaxis(nPlots, 1, iPlot, ...
            'MarginTop', 0, ...
            'MarginBottom', 0.055, ...
            'MarginRight', 0.015, ...
            'MarginLeft', 0.065, ...
            'Spacing', 0, ...
            'PaddingTop', 0.03 ...
            );
    ax = gca;
    [~, ax, ~] = plot_behavior_summary_2D(myData, annotArr{iPlot}, ax, titleStrings{iPlot}, trialGroups);
    ax.FontSize = 11;
    ax.Title.FontSize = 11;
    ax.XLabel.FontSize = 13;
    if iPlot ~= nPlots
        ax.XLabel = [];
        ax.XTickLabels = [];
    end 
    
end

if saveFig
    % Create analysis directory if necessary
    saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    fileName = expDate;
    for iPlot = 1:nPlots
        fileName = [fileName, '_', regexprep(plotNames{iPlot}, ' ', '')];
    end
    fileName = [fileName, '_Summary ', fileNameSuffix];
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        savefig(f, fullfile(saveDir, fileName));
        export_fig(fullfile(saveDir, fileName), '-png', f);

    end
end%if

end%iFold

%% PLOT 1-D VISUALIZATION OF BEHAVIOR DATA ANNOTATIONS

%----- Plot 1D trial-averaged movement data -----
s = myData.stimSepTrials;

saveFig = 1;
fileNameSuffix = '_AllTrials_Locomotion';%'_OdorAvsOdorBvsAirPulse_Locomotion'; %
actionLabel = [ 2  ]; % locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
trialGroups = []; %[[s.OdorA + 2 * s.OdorB + 3 * s.OdorBPulse] .* goodTrials]; %  
figTitle = 'Fly locomotion throughout trial (green = odor)';
plotNames = {'ACV', 'EmptyVial', 'AirPulses'};

stimShadingColors = {'red', 'green'};

for iFold = 1
    
% Get approximate actual ball stopping times
ballStoppedSum = sum(bStopAnnotArr ./ 4);
ballStoppedFrames = ballStoppedSum > (nTrials * 0.25);
startFrames = regexp(num2str(ballStoppedFrames, '%d'), '01');
endFrames = regexp(num2str(ballStoppedFrames, '%d'), '10');
bStopFrames = [startFrames' endFrames'];

% Get odor stim times
odorTimes = [myData.odorStartTimes', myData.odorEndTimes'];
odorFrames = floor(odorTimes * FRAME_RATE);

stimShading = {bStopFrames, odorFrames};

% Create array of annotation data
behavAnnotations = behaviorAnnotArr; % --> [trial, frame]

f = figure(2); clf;
f.Position = [100 100 1600 500];
f.Color = [1 1 1];

if isempty(trialGroups)
    
    % Plot summed movement data
    annotArrSum = sum(ismember(behavAnnotations, actionLabel), 1);
    ax = gca();
    plot_behavior_summary_1D(myData, annotArrSum(2:end-1), ax, figTitle);
    
    % Add shading during stimulus presentations
    yL = ylim();
    for iType = 1:numel(stimShading)
        for iStim = 1:size(stimShading{iType}, 1)
            stimStart = stimShading{iType}(iStim, 1);
            stimLength = stimShading{iType}(iStim, 2) - stimShading{iType}(iStim, 1);
            rectPos = [stimStart, yL(1), stimLength, diff(yL)]; % [x y width height]
            rectangle('Position', rectPos, 'FaceColor', [rgb(stimShadingColors{iType}), 0.1], 'EdgeColor', 'none');
            ylim(yL);
        end
    end
else
    annotArrSum = [];
    yLimsAll = [];
    ax = [];
    for iGroup = 1:length(unique(trialGroups(trialGroups ~= 0)))
        
        % Plot summed movement data
        f.Position = [100 50 1000 950];
        ax{iGroup} = subplot(3, 1, iGroup);
        annotArrSum = sum(ismember(behavAnnotations(trialGroups == iGroup, :), actionLabel), 1);
        plot_behavior_summary_1D(myData, annotArrSum, ax{iGroup}, plotNames{iGroup});
        
        if iGroup ~= length(unique(trialGroups))
            xlabel('');
        end
        
        % Add shading during stimulus presentations
        yL = ylim();
        yLimsAll(iGroup, :) = yL;
        for iType = 1:numel(stimShading)
            for iStim = 1:size(stimShading{iType}, 1)
                stimStart = stimShading{iType}(iStim, 1);
                stimLength = stimShading{iType}(iStim, 2) - stimShading{iType}(iStim, 1);
                rectPos = [stimStart, yL(1), stimLength, diff(yL)]; % [x y width height]
                rectangle('Position', rectPos, 'FaceColor', [rgb(stimShadingColors{iType}), 0.1], 'EdgeColor', 'none');
                ylim(yL);
            end
        end
    end%iGroup
    
    % Make sure all plots use the same yLims
    yLimMax = max(yLimsAll(:));
    for iGroup = 1:length(unique(trialGroups(trialGroups~=0)))
        ylim(ax{iGroup}, [yLimsAll(iGroup, 1), yLimMax]);
    end
    
    suptitle(figTitle);
end

if saveFig
    % Create analysis directory if necessary
    saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    fileName = expDate;
    fileName = [fileName, '_Summed_Movement', fileNameSuffix];
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        savefig(f, fullfile(saveDir, fileName));
        export_fig(fullfile(saveDir, fileName), '-png', f);

    end
end%if

end%iFold


%% SEPARATE TRIALS BASED ON EVENT INTERACTIONS
analysisWindows = []; overshoots = [];


% Display available annotation event types

disp(annotationTypeSummary)

activeEventTypes = [1 7];
disp(annotationTypeSummary(activeEventTypes, 2))

% Odor
analysisWindows(end+1,:) = [ 1  1 ];
overshoots(end+1)        = 0;

% All behavior
analysisWindows(end+1,:) = [ 1  1 ];
overshoots(end+1)        = 0;

% % Ball stopping
% analysisWindows(end+1,:) = [ 1.5 1.5 ];
% overshoots(end+1)        = 0;

% Display active filter types
activeFilterTypes = [1 7];
disp(annotationTypeSummary(activeFilterTypes, 2))

% Create filters for different condition components
allFilts = []; allFiltNames = []; filtWindows = [];

% Odor
withOdor =  [ 0  1  0 ];
noOdor =    [-1 -1  0 ];
anyOdor =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  1   1  ];
allFilts{end+1} = [withOdor; noOdor; anyOdor];
allFiltNames{end+1} = {'WithOdor', 'NoOdor', 'AnyOdor'};

% All behavior
startMove = [-1  1  0 ];
endMove =   [ 1  0 -1 ];
contMove =  [ 1  1  0 ];
noMove =    [-1 -1 -1 ];
anyMove =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  1  1  ];
allFilts{end+1} = [startMove; endMove; contMove; noMove; anyMove];
allFiltNames{end+1} = {'StartMove', 'EndMove', 'contMove', 'NoMove', 'AnyMove'};

% % Ball stopping
% bStopped =  [ 0  1  0 ];
% bRelease =  [ 1  0 -1 ];
% noBall =    [-1 -1 -1 ];
% anyBall =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  2   2  ];
% allFilts{end+1} = [bStopped; bRelease; noBall; anyBall];
% allFiltNames{end+1} = {'BallStopped', 'BallRelease', 'NoBall', 'AnyBall'};

for iFold = 1

%===================================================================================================
%
%                                                 |------------event------------|          
%       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
%       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
%
%===================================================================================================
    
% Get event vols for each filter type
filterEventVols = [];
for iType = 1:numel(annotationTypes)
    filterEventVols(:,:,iType) = annotationTypes{iType}.eventVols;    
end

% Compile annotation info for active events
primaryEventNames = [];
for iType = 1:numel(activeEventTypes)
    currType = activeEventTypes(iType);
    eventLists{iType} = annotationTypes{currType}.eventList;
    primaryEventNames{iType} = annotationTypes{currType}.name;
end

nEventTypes = numel(primaryEventNames);
allCondFilters = []; allCondNames = []; activeFilterEventVols = []; onsetFilterVecs = []; offsetFilterVecs = [];
for iType = 1:nEventTypes
    
    
    % Identify any identical primary-secondary filter pairs
    matchedFilterTypes = (activeFilterTypes == activeEventTypes(iType)); 
    
    % Eliminate identical primary-secondary filter pairs
    currFilts = allFilts(~matchedFilterTypes);
    currFiltNames = allFiltNames(~matchedFilterTypes);
    
    % Combine primary events and filter types into all filter conditions
    [allCondFilters{iType}, allCondNames{iType}] =  create_filter_conditions(primaryEventNames{iType}, currFilts, currFiltNames);
    
    % Select the correct eventVols for each primary event type
    activeFilterEventVols = filterEventVols(:,:,activeFilterTypes(~matchedFilterTypes));    
    
    % Get onset filter vecs for each condition
    for iCond = 1:numel(allCondNames{iType})
        onsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, volumeRate, 'overshoot', overshoots(iType));
    end
    
    % Get offset filter vecs for each condition
    for iCond = 1:numel(allCondNames{iType})
        offsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, volumeRate, 'overshoot', overshoots(iType), 'offsetAlign', 1);
    end
    
end

% ----------------------------------------------------------------------------------------------
% Calculate dF/F
% ----------------------------------------------------------------------------------------------

onsetDffAvg = []; offsetDffAvg = []; onsetDff = []; offsetDff = []; combinedDff = []; combinedDffAvg =[];
onsetCondSummaries = []; offsetCondSummaries = []; allCondSummaries = [];
for iType = 1:nEventTypes
    
    primaryFiltName = primaryEventNames{iType};
    analysisWindow = analysisWindows(iType, :);
    nConds = numel(allCondNames{iType});
    eventList = eventLists{iType};
    
    baselineDur = analysisWindow(1);
    respDur = analysisWindow(2);
    onsetDff{iType} = zeros([sessionSize(1:3), (sec2vols(respDur, volumeRate) + sec2vols(baselineDur, volumeRate) + 2), nConds]);
    offsetDff{iType} = onsetDff{iType};
    onsetDffAvg{iType} = zeros([sessionSize(1:3), nConds]);
    offsetDffAvg{iType} = onsetDffAvg{iType};
    
    for iCond = 1:nConds
        
        disp(['Calculating dF/F for ', primaryFiltName, ' cond #', num2str(iCond), ' of ', num2str(nConds), '...'])
        
        % Calculate dF/F for event onsets
        if sum(onsetFilterVecs{iType}(:,iCond)) > 0

            [baselineData, respData] = extract_event_volumes(eventList, onsetFilterVecs{iType}(:,iCond), baselineDur, respDur, myData, ...
                'offsetAlign', 0); % --> [y, x, plane, volume, event]
            
            baselineDffAvg = mean(mean(baselineData, 5), 4);                            % --> [y, x, plane]
            baselineDffRep = repmat(baselineDffAvg, 1, 1, 1, size(baselineData, 4));    % --> [y, x, plane, volume]
            baselineDff = calc_dFF(baselineData, baselineDffRep, 5);                    % --> [y, x, plane, volume]
            
            currDff = calc_dFF(respData, baselineData, 5);                              % --> [y, x, plane, volume]
            combinedVolsDff = cat(4, baselineDff, currDff);                             % --> [y, x, plane, volume]
            currDffAvg = calc_dFF(respData, baselineData, [4 5]);                       % --> [y, x, plane]
            
            onsetDff{iType}(:,:,:,:, iCond) = combinedVolsDff;                          % --> {eventType}[y, x, plane, volume, condition]
            onsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                              % --> {eventType}[y, x, plane, condition]
            
        end
        
        % Calculate dF/F for event offsets
        if sum(offsetFilterVecs{iType}(:,iCond)) > 0
            
            [baselineData, respData] = extract_event_volumes(eventList, offsetFilterVecs{iType}(:,iCond), baselineDur, respDur, myData, ...
                'offsetAlign', 1); % --> [y, x, plane, volume, event]
            
            baselineDffAvg = mean(mean(baselineData, 5), 4);                            % --> [y, x, plane]
            baselineDffRep = repmat(baselineDffAvg, 1, 1, 1, size(baselineData, 4));    % --> [y, x, plane, volume]
            baselineDff = calc_dFF(baselineData, baselineDffRep, 5);                    % --> [y, x, plane, volume]
            
            currDff = calc_dFF(respData, baselineData, 5);                              % --> [y, x, plane, volume]
            combinedVolsDff = cat(4, baselineDff, currDff);                             % --> [y, x, plane, volume]
            currDffAvg = calc_dFF(respData, baselineData, [4 5]);                       % --> [y, x, plane]
            
            offsetDff{iType}(:,:,:,:, iCond) = combinedVolsDff;                         % --> {eventType}[y, x, plane, volume, condition]
            offsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                             % --> {eventType}[y, x, plane, condition]
        end
    end% iCond
    
    combinedDff{iType} = cat(5, onsetDff{iType}, offsetDff{iType});                     % --> {eventType}[y, x, plane, volume, condition]
    combinedDffAvg{iType} = cat(4, onsetDffAvg{iType}, offsetDffAvg{iType});            % --> {eventType}[y, x, plane, condition]
    
    % Create summary table for onset conditions
    condCountCol = (sum(onsetFilterVecs{iType})');
    alignCol = repmat({'onset'}, nConds, 1);
    baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
    respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
    exWinCol = repmat(filtWindows(iType, :), nConds, 1);
    rowNames = cellfun(@num2str, num2cell(1:nConds), 'uniformOutput', 0);
    varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp', 'ExcludeWin'};
    onsetCondSummaries{iType} = table(condCountCol, allCondNames{iType}, alignCol, baselineCol, respCol, exWinCol, 'RowNames', rowNames, 'VariableNames', varNames);
    
    % Create summary table for offset conditions
    condCountCol = (sum(offsetFilterVecs{iType})');
    alignCol = repmat({'offset'}, nConds, 1);
    baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
    respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
    exWinCol = repmat(filtWindows(iType, :), nConds, 1);
    rowNames = cellfun(@num2str, num2cell((1:nConds) + nConds), 'uniformOutput', 0);
    varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp', 'ExcludeWin'};
    offsetCondSummaries{iType} = table(condCountCol, allCondNames{iType}, alignCol, baselineCol, respCol, exWinCol, 'RowNames', rowNames, 'VariableNames', varNames);
    
    % Concatenate onset and offset summary tables
    allCondSummaries{iType} = vertcat(onsetCondSummaries{iType}, offsetCondSummaries{iType});
    disp(allCondSummaries{iType})
end% iType
disp('dF/F calculation complete')

end%iFold


%% =================================================================================================
%            ODOR STIM ANALYSES                                   
%%%=================================================================================================
for iFoldOut = 1
    
    %% PLOT ODOR ONSET/OFFSET HEATMAPS FOR SOME TRIAL CONDITIONS
    
    % Show summary again
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, 'odor'));
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
    
    condNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [4 5 10];
    sigma = [0.6];   
    rangeType = 'Max';
    rangeScalar = 0.15;
    makeVid = 0;
    saveDir = [];
    fileName = 'Odor_Response_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(condNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir);
     
end%iFoldOut
    
    %% PLOT MEAN dF/F RESPONSES WITHIN ROIs FOR SOME TRIAL CONDITIONS
    
%% =================================================================================================
%           BEHAVIOR ANALYSES                                   
%%%=================================================================================================
for iFold = 1
    %% CALCULATE AND PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES

    locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
    actionLabel = [2];
    baselineLabel = [0];

    smoothingSigma = [0.6]; 
    rangeType = 'max';
    rangeScalar = 0.1;
    makeVid = 0;
    saveDir = [];
    fileName = 'Locomotion_Plane_Heatmaps';
    titleStr = {'dF/F - Locomotion vs. Quiescence'};
    
for iFoldIn = 1
    % Identify behavioral state during each volume
    actionVols = zeros(myData.nTrials, myData.nVolumes); stoppedVols = actionVols;
    for iTrial = 1:myData.nTrials
        if myData.goodTrials(iTrial)

            % Pull out action numbers for each volume
            currActions = myData.trialAnnotations{iTrial}.actionNums;
            volActions = currActions(volFrames);

            % Identify volume actions
            actionVols(iTrial, :) =  ismember(volActions, actionLabel);   %--> [trial, vol]
            stoppedVols(iTrial, :) = ismember(volActions, baselineLabel); %--> [trial, vol]
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
    imgData = myData.wholeSession; %--> [y, x, plane, volume, trial]
    for iTrial = 1:myData.nTrials
       currImgData = imgData(:,:,:,:,iTrial);
       currActionVols = logical(actionVols(iTrial,:));
       currStoppedVols = logical(stoppedVols(iTrial,:));
       % Exclude volumes that occurred during or just after the wind stim
       if myData.stimSepTrials.windTrials(iTrial)
           stimStartVol = ceil(stimStart * volumeRate);
           stimEndVol = floor(stimEnd * volumeRate);
           excludeVols = stimStartVol:(stimEndVol + ceil(volumeRate));
           currImgData(:,:,:,excludeVols) = [];
           currActionVols(excludeVols) = [];
           currStoppedVols(excludeVols) = [];
       end
       
       % Pull out running volumes, if any exist, from the current trial
       if sum(currActionVols) > 0
           meanActionVols(:,:,:,end+1) = mean(currImgData(:,:,:,currActionVols),4);    %--> [y, x, plane, trial]
       end
       
       % Pull out stopping volumes, if any exist, from the current trial
       if sum(currStoppedVols) > 0
           meanStoppedVols(:,:,:,end+1) = mean(currImgData(:,:,:,currStoppedVols),4);  %--> [y, x, plane, trial]
       end

    end
    actionMean = mean(meanActionVols, 4);   %--> [y, x, plane]
    stoppedMean = mean(meanStoppedVols, 4); %--> [y, x, plane] 

    % Get dF/F values for action relative to quiescence
    actionDff = (actionMean - stoppedMean) ./ stoppedMean; % --> [y, x, plane]

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(actionDff, rangeScalar, rangeType);

    % Plot figures
    [f, ~] = plot_heatmaps(actionDff, myData, range, titleStr, smoothingSigma, 'fileName', fileName, 'makeVid', makeVid, ...
                           'saveDir', saveDir);

end%iFoldIn

    %% PLOT INTERACTION HEATMAPS FOR SOME TRIAL CONDITIONS
    
    
    % Show summary again
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, 'move'));
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
                 
    condNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [2 5];
    sigma = [0.6]; 
    rangeType = 'Max';
    rangeScalar = 0.5;
    makeVid = 0;
    saveDir = [];
    fileName = 'Behavior_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(condNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 


            %% CREATE VIDEO OF MEAN dF/F FOR EACH PLANE THROUGHOUT MOVEMENT ONSET
        
            rangeScalar = [0.75];
            sigma = 0.75;
            offsetAlign = 0;

for iFold = 1
    
            % Select correct alignment
            if offsetAlign
                combinedVolsDff = combinedVolsOffsetDff;
                baselineDff = offsetBaselineDff;
                respDff = behavOffsetDff;
                titleStr = 'offset';
                 fileName = ['Behavior_Offset_Response_Heatmap_Vid_', num2str(baselineDur), '_', num2str(respDur)];
            else
                combinedVolsDff = combinedVolsOnsetDff;
                baselineDff = onsetBaselineDff;
                respDff = behavOnsetDff;
                titleStr = 'onset';
                fileName = ['Behavior_Onset_Response_Heatmap_Vid', num2str(baselineDur), '_', num2str(respDur)];
            end

            % Calculate volume times in seconds relative to odor onset
            baselineVolTimes = -(1:size(onsetBaselineDff, 4)) / volumeRate;
            respVolTimes = (1:size(behavOnsetDff, 4)) / volumeRate;
            relTimes = [baselineVolTimes(end:-1:1), respVolTimes];

            % Create cell array with titles for each frame
            titleStrings = [];
            for iVol = 1:size(combinedVolsOnsetDff, 4)
                if iVol <= size(onsetBaselineDff, 4)
                    titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before behavior ', titleStr];
                else
                    titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after behavior ', titleStr];
                end
            end

            % Create video
            range = calc_range(combinedVolsDff, rangeScalar);
            savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
            make_heatmap_vid(combinedVolsDff, myData, range, fileName, titleStrings, savePath, [], [], [], sigma);  
end%iFold            

 
end%iFold

%% =================================================================================================
%           BALL STOPPING ANALYSES                                   
%%%=================================================================================================
for iFold = 1           

    %% PLOT BALL STOP/RELEASE HEATMAPS FOR SOME TRIAL CONDITIONS
    
    % Show summary again
    disp(bStopCondSummary)
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [1 2 3 4];
    sigma = [0.6]; 
    rangeType = 'max';
    rangeScalar = .002;
    makeVid = 1;
    saveDir = [];
    fileName = 'Ball_Stopping_Interaction_Heatmaps_No_Filter';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [bStopCondNames{iCond}, '  (n = ', num2str(bStopCondSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = bStopCondDffAvg(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 

end%iFold

%% =================================================================================================
%           OTHER ANALYSES                                   
%%%=================================================================================================
for iFold = 1

%% CREATE VIDEO OF MEAN dF/F THROUGHOUT ENTIRE TRIAL

%++++++++++ Calculate whole-trial dF/F using overall mean as baseline ++++++++++

allWindTrials = myData.wholeSession;%(:,:,:,:,myData.stimSepTrials.windTrials);   % --> [y, x, plane, volume, trial]   

% Calculate dF/F using whole trial average as baseline
baseline = mean(mean(allWindTrials, 5), 4);                           % --> [y, x, plane]
baselineMeanRep = repmat(baseline, 1, 1, 1, size(allWindTrials, 4));  % --> [y, x, plane, volume]
trialAvg = mean(allWindTrials, 5);                                    % --> [y, x, plane, volume]
wholeTrialDff = (trialAvg - baselineMeanRep) ./ baselineMeanRep;      % --> [y, x, plane, volume]

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

windTrials = myData.wholeSession(:,:,:,:,~logical(myData.stimSepTrials.windTrials));   % --> [y, x, plane, volume, trial]


%% PCA PLOTTING

% Pull out data for one plane
planeNum = 11;
planeData = squeeze(windTrials(:,:,planeNum,:,:)); % --> [y, x, volume, trial]
[w,x,y,z] = size(planeData);
planeDataRS = reshape(planeData, [w, x, y*z]); % --> [y, x, volume]


pcaData = mean(squeeze(myData.wholeSession(:,:,planeNum,:,:)),4); % --> [y, x, volume]
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