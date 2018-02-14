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
    myData.MAX_INTENSITY = 2000; MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots
    myData.FRAME_RATE = 25; FRAME_RATE = 25; % This is the frame rate of the behavior video, not the GCaMP imaging
    if isempty(nFrames)
        nFrames = sum(trialDuration) * FRAME_RATE;
    end
    volTimes = (1:nVolumes)' ./ volumeRate;
    frameTimes = (1:nFrames)' ./ FRAME_RATE;
    
end%iFold

%% VIEW RAW DATA FOR A SINGLE TRIAL AND PLANE

planeNum = 10;
trialNum = 2; % Does not account for any skipped trials

preview_trial_movie(myData.wholeSession, planeNum, trialNum, [], [], []);

%% INITIAL DATA PROCESSING STEPS

skipTrials = [28 52]; 
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


%% LOAD ROI DATA

parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sid, '\Analysis'];
fileName = 'ROI_Data.mat';

load(fullfile(parentDir, fileName));
myData.ROIdata = ROIdata;


%% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS

saveFig = 0;
plotTypes = [2 1]; % 1 = odor stims, 2 = behavior, 3 = ball stopping
s = myData.stimSepTrials;
trialGroups = [[s.OdorA + 2 * s.OdorB] .* goodTrials];%[s.OdorA + 2 * s.OdorB + 3 * s.NoOdor] .* goodTrials; %[s.odorTrials + 2 * (~s.odorTrials)]; % 
plotTitleSuffix = '(ACV vs Control)';%''%
fileNameSuffix = '_OdorAvsOdorB';%'_AllTrials';%

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
fileNameSuffix = '_OdorAvsOdorB_Locomotion'; %'_AllTrials_Locomotion';%
actionLabel = [ 2  ]; % locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
trialGroups = [[s.OdorA + 2 * s.OdorB] .* goodTrials]; % 
figTitle = 'Fly locomotion throughout trial (green = odor)';
plotNames = {'ACV', 'EmptyVial'};

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
    for iGroup = 1:length(unique(trialGroups(trialGroups ~= 0)))
        
        % Plot summed movement data
        f.Position = [100 50 1000 950];
        ax = subplot(3, 1, iGroup);
        annotArrSum = sum(ismember(behavAnnotations(trialGroups == iGroup, :), actionLabel), 1);
        plot_behavior_summary_1D(myData, annotArrSum, ax, plotNames{iGroup});
        
        if iGroup ~= length(unique(trialGroups))
            xlabel('');
        end
        
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
overshoots(end+1)        = 1;

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
filtWindows(end+1,:)     = [  1   1  ];
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
    
    combinedDff{iType} = cat(5, onsetDff{iType}, offsetDff{iType});
    combinedDffAvg{iType} = cat(4, onsetDffAvg{iType}, offsetDffAvg{iType});
    
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
    currConds = [6 8 9];
    sigma = [0.6];   
    rangeType = 'max';
    rangeScalar = 0.6;
    makeVid = 0;
    saveDir = [];
    fileName = 'Odor_Onset_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(condNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 
    
    %% CALCULATE dF/F FOR ODOR STIM ONSET AND OFFSET 
    
    baselineDur = 1;
    respDur = 1;
    eventExclusionWindow = [1 1];
    overshoot = 1;
    
    % Select events
    eventList = odorEventList;
    odorFilter      = [0 0 0];
    odorWindow      = [1 1];
    ballStopFilter  = [0 0 0];
    ballStopWindow  = eventExclusionWindow;
    behavFilter     = [-1 -1 -1];
    behavWindow     = eventExclusionWindow;
    
for iFold = 1
    filterEventVols = cat(3, odorVols, ballStoppingVols, behavVols);
    filterEventWindows = [odorWindow; ballStopWindow; behavWindow];
    filterDirections = [odorFilter; ballStopFilter; behavFilter];
    analysisWindow = [2 1];
    
    filterVec = filter_event_data(eventList, filterEventVols, analysisWindow, filterEventWindows, filterDirections, volumeRate);
    
    [onsetBaselineData, onsetRespData] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData);                     % --> [y, x, plane, volume, event]
    [offsetBaselineData, offsetRespData] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData, 'offsetAlign', 1); % --> [y, x, plane, volume, event]
    
    odorOnsetDff = calc_dFF(onsetRespData, onsetBaselineData, 5);                               % --> [y, x, plane, volume]
    onsetBaselineDffAvg = mean(mean(onsetBaselineData, 5), 4);                                  % --> [y, x, plane]
    onsetBaselineDffRep = repmat(onsetBaselineDffAvg, 1, 1, 1, size(onsetBaselineData, 4));     % --> [y, x, plane, volume]
    onsetBaselineDff = calc_dFF(onsetBaselineData, onsetBaselineDffRep, 5);                     % --> [y, x, plane, volume]
    
    odorOffsetDff = calc_dFF(offsetRespData, offsetBaselineData, 5);                            % --> [y, x, plane, volume]
    offsetBaselineDffAvg = mean(mean(offsetBaselineData, 5), 4);                                % --> [y, x, plane]
    offsetBaselineDffRep = repmat(offsetBaselineDffAvg, 1, 1, 1, size(offsetBaselineData, 4));  % --> [y, x, plane, volume]
    offsetBaselineDff = calc_dFF(offsetBaselineData, offsetBaselineDffRep, 5);                  % --> [y, x, plane, volume]
        
    combinedVolsOnsetDff = cat(4, onsetBaselineDff, odorOnsetDff);
    combinedVolsOffsetDff = cat(4, onsetBaselineDff, odorOffsetDff);

    odorOnsetDffAvg = calc_dFF(onsetRespData, onsetBaselineData, [4 5]);            % --> [y, x, plane]
    odorOffsetDffAvg = calc_dFF(offsetRespData, offsetBaselineData, [4 5]);         % --> [y, x, plane]

    % Get dF/F averages for individual odors as well
    eventList = odorAEventList;
    filterVec = filter_event_data(eventList, filterEventVols, analysisWindow, filterEventWindows, filterDirections, volumeRate);
    [onsetBaselineData_A, onsetRespData_A] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData);                     % --> [y, x, plane, volume, event]
    [offsetBaselineData_A, offsetRespData_A] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData, 'offsetAlign', 1); % --> [y, x, plane, volume, event]
    odorOnsetDffAvg_A = calc_dFF(onsetRespData_A, onsetBaselineData_A, [4 5]);      % --> [y, x, plane]
    odorOffsetDffAvg_A = calc_dFF(offsetRespData_A, offsetBaselineData_A, [4 5]);   % --> [y, x, plane]
    
    eventList = odorBEventList;
    filterVec = filter_event_data(eventList, filterEventVols, analysisWindow, filterEventWindows, filterDirections, volumeRate);
    [onsetBaselineData_B, onsetRespData_B] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData);                     % --> [y, x, plane, volume, event]
    [offsetBaselineData_B, offsetRespData_B] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData, 'offsetAlign', 1); % --> [y, x, plane, volume, event]
    odorOnsetDffAvg_B = calc_dFF(onsetRespData_B, onsetBaselineData_B, [4 5]);      % --> [y, x, plane]
    odorOffsetDffAvg_B = calc_dFF(offsetRespData_B, offsetBaselineData_B, [4 5]);   % --> [y, x, plane]
end


%% PLOT ODOR ONSET AND OFFSET dF/F HEATMAPS
    
    rangeType = 'max';
    rangeScalar = 1;
    sigma = 0.6;
    plotTitles = {'Odor stim onset dF/F', 'Odor stim offset dF/F'};
    makeVid = 0;
    fileName = ['Odor_Only_Heatmaps_', num2str(baselineDur), '_', num2str(respDur)];
    
    plotData = cat(4, odorOnsetDffAvg, odorOffsetDffAvg);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    plot_heatmaps(plotData, myData, range, plotTitles, sigma, 'makeVid', makeVid, 'fileName', fileName);
    
    %% COMPARE ODOR A vs B ONSET/OFFSET RESPONSES
    
    rangeType = 'max';
    rangeScalar = 0.5;
    sigma = 0.6;
    plotTitles = {'Odor A onset dF/F', 'Odor B onset dF/F', 'Odor A stim offset dF/F', 'Odor B stim offset dF/F'};
    makeVid = 0;
    fileName = ['Odor_Only_AvsB_Heatmaps_', num2str(baselineDur), '_', num2str(respDur)];
    
    plotData = cat(4, odorOnsetDffAvg_A, odorOffsetDffAvg_A, odorOnsetDffAvg_B, odorOffsetDffAvg_B);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    plot_heatmaps(plotData, myData, range, plotTitles, sigma, 'makeVid', makeVid, 'fileName', fileName);
    
            %% CREATE VIDEO OF MEAN dF/F THROUGHOUT ALL ODOR RESPONSES
    
    rangeType = 'stdDev';
    rangeScalar = [0.6];
    sigma = 0.6;
    offsetAlign = 0;
    
 for iFold = 1
    % Select correct alignment
    if offsetAlign
        combinedVolsDff = combinedVolsOffsetDff;
        baselineDff = offsetBaselineDff;
        respDff = odorOffsetDff;
        titleStr = 'offset';
         fileName = ['Odor_Offset_Response_Heatmap_Vid_', num2str(baselineDur), '_', num2str(respDur)];
    else
        combinedVolsDff = combinedVolsOnsetDff;
        baselineDff = onsetBaselineDff;
        respDff = odorOnsetDff;
        titleStr = 'onset';
        fileName = ['Odor_Onset_Response_Heatmap_Vid', num2str(baselineDur), '_', num2str(respDur)];
    end
    
    % Calculate volume times in seconds relative to odor onset
    baselineVolTimes = -(1:size(onsetBaselineDff, 4)) / volumeRate;
    respVolTimes = (1:size(odorOnsetDff, 4)) / volumeRate;
    relTimes = [baselineVolTimes(end:-1:1), respVolTimes];
    
    % Create cell array with titles for each frame
    titleStrings = [];
    for iVol = 1:size(combinedVolsOnsetDff, 4)
        if iVol <= size(onsetBaselineDff, 4)
            titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before odor ', titleStr];
        else
            titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after odor ', titleStr];
        end
    end
    
    % Create video
    range = calc_range(combinedVolsDff, rangeScalar, rangeType);
    savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    make_heatmap_vid(combinedVolsDff, myData, range, fileName, titleStrings, savePath, [], [], [], sigma);
 end
 
end%iFoldOut

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
    rangeScalar = 1;
    makeVid = 1;
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
    currConds = [1 2 5];
    sigma = [0.6]; 
    rangeType = 'Max';
    rangeScalar = 0.8;
    makeVid = 1;
    saveDir = [];
    fileName = 'Behavior_Interaction_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(condNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 
   

            %% PLOT BEHAVIOR ONSET/OFFSET HEATMAPS FOR EACH PLANE
    
            rangeType = 'max';
            rangeScalar = 1;
            sigma = 0.6;
            plotTitles = {'Behavior onset dF/F', 'Behavior offset dF/F'};
            makeVid = 0;
            fileName = ['Behavior_Onset_Offset_Heatmaps_No_Odor', num2str(baselineDur), '_', num2str(respDur)];

            plotData = cat(4, behavOnsetDffAvg, behavOffsetDffAvg);
            range = calc_range(plotData, rangeScalar, rangeType);
            plot_heatmaps(plotData, myData, range, plotTitles, sigma, 'makeVid', makeVid, 'fileName', fileName);
    

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

    %% CREATE COMBINED PLOTTING VIDEO OF dF/F FOR ALL INDIVIDUAL MOVEMENT BOUT ONSETS

    baselineLength = 12; % Number of baseline volumes to analyze/plot for each movement bout
    rangeScalar = [0.5];
    
for iFoldIn = 1

    % Select behavior events
    eventExclusionWindow = [2 2];
    eventList = behavEventList;
    odorFilter      = [0 0 0];
    odorWindow      = eventExclusionWindow;
    ballStopFilter  = [0 0 0];
    ballStopWindow  = eventExclusionWindow;
    behavFilter     = [-1 0 0];
    behavWindow     = eventExclusionWindow;
                   
    filterEventVols = cat(3, odorVols, ballStoppingVols, behavVols);
    filterEventWindows = [odorWindow; ballStopWindow; behavWindow];
    filterDirections = [odorFilter; ballStopFilter; behavFilter];
    
    filterVec = filter_event_data(eventList, filterEventVols, filterEventWindows, filterDirections);
    allBouts = behavEventList(filterVec, :);

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

        % Calculate dF/F for the current bout ****Is this the best way to calculate dF/F?*****
        startInd = allBouts(iBout, 1) - baselineLength - 1;
        if startInd < 1
            baselineLength = baselineLength + startInd;
            startInd = 1;
        end
        endInd = allBouts(iBout, 2);
        fullLength = endInd - startInd;
        boutData = double(squeeze(myData.wholeSession(:,:,:, (startInd + 1):endInd, allBouts(iBout,3))));   % --> [y, x, plane, volume]
        boutBaseline = mean(boutData(:,:,:,1:baselineLength),4);                                      % --> [y, x, plane]
        boutBaselineRep = repmat(boutBaseline,1,1,1,fullLength);     % --> [y, x, plane, volume]
        boutDff = (boutData - boutBaselineRep) ./ boutBaselineRep;                                    % --> [y, x, plane, volume]

        boutDff(isinf(boutDff)) = 0; % To eliminate inf values from dividiing by zero above
        boutDff(isnan(boutDff)) = 0;

        % Average each volume with its neighbors to improve SNR
        boutDff = movmean(boutDff, 3, 3);

        % Calculate dF/F value range
        range = calc_range(boutDff, rangeScalar);

        % Load behavior vid for the current trial
        vidDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', myData.expDate, '_Movies');
        trialStr = ['sid_', num2str(myData.sid), '_tid_', sprintf('%03.f',(allBouts(iBout, 3)))];
        myMovie = [];
        flyVid = VideoReader(fullfile(vidDir, [trialStr, '.mp4']));
        while hasFrame(flyVid)
            currFrame = readFrame(flyVid);
            myMovie(:,:,end+1) = rgb2gray(currFrame);
        end
        myMovie = uint8(myMovie(:,:,2:end)); % Adds a black first frame for some reason, so drop that
        
        % Calculate frame times in seconds relative to behavior bout onset
        boutFrames = volFrames(allBouts(iBout,1)):volFrames(allBouts(iBout,2));
        baselineLengthFrames = volFrames(baselineLength);
        plotFrames = (volFrames(allBouts(iBout, 1)) - baselineLengthFrames):volFrames(allBouts(iBout, 2));
        baselineFrameTimes = frameTimes(baselineLengthFrames:-1:1);
        boutFrameTimes = frameTimes(1:(numel(boutFrames)));
        relFrameTimes = [baselineFrameTimes; boutFrameTimes];
        for iFrame = 1:numel(plotFrames)
            if iFrame <= baselineLengthFrames
                relTimeStr = 'before';
            else
                relTimeStr = 'after';
            end
            titleStrings{iFrame} = ['Time = ', sprintf('%05.2f', relFrameTimes(iFrame)), ' sec ', relTimeStr, ' movement onset'];
        end

        for iFrame = 1:numel(plotFrames)

            % Create fig
            f = figure(1); clf
            f.Position = [50 45, 1620, 855];

            for iPlane = myData.nPlanes:-1:1 % Planes going from dorsal --> ventral

                % Plot dF/F for each plane
                ax = subaxis(3, 4, iPlane, 'Spacing', 0, 'MB', 0.025);

                if iPlane == 1
                    imshow(myMovie(:,:, plotFrames(iFrame)));
                    axis image; axis off
                else
                    [~, currVol] = min(abs(volFrames - iFrame));
                    imagesc(boutDff(:, :, iPlane, currVol))
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
    %% CALCULATE MEAN dF/F AROUND BALL STOPPING AND PLOT HEATMAPS

    baselineDur = 2;
    respDur = 2;
    eventExclusionWindow = [2 2];
    overshoot = 0;
    
    rangeType = 'stdDev';
    rangeScalar = 1;
    sigma = 0.6;
    plotTitles = {'Ball stopping dF/F', 'Ball release dF/F', 'Ball release dF/F (pre-stop baseline)'};
    makeVid = 0;
    fileName = ['Ball_Stopping_Heatmaps', num2str(baselineDur), '_', num2str(respDur)];
    
    % Select events 
    eventList = ballStoppingEventList;
    odorFilter      = [0 0 0];
    odorWindow      = eventExclusionWindow;
    ballStopFilter  = [0 0 0];
    ballStopWindow  = eventExclusionWindow;
    behavFilter     = [0 0 0];
    behavWindow     = eventExclusionWindow;
    
for iFold = 1
    filterEventVols = cat(3, odorVols, ballStoppingVols, behavVols);
    filterEventWindows = [odorWindow; ballStopWindow; behavWindow];
    filterDirections = [odorFilter; ballStopFilter; behavFilter];
    
    filterVec = filter_event_data(eventList, filterEventVols, filterEventWindows, filterDirections);
    
    % Append number of events to plot titles
    titleStr = ['nTrials = ', num2str(sum(filterVec))];
    disp(titleStr)
    for iPlot = 1:numel(plotTitles)
        plotTitles{iPlot} = {[plotTitles{iPlot}, ' (', titleStr, ')']};
    end
    
    [onsetBaselineData, onsetRespData] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData, 'overshoot', overshoot);                     % --> [y, x, plane, volume, event]
    [offsetBaselineData, offsetRespData] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, myData, 'offsetAlign', 1, 'overshoot', overshoot); % --> [y, x, plane, volume, event]
    
    ballStoppingDff = calc_dFF(onsetRespData, onsetBaselineData, 5);                            % --> [y, x, plane, volume]
    onsetBaselineDffAvg = mean(mean(onsetBaselineData, 5), 4);                                  % --> [y, x, plane]
    onsetBaselineDffRep = repmat(onsetBaselineDffAvg, 1, 1, 1, size(onsetBaselineData, 4));     % --> [y, x, plane, volume]
    onsetBaselineDff = calc_dFF(onsetBaselineData, onsetBaselineDffRep, 5);                     % --> [y, x, plane, volume]
    
    ballReleaseDff = calc_dFF(offsetRespData, offsetBaselineData, 5);                           % --> [y, x, plane, volume]
    offsetBaselineDffAvg = mean(mean(offsetBaselineData, 5), 4);                                % --> [y, x, plane]
    offsetBaselineDffRep = repmat(offsetBaselineDffAvg, 1, 1, 1, size(offsetBaselineData, 4));  % --> [y, x, plane, volume]
    offsetBaselineDff = calc_dFF(offsetBaselineData, offsetBaselineDffRep, 5);                  % --> [y, x, plane, volume]
        
    combinedVolsStoppingDff = cat(4, onsetBaselineDff, ballStoppingDff);
    combinedVolsReleaseDff = cat(4, onsetBaselineDff, ballReleaseDff);
    
    
    ballStoppingDffAvg = calc_dFF(onsetRespData, onsetBaselineData, [4 5]);         % --> [y, x, plane]
    ballReleaseDffAvg = calc_dFF(offsetRespData, offsetBaselineData, [4 5]);        % --> [y, x, plane]
    ballReleaseRelDffAvg = calc_dFF(offsetRespData, onsetBaselineData, [4 5]);      % --> [y, x, plane]
    
    plotData = cat(4, ballStoppingDffAvg, ballReleaseDffAvg, ballReleaseRelDffAvg); % --> [y, x, plane, plot]
    range = calc_range(plotData(:,:,2:end-1,:), rangeScalar, rangeType);
    plot_heatmaps(plotData, myData, range, plotTitles, sigma, 'makeVid', makeVid, 'fileName', fileName);    
end
            

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
planeNum = 6;
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
%     %---------- Get trial averaged baseline and stimulus data ----------                                                                                              % myData.wholeSession = [y, x, plane, volume, trial]                                                                                                          
%     wholeTrialAvg(:,:,:,:,iStim) = mean(myData.wholeSession(:,:,:,:,stimSepTrials.(stimTypes{iStim})), 5);                                                            % --> [y, x, plane, volume, StimType]
%     
%     % Pre-stim
%     if floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate) > 0 
%         baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimStart*volumeRate)-floor((stimEnd-stimStart)*volumeRate):floor(stimStart*volumeRate),iStim), 4); % --> [y, x, plane, StimType]
%     else
%         % if stimDuration > preStimDuration, start baseline one second after beginning of trial
%         baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),iStim), 4);                                                 % --> [y, x, plane, StimType]
%     end
%     
%     % Stim period
%     stimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,ceil(stimStart*volumeRate):floor(stimEnd*volumeRate),iStim), 4);                                                  % --> [y, x, plane, StimType]
%     
%     % Post stim
%     if floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate) < size(wholeTrialAvg, 4)
%         postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):floor(stimEnd*volumeRate) + floor((stimEnd-stimStart)*volumeRate),iStim), 4);   % --> [y, x, plane, StimType]
%     else
%         % if stimEnd + stimLength > fullTrialDuration, end post-stim period at end of trial
%         postStimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:,floor(stimEnd*volumeRate):size(wholeTrialAvg, 4)-floor(volumeRate),iStim), 4);                            % --> [y, x, plane, StimType]
%     end
%     
% end%for
% 
% % Calculate dF/F values
% dffAvg = (stimAvg - baselineAvg) ./ baselineAvg; % --> [y, x, plane, StimType]
% dffAvgPost = (postStimAvg - stimAvg) ./ stimAvg; % --> [y, x, plane, StimType]
% 
% %% CALCULATE MEAN dF/F AROUND WIND RESPONSES
% stimSepTrials = [];
% 
% % Separate out  wind trials
% for iStim = 1:length(stimTypes)
%     stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));  
% end 
% windTrials = myData.wholeSession(:,:,:,:,logical(myData.stimSepTrials.windTrials));  % --> [y, x, plane, volume, trial]
% 
% % Calculate dF/F before and after wind onset using an equal period before onset as baseline
% stimLength = stimEnd - stimStart;
% stimLengthVols = floor(stimLength * volumeRate);
% 
% % Pre-stim
% if floor(stimStart*volumeRate) - stimLengthVols > 0
%     baselineVols = windTrials(:,:,:,floor(stimStart*volumeRate) - stimLengthVols:floor(stimStart*volumeRate),:);                                   % --> [y, x, plane, volume, trial]
% else
%     % If stimDuration > preStimDuration, start baseline one second after beginning of trial
%     baselineVols = windTrials(:,:,:,floor(volumeRate):floor(stimStart*volumeRate),:);                                                              % --> [y, x, plane, volume, trial]
% end
% 
% % Stim period + 3 seconds after the stim offset
% stimVols = windTrials(:,:,:,ceil((stimStart*volumeRate)-size(baselineVols, 4)):floor((stimStart*volumeRate) + (3*volumeRate) + stimLengthVols),:); % --> [y, x, plane, volume, trial]  
% 
% stimVolsMean = mean(stimVols, 5);                                                                                                                  % --> [y, x, plane, volume]
% baselineMean = mean(mean(baselineVols, 5), 4);                                                                                                     % --> [y, x, plane]
% baselineMeanRep = repmat(baselineMean, 1, 1, 1, size(stimVols, 4));                                                                                % --> [y, x, plane, volume]
% windDffVols = (stimVolsMean - baselineMeanRep) ./ baselineMeanRep;                                                                                 % --> [y, x, plane, volume]
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
%     boutData = double(squeeze(myData.wholeSession(:,:,:, allBouts(iBout,2):allBouts(iBout,3),allBouts(iBout,1))));   % --> [y, x, plane, volume]
%     boutBaseline = mean(boutData(:,:,:,1:baselineLength),4);                                                         % --> [y, x, plane]
%     boutBaselineRep = repmat(boutBaseline,1,1,1,length(allBouts(iBout,2):allBouts(iBout,3)));                        % --> [y, x, plane, volume]
%     boutDff = (boutData - boutBaselineRep) ./ boutBaselineRep;                                                       % --> [y, x, plane, volume]
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
% 
% %% PLOT AND SAVE VISUALIZATION OF BEHAVIOR DATA ANNOTATIONS
% 
% saveFigs = 0;
% trialGroups = []; %myData.stimSepTrials.windTrials + 2 * (~myData.stimSepTrials.windTrials);
% 
% for iFold = 0
%     
% % Create array of annotation data (row = trial, col = frame)
% annotationArr = [];
% annotTrials = 1:myData.nTrials;
% for iTrial = annotTrials(goodTrials) % Skip any trials with dropped frames 
%    annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.actionNums; %--> [trial, frame]
% end
% 
% %----------Plot 2D data with seconds on X-axis----------
% titleStr = [regexprep(expDate, '_', '\\_'), '  Behavior Summary']; % regex to add escape characters
% [~, ~, f] = plot_behavior_summary_2D(myData, annotationArr, [], titleStr, trialGroups);
% 
% %----- Plot 1D trial-averaged movement data -----
% 
% % Average across trials for wind and control trials
% annotArrLog = annotationArr ~= 0;
% annotArrSum_wind = sum(annotArrLog(myData.stimSepTrials.windTrials,:), 1);
% annotArrSum_noWind = sum(annotArrLog(~myData.stimSepTrials.windTrials,:), 1);
% 
% % Plot stimulus trial data
% h = figure(2); clf;
% h.Position = [100 50 1600 950];
% subplot(2,1,1); ax1 = gca();
% [~, ~, ~] = plot_behavior_summary_1D(myData, annotArrSum_wind, ax1, 'Wind Trials');
% 
% % Add shading during stimulus presentation
% yL = ylim();
% rectPos = [stimStart*FRAME_RATE, yL(1), (stimEnd-stimStart)*FRAME_RATE, diff(yL)]; % [x y width height]
% rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
% ylim(yL);
% 
% % Plot control trial data
% subplot(2,1,2); ax2 = gca();
% [~, ~, ~] = plot_behavior_summary_1D(myData, annotArrSum_noWind, ax2, 'Control Trials');
% 
% suptitle([regexprep(expDate, '_', '\\_'), '  Summed Movement Frames']) % regex to add escape characters
% 
% % ----------Save data as a .fig file and a .png file for each figure----------
% if saveFigs
%     
%     % Create analysis directory if necessary
%     saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', myData.sid, '\Analysis'];
%     if ~isdir(saveDir)
%         mkdir(saveDir);
%     end
%     
%     % Warn user and offer to cancel save if this will overwrite existing files
%     overwrite = 1;
%     if exist(fullfile(saveDir, [expDate, '_Behavior_Summary.fig']), 'file') ~= 0 || ...
%             exist(fullfile(saveDir, [expDate, '_Behavior_Summary.png']), 'file') ~= 0 || ...
%             exist(fullfile(saveDir, [expDate, '_Summed_Movement_Frames.fig']), 'file') ~= 0 || ...
%             exist(fullfile(saveDir, [expDate, '_Summed_Movement_Frames.png']), 'file') ~= 0
%         dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
%         if strcmp(dlgAns, 'No')
%             overwrite = 0;
%             disp('Saving cancelled')
%         end
%     end
%     
%     % Save figure files
%     if overwrite
%         savefig(f, fullfile(saveDir, [expDate, '_Behavior_Summary']));
%         export_fig(fullfile(saveDir, [expDate, '_Behavior_Summary']), '-png', f);
%         savefig(h, fullfile(saveDir, [expDate, '_Summed_Movement_Frames']));
%         export_fig(fullfile(saveDir, [expDate, '_Summed_Movement_Frames']), '-png', h);
%     end
% end%if
% end%iFold
%%% =================================================================================================
%            WIND STIMULUS ANALYSES                                   
%%%=================================================================================================
% for iFold = 1
%     %% CALCULATE MEAN dF/F FOR WIND STIM ONSET AND OFFSET
% 
%     combineStimTypes = 0;
%     baselineLength = stimDuration(2);
%     respLength = stimDuration(2);
% 
%     for iFold = 1;
%         % Divide data into different stim types
%         stimTypeData = sep_stim_types(myData, combineStimTypes); % --> [stimType]{x, y, plane, volume, trial}
% 
%         % Figure out which volumes are needed for each period
%         onsetBaselineVols = ceil((stimStart - baselineLength) * volumeRate):floor(stimStart * volumeRate);
%         offsetBaselineVols = ceil((stimEnd - baselineLength) * volumeRate):floor(stimEnd * volumeRate);
%         onsetRespVols = ceil(stimStart * volumeRate):floor(stimEnd * volumeRate);
%         offsetRespVols = ceil(stimEnd * volumeRate):floor((stimEnd + respLength) * volumeRate);
% 
%         % Initialize arrays
%         onsetBaseline = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(onsetBaselineVols), nTrials, length(stimTypeData)]);
%         offsetBaseline = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(offsetBaselineVols), nTrials, length(stimTypeData)]);
%         onsetResp = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(onsetRespVols), nTrials, length(stimTypeData)]);
%         offsetResp = zeros([size(stimTypeData{iStim}(:,:,:,1,1)), length(offsetRespVols), nTrials, length(stimTypeData)]);
%         odorOnsetDff = zeros(size(onsetBaseline(:,:,:,:,1,:)));
%         offsetDff = zeros(size(offsetBaseline(:,:,:,:,1,:)));
%         onsetDffAvg = zeros(size(onsetBaseline(:,:,:,1,1,:)));
%         offsetDffAvg = zeros(size(offsetBaseline(:,:,:,1,1,:)));
% 
%         for iStim = 1:length(stimTypeData)
% 
%             % Extract baseline and response period volumes from the each stim type
%             disp('Extracting volumes...')
%             onsetBaseline(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, onsetBaselineVols,:);                           % --> [y, x, plane, volume, trial, stimType]
%             offsetBaseline(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, offsetBaselineVols,:);                         % --> [y, x, plane, volume, trial, stimType]
%             onsetResp(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, onsetRespVols,:);                                   % --> [y, x, plane, volume, trial, stimType]
%             offsetResp(:,:,:,:,:,iStim) = stimTypeData{iStim}(:,:,:, offsetRespVols,:);                                 % --> [y, x, plane, volume, trial, stimType]
% 
%             % Calculate dF/F for trial-averaged data
%             disp('Calculating trial-averaged dF/F...')
%             odorOnsetDff(:,:,:,:,iStim) = calc_dFF(onsetResp(:,:,:,:,:,iStim), onsetBaseline(:,:,:,:,:,iStim), 5);          % --> [y, x, plane, volume, stimType]
%             offsetDff(:,:,:,:,iStim) = calc_dFF(offsetResp(:,:,:,:,:,iStim), offsetBaseline(:,:,:,:,:,iStim), 5);       % --> [y, x, plane, volume, stimType]
% 
%             % Avgerage dF/F data across volumes as well
%             disp('Calculating trial- and volume-averaged dF/F...')
%             onsetDffAvg(:,:,:,iStim) = calc_dFF(onsetResp(:,:,:,:,:,iStim), onsetBaseline(:,:,:,:,:,iStim), [4 5]);     % --> [y, x, plane, stimType]
%             offsetDffAvg(:,:,:,iStim) = calc_dFF(offsetResp(:,:,:,:,:,iStim), offsetBaseline(:,:,:,:,:,iStim), [4 5]);  % --> [y, x, plane, stimType]
%             disp('Calculations complete');
%         end
%     end%iFold
% 
% 
%                     %% CREATE VIDEO OF MEAN dF/F THROUGHOUT WIND RESPONSES
% 
%     fileName = ['Wind_Stim_Responses'];
%     sid = 1;
%     smoothingSigma = 0.5;
%     scalingFactor = 0.5;
%     offsetAlign = 0;
% 
%     % Make sure wind stim types were combined for dF/F calculation
%     if combineStimTypes
% 
%         % Select correct alignment
%         if offsetAlign
%             windDffVols = offsetDff(:,:,:,:, 1);
%             baselineVols = offsetBaselineVols;
%             respVols = offsetRespVols;
%             titleStr = 'offset';
%         else
%             windDffVols = odorOnsetDff(:,:,:,:, 1)
%             baselineVols = onsetBaselineVols;
%             respVols = onsetRespVols;
%             titleStr = 'onset';
%         end
% 
% 
%         % Calculate absolute max dF/F value across all planes and action states
%         range = calc_range(windDffVols,scalingFactor);
% 
%         % Calculate volume times in seconds relative to wind onset
%         baselineVolTimes = -(1:length(baselineVols)) / myData.volumeRate;
%         stimVolTimes = 1:length(respVols) / myData.volumeRate;
%         relTimes = [baselineVolTimes(end:-1:1), stimVolTimes];
% 
%         % Create cell array with titles for each frame
%         titleStrings = [];
%         for iVol = 1:size(windDffVols, 4)
%             if iVol <= size(baselineVols, 4)
%                 titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before wind ', titleStr];
%             else
%                 titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after wind ', titleStr];
%             end
%         end
% 
%         % Create video
%         savePath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
%         make_heatmap_vid(windDffVols, myData, range, fileName, titleStrings, savePath, [], [], [], smoothingSigma);
%     else
%         disp('Error: wind stim types must be combined to create this video')
%     end
% end%iFold
% 
%     %% CALCULATE MEAN dF/F AROUND BEHAVIOR ONSET
% for iFoldIn = 1 
%     
%     %---------- Identify behavioral state during each volume ----------
%     runVols = zeros(myData.nTrials, myData.nVolumes); 
%     stoppedVols = runVols; legMoveVols = runVols; volActions = runVols;
%     onsetVols = [];
%     for iTrial = 1:myData.nTrials
% 
%         if myData.goodTrials(iTrial)
% 
%             % Pull out action numbers for each volume
%             currActions = myData.trialAnnotations{iTrial}.actionNums;
%             volActions(iTrial,:) = currActions(volFrames);
% 
%             % Identify volume actions
%             locomotionLabel = 2; noActionLabel = 0; isoMovementLabel = 4;
%             runVols(iTrial, :) = (volActions(iTrial,:) == locomotionLabel);       % [trial, vol]
%             stoppedVols(iTrial, :) = (volActions(iTrial,:) == noActionLabel);     % [trial, vol]
%             legMoveVols(iTrial, :) = (volActions(iTrial,:) == isoMovementLabel);  % [trial, vol]
% 
%             % Find onsets of running bouts >respLen volumes in duration and preceded by >baseLen volumes of quiescence
%             baseLen = floor(baseLenSec * volumeRate);
%             respLen = floor(respLenSec * volumeRate);
%             actionPattern = [zeros(1,baseLen), (ones(1,respLen) * locomotionLabel)];
%             patternLen = length(actionPattern);
%             currTrialOnsets = strfind(volActions(iTrial, :), actionPattern);
% 
%             % Discard any bouts that occurred during or just after the wind stim
%             if myData.stimSepTrials.windTrials(iTrial)
%                 stimStartVol = ceil(stimStart * volumeRate);
%                 stimEndVol = floor(stimEnd * volumeRate);
%                 currTrialOnsets(currTrialOnsets > (stimStartVol - patternLen) & currTrialOnsets < (stimEndVol + ceil(volumeRate))) = [];
%             end
% 
%             onsetVols{iTrial} = currTrialOnsets;
% 
%         else
%             % So data from invalid trials won't ever be matched to an action state
%             runVols(iTrial, :) = 0;
%             stoppedVols(iTrial, :) = 0;
%             legMoveVols(iTrial, :) = 0;
%             volActions(iTrial, :) = -1000;
%             onsetVols{iTrial} = [];
%         end
%     end
% 
%     %---------- Get imaging data for running onsets ----------
%     onsetData = [];
%     for iTrial = 1:myData.nTrials
%         % myData.wholeSession = [y, x, plane, volume, trial]
%         if ~isempty(onsetVols{iTrial})
%             onsets = onsetVols{iTrial};
%             for iOnset = 1:length(onsets)
%                 volIdx = onsets(iOnset):onsets(iOnset) + patternLen-1;
%                 onsetData(:,:,:,:,end+1) =  myData.wholeSession(:,:,:, volIdx, iTrial); % --> [y, x, plane, onsetVolume, onsetNum]
%             end
%         end
%     end
% 
%     % Calculate dF/F before and after movment onset using pre-movement period as baseline
%     onsetBaselines = onsetData(:,:,:, 1:baseLen,:);                            % --> [y, x, plane, volume, onsetNum]
%     onsetBaselineMean = mean(mean(onsetBaselines, 5), 4);                      % --> [y, x, plane]
%     onsetBaselineMeanRep = repmat(onsetBaselineMean, 1, 1, 1, patternLen);     % --> [y, x, plane, volume]
%     onsetMean = mean(onsetData, 5);                                            % --> [y, x, plane, volume]
%     onsetDffVols = (onsetMean - onsetBaselineMeanRep) ./ onsetBaselineMeanRep; % --> [y, x, plane, volume]
%     onsetMeanDff = mean(onsetDffVols, 4);                                      % --> [y, x, plane]
% end%iFoldIn
%                     %% PLOT MOVEMENT ONSET HEATMAPS FOR EACH PLANE
% 
%     smoothingSigma = [0.5];    
%     scalingFactor = [0.75];
%     makeVid = 0;
%     saveDir = [];
%     fileName = 'Locomotion_Onset_Plane_Heatmaps';
%     plotTitle = {['dF/F - Locomotion onset - ', num2str(baseLenSec), ' sec baseline, ', num2str(respLenSec), ' sec response period']};
%     
%     % Calculate dF/F range
%     range = calc_range(onsetMeanDff, scalingFactor);
% 
%     % Plot figures
%     [~, ~] = plot_heatmaps(onsetMeanDff, myData, range, plotTitle, smoothingSigma, fileName, 'makeVid', makeVid, 'saveDir', saveDir);
