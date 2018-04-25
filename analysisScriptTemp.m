%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA

  
% Load .mat file containing trial data
myData = load_imaging_data();

% myData.goodTrials = [0, 0, 0, myData.goodTrials];
% myData.trialAnnotations = [{[], [], []}, myData.trialAnnotations];

% myData.nTrials = 121;
% myData.origFileNames(22) = [];
% myData.stimDurs(22) = [];
% myData.stimOnsetTimes(22) = [];
% myData.trialType(22) = [];
% myData.stimSepTrials.OdorA(22) = [];
% myData.stimSepTrials.OdorB(22) = [];
% myData.wholeSession(:,:,:,:, 22) = [];

for iFold = 1    
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
    stimOnsetTimes = myData.stimOnsetTimes;
    stimDurs = myData.stimDurs;
    trialDuration = myData.trialDuration;
    volumeRate = myData.volumeRate;
    volFrames = myData.volFrames;
    goodTrials = myData.goodTrials;
        
    % Create hardcoded parameters
    myData.ROIdata = [];
    myData.MAX_INTENSITY = 800; MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots
    myData.FRAME_RATE = 25; FRAME_RATE = 25; % This is the frame rate of the behavior video, not the GCaMP imaging
    if isempty(nFrames)
        nFrames = sum(trialDuration) * FRAME_RATE;
    end
    volTimes = (1:nVolumes)' ./ volumeRate;
    frameTimes = (1:nFrames)' ./ FRAME_RATE;
    
end%iFold

%% VIEW RAW DATA FOR A SINGLE TRIAL AND PLANE

planeNum = 10;
trialNum = 4; % Does not account for any skipped trials

preview_trial_movie(myData.wholeSession, planeNum, trialNum, [], [], []);
clear planeNum trialNum

%% MAKE AVERAGE FLUORESCENCE VIDEO FROM REGISTERED DATA

saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
fileName = 'postReg_average_fluorescence_by_trial';

% Setup
myVid = VideoWriter(fullfile(saveDir, fileName));
myVid.FrameRate = 1;
open(myVid);
volAvgData = squeeze(mean(wholeSession, 4)); % --> [y, x, plane, trial]
for iTrial = 1:nTrials
    
    % Create fig
    f = figure(1);clf
    f.Position = [50 45, 1620, 950];
    
    % Figure out how many subplots are needed
    nPlots = numSubplots(nPlanes);
    
    for iPlane = nPlanes:-1:1 % Reverse order so planes go from dorsal --> ventral
        
        % Plot averaged image for each plane
        ax = subaxis(nPlots(1), nPlots(2), iPlane, 'Spacing', 0, 'MB', 0.025);
        imshow(volAvgData(:,:,iPlane, iTrial), [])
        
        % Label postions
        if iPlane == nPlanes
            title('Ventral')
        elseif iPlane == 1
            title('Dorsal')
        end
        
    end%iPlane
    
    % Write frame to video
    writeFrame = getframe(f);
    writeVideo(myVid, writeFrame);
    close(f);
end
close(myVid);
clear saveDir fileName myVid volAvgData f nPlots ax writeFrame

%% INITIAL DATA PROCESSING STEPS

skipTrials = [18 19 54 55];
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
noStimAnnotArr = odorAnnotArr;
odorStims = {'OdorA', 'OdorB'};
myData.stimSepTrials.odorTrials = logical(zeros(nTrials, 1));
for iStim = 1:numel(odorStims)
    myData.stimSepTrials.odorTrials(myData.stimSepTrials.(odorStims{iStim})) = 1;
end
for iTrial = 1:nTrials
    onsetFrame = round(myData.stimOnsetTimes(iTrial)) * FRAME_RATE;
    offsetFrame = (round(myData.stimOnsetTimes(iTrial)) + round(myData.stimDurs(iTrial))) * FRAME_RATE;
    if myData.stimSepTrials.odorTrials(iTrial)
        odorAnnotArr(iTrial, onsetFrame:offsetFrame) = 4;   %--> [trial, frame]
    else
        noStimAnnotArr(iTrial, onsetFrame:offsetFrame) = 4; %--> [trial, frame]
    end
end
goodOdorTrials = logical(myData.stimSepTrials.odorTrials .* goodTrials');
odorAnnotArr(~goodOdorTrials, :) = 0;
noStimAnnotArr(~goodTrials', :) = 0;

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

% No stim "lack of events"
annotArr_NoStim = noStimAnnotArr;
annotArr_NoStim(~myData.stimSepTrials.NoStim, :) = 0;
odorAnnotations_NoStim = annotationType(myData, annotArr_NoStim, skipTrials, 'NoStim');
odorAnnotations_NoStim = get_event_vols(odorAnnotations_NoStim, '04', '40');
annotationTypes{end + 1} = odorAnnotations_NoStim;


% % Carrier stream stopping events
% annotArr_CarrierStreamStop = odorAnnotArr;
% annotArr_CarrierStreamStop(~myData.stimSepTrials.CarrierStreamStop, :) = 0;
% odorAnnotations_CarrierStreamStop = annotationType(myData, annotArr_CarrierStreamStop, skipTrials, 'carrier_stream_stop');
% odorAnnotations_CarrierStreamStop = get_event_vols(odorAnnotations_CarrierStreamStop, '04', '40');
% annotationTypes{end + 1} = odorAnnotations_CarrierStreamStop;

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
% IR laser
% ----------------------------------------------------------------------------------------------

% % Add IR laser stim frames
% laserAnnotArr = zeros(nTrials, nFrames);
% for iTrial = 1:nTrials
%     if myData.stimSepTrials.Laser(iTrial)
%         onsetFrame = round(myData.stimOnsetTimes(iTrial)) * FRAME_RATE;
%         offsetFrame = (round(myData.stimOnsetTimes(iTrial)) + round(myData.stimDurs(iTrial))) * FRAME_RATE;
%         laserAnnotArr(iTrial, onsetFrame:offsetFrame) = 3;
%     end
% end
% laserAnnotArr(~goodTrials, :) = 0;
% 
% laserAnnotations = annotationType(myData, laserAnnotArr, skipTrials, 'laser');
% laserAnnotations = get_event_vols(laserAnnotations, '03', '30');
% annotationTypes{end + 1} = laserAnnotations;

% ==================================================================================================

% Make list of all annotation type names
for iType = 1:numel(annotationTypes)
    annotationTypeNames{iType} = annotationTypes{iType}.name;
end
annotationTypeSummary = table((1:numel(annotationTypeNames))', annotationTypeNames', 'VariableNames', {'Index', 'AnnotationType'})


clear odorStims onsetFrame offsetFrame annotTrials onsetRegExpStr offsetRegExpStr

%% SEPARATE TRIALS BASED ON EVENT INTERACTIONS
analysisWindows = []; overshoots = []; filterMatches = [];

% Display available annotation event types
disp(annotationTypeSummary)

% Choose active alignment events
activeEventTypes = [5 7];
disp(annotationTypeSummary(activeEventTypes, 2))

% Choose active filter events
activeFilterTypes = [2 3];
disp(annotationTypeSummary(activeFilterTypes, 2))


% --------- ALIGNMENT EVENTS -------------

% % Odor A
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];
% 
% % Odor B
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];

% % All behavior
% analysisWindows(end+1,:) = [ 1  2 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [2];

% Locomotion
analysisWindows(end+1,:) = [ 2  3 ];
overshoots(end+1)        = 0;
filterMatches{end+1} = [];

% % Isolated Movement
% analysisWindows(end+1,:) = [ 2  2 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [3];

% Grooming
analysisWindows(end+1,:) = [ 2  3 ];
overshoots(end+1)        = 0;
filterMatches{end+1} = [];

% ------------- FILTER EVENTS -----------------

% Create filters for different condition components
allFilts = []; allFiltNames = []; filtWindows = [];

% % Odor
% withOdor =  [ 0  1  0 ];
% noOdor =    [-1 -1  0 ];
% anyOdor =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  1  0  ];
% allFilts{end+1} = [withOdor; noOdor]; %
% allFiltNames{end+1} = {'WithOdor', 'NoOdor'}; %

% Odor A
withOdorA =  [ 0  1  0 ];
noOdorA =    [-1 -1  0 ];
anyOdorA =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  1  0  ];
allFilts{end+1} = [withOdorA; noOdorA]; %
allFiltNames{end+1} = {'WithOdorA', 'NoOdorA'}; %

% Odor B
withOdorB =  [ 0  1  0 ];
noOdorB =    [-1 -1  0 ];
anyOdorB =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  1  0  ];
allFilts{end+1} = [withOdorB; noOdorB]; %
allFiltNames{end+1} = {'WithOdorB', 'NoOdorB'}; %

% % Locomotion
% startLoc = [-1  1  0 ];
% endLoc =   [ 1  0 -1 ];
% contLoc =  [ 1  1  0 ];
% noLoc =    [-1 -1 -1 ];
% anyLoc =   [ 0  0  0 ];
% withLoc =  [ 0  1  0 ];
% filtWindows(end+1,:)     = [  1 1  ];
% allFilts{end+1} = [noLoc; startLoc; contLoc]; %endMove; anyMove;startLoc; 
% allFiltNames{end+1} = { 'NoLoc', 'startLoc', 'contLoc'}; %'EndMove', 'AnyMove','StartLoc',

% % Isolated movement
% startIsoMove = [-1  1  0 ];
% endIsoMove =   [ 1  0 -1 ];
% contIsoMove =  [ 1  1  0 ];
% noIsoMove =    [-1 -1 -1 ];
% anyIsoMove =   [ 0  0  0 ];
% withIsoMove =  [ 0  1  0 ];
% afterIsoMove = [ 1  0  0 ];
% beforeIsoMove =[-1  1  0 ];
% filtWindows(end+1,:)     = [  1  1  ];
% allFilts{end+1} = [noIsoMove; beforeIsoMove; afterIsoMove]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoIsoMove', 'beforeIsoMove', 'afterIsoMove'}; %'EndMove', , 'AnyMove'

% % Grooming
% startGroom = [-1  1  0 ];
% endGroom =   [ 1  0 -1 ];
% contGroom =  [ 1  1  0 ];
% noGroom =    [-1 -1 -1 ];
% anyGroom =   [ 0  0  0 ];
% withGroom =  [ 0  1  0 ];
% filtWindows(end+1,:)     = [  1  1  ];
% allFilts{end+1} = [noGroom; withGroom]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoGroom', 'withGroom'}; %'EndMove', , 'AnyMove'
% % 
% % All behavior
% startMove = [-1  1  0 ];
% endMove =   [ 1  0 -1 ];
% contMove =  [ 1  1  0 ];
% noMove =    [-1 -1 -1 ];
% anyMove =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [ 0  0 ];
% allFilts{end+1} = [noMove; startMove; contMove]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoMove', 'StartMove', 'ContMove'}; %'EndMove', , 'AnyMove'

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
allCondFilters = []; condNames = []; activeFilterEventVols = []; onsetFilterVecs = []; offsetFilterVecs = [];
for iType = 1:nEventTypes
    
    % Eliminate filter types that overlap with the current primary event
    matchedFilterTypes = zeros(1, numel(allFilts));
    matchedFilterTypes(filterMatches{iType}) = 1;
    currFilts = allFilts(~matchedFilterTypes);
    currFiltNames = allFiltNames(~matchedFilterTypes);
    
    % Combine primary events and filter types into all filter conditions
    [allCondFilters{iType}, condNames{iType}] =  create_filter_conditions(primaryEventNames{iType}, currFilts, currFiltNames);
    
    % Select the correct eventVols for each primary event type
    activeFilterEventVols = filterEventVols(:,:,activeFilterTypes(~matchedFilterTypes));    
    
    % Get onset filter vecs for each condition
    for iCond = 1:numel(condNames{iType})
        onsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, volumeRate, 'overshoot', overshoots(iType));
    end
    
    % Get offset filter vecs for each condition
    for iCond = 1:numel(condNames{iType})
        offsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, volumeRate, 'overshoot', overshoots(iType), 'offsetAlign', 1);
    end
    
end

clear withOdor noOdor anyOdor startMove endMove contMove noMove anyMove withLaser noLaser anyLaser currType matchedFilterTypes currFilts currFiltNames

% ----------------------------------------------------------------------------------------------
% Calculate dF/F
% ----------------------------------------------------------------------------------------------

onsetDffAvg = []; offsetDffAvg = []; onsetDff = []; offsetDff = []; combinedDff = []; combinedDffAvg =[]; respRawF = []; baselineRawF = [];
onsetCondSummaries = []; offsetCondSummaries = []; allCondSummaries = [];
for iType = 1:nEventTypes
    
    primaryFiltName = primaryEventNames{iType};
    analysisWindow = analysisWindows(iType, :);
    nConds = numel(condNames{iType});
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
            
            baselineAvg = mean(mean(baselineData, 5), 4);                            % --> [y, x, plane]
            baselineRep = repmat(baselineAvg, 1, 1, 1, size(baselineData, 4));       % --> [y, x, plane, volume]
            baselineDff = calc_dFF(baselineData, baselineRep, 5);                    % --> [y, x, plane, volume]
            
            currDff = calc_dFF(respData, baselineData, 5);                           % --> [y, x, plane, volume]
            combinedVolsDff = cat(4, baselineDff, currDff);                          % --> [y, x, plane, volume]
            currDffAvg = calc_dFF(respData, baselineData, [4 5]);                    % --> [y, x, plane]
            
            baselineRawF{iType}{iCond} = baselineData;                               % --> {eventType}{Condition}[y, x, plane, volume, event]
            respRawF{iType}{iCond} = respData;                                       % --> {eventType}{Condition}[y, x, plane, volume, event]
            
            onsetDff{iType}(:,:,:,:, iCond) = combinedVolsDff;                       % --> {eventType}[y, x, plane, volume, condition]
            onsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                           % --> {eventType}[y, x, plane, condition]
            
        end
        
        % Calculate dF/F for event offsets
        if sum(offsetFilterVecs{iType}(:,iCond)) > 0
            
            [baselineData, respData] = extract_event_volumes(eventList, offsetFilterVecs{iType}(:,iCond), baselineDur, respDur, myData, ...
                'offsetAlign', 1); % --> [y, x, plane, volume, event]
            
            baselineAvg = mean(mean(baselineData, 5), 4);                            % --> [y, x, plane]
            baselineRep = repmat(baselineAvg, 1, 1, 1, size(baselineData, 4));       % --> [y, x, plane, volume]
            baselineDff = calc_dFF(baselineData, baselineRep, 5);                    % --> [y, x, plane, volume]
            
            currDff = calc_dFF(respData, baselineData, 5);                              % --> [y, x, plane, volume]
            combinedVolsDff = cat(4, baselineDff, currDff);                             % --> [y, x, plane, volume]
            currDffAvg = calc_dFF(respData, baselineData, [4 5]);                       % --> [y, x, plane]
            
            offsetDff{iType}(:,:,:,:, iCond) = combinedVolsDff;                         % --> {eventType}[y, x, plane, volume, condition]
            offsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                             % --> {eventType}[y, x, plane, condition]
        end
    end% iCond
    clear respData baselineData
    
    combinedDff{iType} = cat(5, onsetDff{iType}, offsetDff{iType});                     % --> {eventType}[y, x, plane, volume, condition]
    combinedDffAvg{iType} = cat(4, onsetDffAvg{iType}, offsetDffAvg{iType});            % --> {eventType}[y, x, plane, condition]
    clear onsetDff onsetDffAvg offsetDff offsetDffAvg
  
    % Create summary table for onset conditions
    condCountCol = (sum(onsetFilterVecs{iType})');
    alignCol = repmat({'onset'}, nConds, 1);
    baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
    respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
    rowNames = cellfun(@num2str, num2cell(1:nConds), 'uniformOutput', 0);
    varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp'};
    onsetCondSummaries{iType} = table(condCountCol, condNames{iType}, alignCol, baselineCol, respCol, 'RowNames', rowNames, 'VariableNames', varNames);
    
    % Create summary table for offset conditions
    condCountCol = (sum(offsetFilterVecs{iType})');
    alignCol = repmat({'offset'}, nConds, 1);
    baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
    respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
    rowNames = cellfun(@num2str, num2cell((1:nConds) + nConds), 'uniformOutput', 0);
    varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp'};
    offsetCondSummaries{iType} = table(condCountCol, condNames{iType}, alignCol, baselineCol, respCol, 'RowNames', rowNames, 'VariableNames', varNames);
    
    % Concatenate onset and offset summary tables
    allCondSummaries{iType} = vertcat(onsetCondSummaries{iType}, offsetCondSummaries{iType});
    allCondNames{iType} = [condNames{iType}; condNames{iType}];
    
    % Remove conditions with one or fewer occurences 
    nullConds = allCondSummaries{iType}.Count <= 1; % using one instead of zero because it causes a bug later
    combFilterVecs{iType} = [onsetFilterVecs{iType}, offsetFilterVecs{iType}]; 
    combFilterVecs{iType}(:,nullConds) = [];
    allCondSummaries{iType}(nullConds, :) = [];
    allCondSummaries{iType}.Properties.RowNames = cellfun(@num2str, num2cell(1:size(allCondSummaries{iType}, 1)), 'UniformOutput', 0);
    allCondNames{iType}(nullConds) = [];
    combinedDff{iType}(:,:,:,:, nullConds) = [];    % --> {eventType}[y, x, plane, volume, condition]
    combinedDffAvg{iType}(:,:,:, nullConds) = [];   % --> {eventType}[y, x, plane, condition]     
    disp(allCondSummaries{iType})
    
end% iType
disp('dF/F calculation complete')

clear primaryFiltName baselineDur respDur baselineAvg baselineRep baselineDff currDff combinedVolsDff currDffAvg baselineRawF respRawF 
clear condCountCol alignCol baselineCol respCol rowNames varNames offsetCondSummaries onsetCondSummaries nullConds analysisWindow eventList
end%iFold


%% =================================================================================================
%           BEHAVIOR SUMMARIES                                   
%%%=================================================================================================
%% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS

saveFig = 0;
plotTypes = [2 1]; % 1 = odor stims, 2 = behavior
s = myData.stimSepTrials;

trialGroups = [];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';
% 
% trialGroups = [[s.OdorA + 2 * s.OdorB + 3 * s.NoStim] .* goodTrials]; 
% plotTitleSuffix = ' - Ethanol\_neat (top) vs. CO2\_e-2 (mid) vs. no Stim (bottom)';%
% fileNameSuffix = '_OdorAvsOdorBvsNoStim';

for iFold = 1

% Create plot titles
nPlots = length(plotTypes);
titleStrings = [];
plotNames = [];
annotArr = [];
for iPlot = 1:nPlots
    if plotTypes(iPlot) == 1
        % Odor stim
        tempAnnotArr = annotArr_OdorA + 0.5 * annotArr_OdorB;
        plotNames{iPlot} = 'Odor Delivery';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = tempAnnotArr;
    elseif plotTypes(iPlot) == 2
         % Behavior
        plotNames{iPlot} = 'Behavior Annotation';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = behaviorAnnotArr;
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
    saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    for iPlot = 1:nPlots
        fileName = [regexprep(plotNames{iPlot}, ' ', '')];
    end
    fileName = [fileName, '_Summary ', fileNameSuffix, '_', expDate];
    
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

clear saveFig plotTypes trialGroups plotTitleSuffix fileNameSuffix s nPlots titleStrings plotNames annotArr f saveDir fileName overwrite dlgAns
end%iFold'

%% PLOT 1-D VISUALIZATION OF BEHAVIOR DATA ANNOTATIONS

%----- Plot 1D trial-averaged movement data -----
s = myData.stimSepTrials;

saveFig = 0
actionLabel = [2]; % locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
figTitle = regexprep([expDate, '  �  Fly locomotion throughout trial (red = odor)'], '_', '\\_');
plotNames = {'Ethanol\_neat', 'CO2\_e-2', 'No Stim'};

% trialGroups =  [goodTrials]; 
% % trialGroups(1:45) = 0;
% % trialGroups(45:end) = 0;
% fileNameSuffix = '_AllTrials_Locomotion';
% 
trialGroups =  [s.OdorA + 2 * s.OdorB + 3 * s.NoStim] .* goodTrials; 
% trialGroups(1:10) = 0;
trialGroups(40:end) = 0;
fileNameSuffix = '_OdorAvsOdorBvsNoStim_Locomotion_LateTrials'; 

stimShadingColors = {'red', 'green'};

for iFold = 1
   
% Get odor stim times
odorOnset = mean(myData.stimOnsetTimes(logical(s.OdorA + s.OdorB)));
odorOffset = odorOnset + mean(myData.stimDurs(logical(s.OdorA + s.OdorB)));
odorTimes = [odorOnset, odorOffset];
odorFrames = floor(odorTimes * FRAME_RATE);

% Get laser stim times
% laserOnset = mean(myData.stimOnsetTimes(logical(s.Laser)));
% laserOffset = laserOnset + mean(myData.stimDurs(logical(s.Laser)));
% laserTimes = [laserOnset, laserOffset];
% laserFrames = floor(laserTimes * FRAME_RATE);

stimShading = {odorFrames};

% Create array of annotation data
f = figure(2); clf;
f.Position = [100 100 1600 500];
f.Color = [1 1 1];

if isempty(trialGroups)
    
    % Plot summed movement data
    annotArrSum = sum(ismember(behaviorAnnotArr, actionLabel), 1) ./ nTrials;
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
        ax{iGroup} = subplot(numel(unique(trialGroups(trialGroups ~= 0))), 1, iGroup);
        annotArrSum = sum(ismember(behaviorAnnotArr(trialGroups == iGroup, :), actionLabel), 1) ./ sum(trialGroups == iGroup);
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
                rectPos = [stimStart, 0, stimLength, 1000]; % using large height value in case yLims increase later
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
    saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    fileName = regexprep(['Summed_Movement', fileNameSuffix, '_', expDate], '_', '\_');
    
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

clear s saveFig fileNameSuffix actionLabel trialGroups figTitle plotNames stimShadingColors odorOnset odorOffset odorTimes odorFrames laserOnset laserOffsetlaserTimes 
clear laserFrames stimShading f annotArrSum ax yL stimStart stimLength rectPos yLimsAll yLimMax saveDir fileName overwrite dlgAns
end%iFold

%% =================================================================================================
%            ODOR STIM ANALYSES                                   
%%%=================================================================================================
for iFoldOut = 1
    
    %% PLOT ODOR ONSET/OFFSET HEATMAPS FOR SOME TRIAL CONDITIONS
    
    % Show summary again
    odorEventName = 'odor_B';
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, odorEventName));
    currSummary = allCondSummaries{eventInd};
    currCondNames = allCondNames{eventInd};
    disp(currSummary)
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [3 6];
    makeVid = 1;
    sigma = [0.6];   
    rangeType = 'Max';
    rangeScalar = 0.8;
    saveDir = [];
    fileName = 'CO2_Response_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir);
     
clear odorEventName eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles dffCurrConds range
end%iFoldOut   
    
%% =================================================================================================
%           BEHAVIOR ANALYSES                                   
%%%=================================================================================================
for iFold = 1
    %% CALCULATE AND PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES

    locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
    actionLabel = [3];
    baselineLabel = [0];

    smoothingSigma = [0.6]; 
    rangeType = 'max';
    rangeScalar = 1;
    makeVid = 1;
    saveDir = [];
    fileName = 'Grooming_Plane_Heatmaps';
    titleStr = {'dF/F - Grooming vs. Quiescence'};
    
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
    allTrialBaselines = [];
    for iTrial = 1:myData.nTrials
        disp(['Trial ', num2str(iTrial)])
        currImgData = wholeSession(:,:,:,:,iTrial);
        currActionVols = logical(actionVols(iTrial,:));
        currStoppedVols = logical(stoppedVols(iTrial,:));
        
        % Pull out running volumes, if any exist, from the current trial
        if sum(currActionVols) > 0
            meanActionVols(:,:,:,end+1) = mean(currImgData(:,:,:,currActionVols),4);    %--> [y, x, plane, trial]
        end
        
        % Pull out stopping volumes, if any exist, from the current trial
        if sum(currStoppedVols) > 0
            meanStoppedVols(:,:,:,end+1) = mean(currImgData(:,:,:,currStoppedVols),4);  %--> [y, x, plane, trial]
        end
        
%         trialDataSorted = sort(wholeSession(:, :, :, :, iTrial), 4);                % --> [y, x, plane, volume]
%         trialBaseline = mean(trialDataSorted(:,:,:,1:round(nVolumes * 0.05)), 4);   % --> [y, x, plane]
%         allTrialBaselines(:, :, :, iTrial) = trialBaseline;                         % --> [y, x, plane, trial]
        
    end    
        
    actionMean = mean(meanActionVols, 4);          %--> [y, x, plane]
    stoppedMean = mean(meanStoppedVols, 4);        %--> [y, x, plane] 
    
    
    % Get dF/F values for action relative to quiescence
    actionDff = (actionMean - stoppedMean) ./ stoppedMean; % --> [y, x, plane]

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(actionDff, rangeScalar, rangeType);

    % Plot figures
    [f, ~] = plot_heatmaps(actionDff, myData, range, titleStr, smoothingSigma, 'fileName', fileName, 'makeVid', makeVid, ...
                           'saveDir', saveDir);

clear actionLabel baselineLabel smoothingSigma rangeType rangeScalar makeVid saveDir fileName titleStr actionVols stoppedVols currActions currActionVols
clear meanActionVols meanStoppedVols currImgData currACtionVols currStoppedVols meanActionVols meanStoppedVols actionMean stoppedMean range volActions 
clear trialDatasorted trialBaseline allTrialBaselines baselineMean
end%iFoldIn

    %% PLOT INTERACTION HEATMAPS FOR SOME TRIAL CONDITIONS
    
    % Show summary again
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, 'locomotion'));
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
                 
    currCondNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [2 3 6 5];
    sigma = [0.6]; 
    rangeType = 'Max';
    rangeScalar = 0.8;
    makeVid = 1;
    saveDir = [];
    fileName = 'Locomotion_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, myData, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 

    clear eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles range

             %% CREATE VIDEO OF MEAN dF/F FOR EACH PLANE THROUGHOUT MOVEMENT ONSET
%         
%             rangeScalar = [0.75];
%             sigma = 0.75;
%             offsetAlign = 0;
% 
% for iFold = 1
%     
%             % Select correct alignment
%             if offsetAlign
%                 combinedVolsDff = combinedVolsOffsetDff;
%                 baselineDff = offsetBaselineDff;
%                 respDff = behavOffsetDff;
%                 titleStr = 'offset';
%                  fileName = ['Behavior_Offset_Response_Heatmap_Vid_', num2str(baselineDur), '_', num2str(respDur)];
%             else
%                 combinedVolsDff = combinedVolsOnsetDff;
%                 baselineDff = onsetBaselineDff;
%                 respDff = behavOnsetDff;
%                 titleStr = 'onset';
%                 fileName = ['Behavior_Onset_Response_Heatmap_Vid', num2str(baselineDur), '_', num2str(respDur)];
%             end
% 
%             % Calculate volume times in seconds relative to odor onset
%             baselineVolTimes = -(1:size(onsetBaselineDff, 4)) / volumeRate;
%             respVolTimes = (1:size(behavOnsetDff, 4)) / volumeRate;
%             relTimes = [baselineVolTimes(end:-1:1), respVolTimes];
% 
%             % Create cell array with titles for each frame
%             titleStrings = [];
%             for iVol = 1:size(combinedVolsOnsetDff, 4)
%                 if iVol <= size(onsetBaselineDff, 4)
%                     titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec before behavior ', titleStr];
%                 else
%                     titleStrings{iVol} = ['Time = ', sprintf('%05.2f', relTimes(iVol)), ' sec after behavior ', titleStr];
%                 end
%             end
% 
%             % Create video
%             range = calc_range(combinedVolsDff, rangeScalar);
%             savePath = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
%             make_heatmap_vid(combinedVolsDff, myData, range, fileName, titleStrings, savePath, [], [], [], sigma);  
% end%iFold            

end%iFold

%% =================================================================================================
%           ROI-BASED ANALYSES                                   
%%%=================================================================================================
for iFold = 0
    
%% PLOT AND SAVE NEW ROIs

ROIselectionGui();
    
%% LOAD ROI DATA

parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(myData.sid), '\Analysis\ROIs'];
fileName = 'ROI_Data.mat';

load(fullfile(parentDir, fileName));
myData.ROIdata = ROIdata;
myData.nROIs = numel(ROIdata); nROIs = myData.nROIs;
disp('ROIs loaded')

clear parentDir filename

%% CALCULATE EVENT-TRIGGERED dF/F WITHIN ROIs

baselineROIDff = []; respROIDff = []; ROIEventDff = []; baselineDff = []; ROIEventDffAvgCell = [];
for iType = 1:nEventTypes
    
    primaryFiltName = primaryEventNames{iType};
    analysisWindow = analysisWindows(iType, :);
    baselineDur = analysisWindow(1);
    respDur = analysisWindow(2);
    
    currCondSum = allCondSummaries{iType};
    nConds = numel(currCondSum.CondName);
    
    eventList = eventLists{iType};

    for iCond = 1:nConds
        
        disp(['Extracting event data for ', primaryFiltName, ' cond #', num2str(iCond), ' of ', num2str(nConds), '...'])
        
        % Get ROI data for event onsets
        offsetAlign = strcmp(currCondSum.Align{iCond}, 'offset');
            [baselineData, respData] = extract_event_volumes(eventList, combFilterVecs{iType}(:,iCond), baselineDur, respDur, myData, ...
                'offsetAlign', offsetAlign); % --> [y, x, plane, volume, event]
                      
            currBaselineROIAvg = []; currRespROIAvg = []; currROIDff = []; baselineDff = [];
            for iROI = 1:nROIs  
                disp(['ROI #', num2str(iROI), ' of ', num2str(nROIs), '...'])
                currMask = myData.ROIdata(iROI).mask;
                currPlane = myData.ROIdata(iROI).plane;
                nEvents = size(baselineData, 5);
                
                % Baseline period
                baselineVols = size(baselineData, 4);
                currPlaneBaselineData = squeeze(baselineData(:,:,currPlane,:,:));                       % --> [y, x, volume, event]
                currPlaneBaselineData(~currMask(:,:, ones(1, baselineVols), ones(1, nEvents))) = nan;   % --> [y, x, volume, event]
                currBaselineLin = reshape(currPlaneBaselineData, size(currPlaneBaselineData, 1)*size(currPlaneBaselineData, 2), ...
                                  baselineVols, nEvents);                                               % --> [pixel, volume, event]
                currBaselineROIAvg = squeeze(mean(currBaselineLin, 1, 'omitnan'));                      % --> [volume, event]
                
                % Response period
                respVols = size(respData, 4);
                currPlaneRespData = squeeze(respData(:,:,currPlane,:,:));                               % --> [y, x, volume, event]
                currPlaneRespData(~currMask(:,:, ones(1, respVols), ones(1, nEvents))) = nan;           % --> [y, x, volume, event]
                currRespLin = reshape(currPlaneRespData, size(currPlaneRespData, 1)*size(currPlaneRespData, 2), ...
                                  respVols, nEvents);                                                   % --> [pixel, volume, event]
                currRespROIAvg = squeeze(mean(currRespLin, 1, 'omitnan'));                              % --> [volume, event]                

                % Get dF/F for entire period
                baselineMean = mean(currBaselineROIAvg, 1);                                             % --> [event]
                baselineRep = repmat(baselineMean, size(currRespROIAvg, 1), 1);                         % --> [volume, event]
                currROIDff(:,:,iROI) = (currRespROIAvg - baselineRep) ./ baselineRep;                   % --> [volume, event, ROI]
                baselineRep = repmat(baselineMean, size(currBaselineROIAvg, 1), 1);
                baselineDff(:,:,iROI) = (currBaselineROIAvg - baselineRep) ./ baselineRep;              % --> [volume, event, ROI]                                  
            end%if
            combDff = cat(1, baselineDff, currROIDff);                                                  % --> [volume, event, ROI]
            ROIEventDff{iType}{iCond} = combDff;                                                        % --> {eventType}{cond}[volume, event, ROI]
            ROIEventDffAvgCell{iType}{iCond} = squeeze(mean(combDff, 2));                               % --> {eventType}{cond}[volume, ROI]       
    
    end%iCond       
    clear respData baselineData currRespLin currRespROIAvg baselineMean baselineRep respVols currPlaneRespData currPlaneBaselineData currBaselineLin baselineDff combDff
    disp('ROI event data extracted');
end% iType
clear primaryFiltName analysisWindow eventList baselineDur respDur currMask currPlan nEvents baselineVols 

% Collapse result to a single-level cell array
ROIEventDffAvg = [];
for iType = 1:nEventTypes
    currDffAvg = ROIEventDffAvgCell{iType};         % --> {cond}[volume, ROI]}
    ROIEventDffAvg{iType} = cat(3, currDffAvg{:});  % --> {eventType}[volume, ROI, cond];
end
clear currDffAvg currDff


%% PLOT EVENT-ALIGNED dF/F WITHIN ROIs

% Show summary again
eventName = 'groom';
shadeDur = 0;
eventInd = ~cellfun(@isempty, strfind(primaryEventNames, eventName));
currSummary = allCondSummaries{eventInd};
disp(currSummary)

currConds = [2 5];
currCondNames = allCondNames{eventInd}(currConds);

currDffData = ROIEventDff{eventInd}(currConds); % --> {cond}[volume, event, ROI]

fileNamePrefix = 'EtOH_responses_';

saveDir = 0;
saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');


ROIlist = 1:size(currDffData{1}, 3);
% ROIlist = [1 2 3];

ax = []; yLims = [];
for iROI = ROIlist
    
    % Create figure
    f = figure(iROI); clf
    f.Position = [-1850 100 1500 800];
    f.Color = [1 1 1];
    
    % Plot reference image and ROI outline
    currPlane = myData.ROIdata(iROI).plane;
    xData = myData.ROIdata(iROI).xi;
    yData = myData.ROIdata(iROI).yi;
    nPlots = numel(currDffData) + 1;
    if nPlots == 3
        plotPos = [2 2]; % Because [1 3] looks bad
    else
        plotPos = numSubplots(nPlots);
    end
    subaxis(plotPos(1), plotPos(2), 1,'ML', 0.05, 'MR', 0.02, 'MT', 0.05, 'MB', 0.08, 'SV', 0.1, 'SH', 0.05)
    hold on
    imshow(myData.refImg{currPlane}, [0 MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
    title(['Plane #', num2str(myData.ROIdata(iROI).plane)])
    
    % Plot dF/F for each condition
    alignStr = currSummary.Align(currConds);
    for iPlot = 1:(nPlots - 1) 
        if strcmp(alignStr{iPlot}, 'onset')
            eventShading = [0, shadeDur];
        else
            eventShading = [-shadeDur, 0];
        end
        plotDff = currDffData{iPlot}(:,:,iROI);
        if nPlots == 3
            subaxis(plotPos(1), plotPos(2), iPlot + 2); % If there's only two plots they look better side-by-side
        elseif nPlots == 5
            subaxis( ((mod(iPlot,2) + iPlot) / 2) + iPlot ); % If there's four they should be in a square
        else
            subaxis(plotPos(1), plotPos(2), iPlot + 1)
        end
        ax{iPlot} = gca;
        volOffset = round(analysisWindows(eventInd, 1) * volumeRate * -1);
        plot_ROI_data(ax{iPlot}, plotDff, 'EventShading', eventShading, 'VolumeRate', volumeRate, 'VolOffset', volOffset, 'OutlierSD', 5); 
        title(regexprep([expDate, ' - ', currCondNames{iPlot}, ' - ', alignStr{iPlot}], '_', '\\_'))
        yLims(iPlot,:) = ylim(ax{iPlot});
    end
    
    % Scale y-axes to match
    yMin = min(yLims(:));
    yMax = max(yLims(:));
    for iPlot = 1:(nPlots - 1)
       ylim(ax{iPlot}, [yMin yMax]); 
    end
    
    % Save figure -------------------------------------------------------------------------------
    if saveDir
        fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_', expDate];
        savefig(f, fullfile(saveDir, fileName));
        export_fig(fullfile(saveDir, fileName), '-png', f);
    end
end


%% EXTRACT SESSION DATA WITHIN ROIs

nROIs = numel(myData.ROIdata);
ROIDataAvg = [];
for iROI = 1:nROIs  
    disp(['Extracting data for ROI #', num2str(iROI), ' of ', num2str(nROIs), '...'])
    currMask = myData.ROIdata(iROI).mask;
    currPlane = myData.ROIdata(iROI).plane;
    currPlaneData = squeeze(wholeSession(:,:,currPlane,:,:));                                   % --> [y, x, volume, trial]
    currPlaneData(~currMask(:,:,ones(1, nVolumes), ones(1, nTrials))) = nan;                    % --> [y, x, volume, trial]
    
    currDataLin = reshape(currPlaneData, size(currPlaneData, 1)*size(currPlaneData, 2), ...
                                nVolumes, nTrials);                                             % --> [pixel, volume, trial, ROI]
    ROIDataAvg(:,:,iROI) = squeeze(mean(currDataLin, 1, 'omitnan'));                            % --> [volume, trial, ROI] 
end
clear('currDataLin', 'currPlaneData');
disp('ROI extraction complete')
myData.ROIDataAvg = ROIDataAvg;


% CALCULATE MEAN dF/F WITHIN ROIs THROUGHOUT ENTIRE EXPERIMENT

    % Using bottom 5% of entire ROI's mean value throughout each trial as baseline
    ROIDataAvgSorted = sort(ROIDataAvg, 1);                                     % --> [volume, trial, ROI] 
    baselineMean = mean(ROIDataAvgSorted(1:round(nVolumes * 0.05), :, :), 1);   % --> [trial, ROI] 
    baselineMeanRep = baselineMean(ones(1, nVolumes), :, :);                    % --> [volume, trial, ROI] 
    ROIDffAvg = (ROIDataAvg - baselineMeanRep) ./ baselineMeanRep;              % --> [volume, trial, ROI] 
    myData.ROIDffAvg = ROIDffAvg;
    
 
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR ONE OR MORE STIM TYPES

    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, ' - EtOH_neat (top), CO2_e-2 (middle), No stim (bottom)'], '_', '\\_');
    fileNamePrefix = 'Whole_Trial_Responses_';
    s = myData.stimSepTrials;
    eventShading = [13 15];
    filterVecs = logical([s.OdorA; s.OdorB; s.NoStim] .* repmat(goodTrials, 3, 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.5;
    stdDevShading = 1;
    
    ROIlist = 1:size(ROIDffAvg, 3);
%     ROIlist = [1 2 3];
    
    yL = [];
    for iROI = ROIlist       

        nPlots = size(filterVecs, 1);
        nRows = nPlots + 1;
        
        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 1000 900];
        f.Color = [1 1 1];
        
        % Plot reference image and ROI outline
        currPlane = myData.ROIdata(iROI).plane;
        xData = myData.ROIdata(iROI).xi;
        yData = myData.ROIdata(iROI).yi;
        subaxis(nRows, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(myData.refImg{currPlane}, [0 MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
        title(['Plane #', num2str(myData.ROIdata(iROI).plane)])

        clear ax
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, filterVecs(iPlot, :), iROI); % --> [volume, trial]
            subaxis(nRows, 3, ( (iPlot*3 + 1):(iPlot*3 + 3) ))
            ax(iPlot) = gca;
            plot_ROI_data(ax(iPlot), currDffAvg, 'EventShading', eventShading, ...
                                                 'SingleTrials', singleTrials, ...
                                                 'SingleTrialAlpha', singleTrialAlpha, ...
                                                 'StdDevShading', stdDevShading, ...
                                                 'OutlierSD', 4);
            yL{iPlot} = ylim(ax(iPlot));
        end
        ax(1).Title.String = plotTitle;

        % Scale y-axes to match
        yMin = min([yL{:}]);
        yMax = max([yL{:}]);
        for iPlot = 1:nPlots
            ylim(ax(iPlot), [yMin yMax]);
        end

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_', expDate];
            savefig(f, fullfile(saveDir, fileName));
            export_fig(fullfile(saveDir, fileName), '-png', f);
        end
    end
    
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR EARLY VS LATE TRIALS

    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    
    plotTitle = regexprep([expDate, ' - No stim - early vs late trials'], '_', '\\_');
    fileNamePrefix = 'NoStim_EarlyVsLateTrials_';
    s = myData.stimSepTrials.NoStim;
    
    eventShading = [13 15];
    
    singleTrials = 1;
    singleTrialAlpha = 0.25;
    stdDevShading = 1;
    outlierSD = 2;
    legendStr = {'Trials 1:30', 'Trials 31:60', 'Trials 61:90'};
    
    trialGroups = ones(1, nTrials);
    trialGroups(1:10) = 0;
    trialGroups(30:end) = 2;
%     trialGroups(80:end) = 3;
    trialGroups = trialGroups(logical(s .* goodTrials));

%     ROIlist = 1:size(ROIDffAvg, 3);
    ROIlist = [1];

    for iROI = ROIlist       

        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 1000 900];
        f.Color = [1 1 1];

        % Plot reference image and ROI outline
        currPlane = myData.ROIdata(iROI).plane;
        xData = myData.ROIdata(iROI).xi;
        yData = myData.ROIdata(iROI).yi;
        subaxis(2, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(myData.refImg{currPlane}, [0 MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
        title(['Plane #', num2str(myData.ROIdata(iROI).plane)])
        
        % Create plot
        currDffAvg = ROIDffAvg(:,logical(s .* goodTrials), iROI); % --> [volume, trial]
        subaxis(2, 3, [4:6])
        ax = gca;
        plot_ROI_data(ax, currDffAvg, 'EventShading', eventShading, ...
                                      'TrialGroups', trialGroups,   ...
                                      'SingleTrials', singleTrials, ...
                                      'SingleTrialAlpha', singleTrialAlpha, ...
                                      'OutlierSD', outlierSD, ...
                                      'Legend', legendStr, ...
                                      'StdDevShading', stdDevShading);
        title(plotTitle);

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, '_ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_', expDate];
            savefig(f, fullfile(saveDir, fileName));
            export_fig(fullfile(saveDir, fileName), '-png', f);
        end
    end

    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL COLOR CODED BY BEHAVIOR
    
    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, ' - EtOH_neat (top), CO2_e-2 (middle), No stim (bottom)'], '_', '\\_');
    fileNamePrefix = 'Behavior_Coded_Whole_Trial_Responses_';
    s = myData.stimSepTrials;
    eventShading = [13 15];
    annotValues = [2 0];    
    
    trialFilterVecs = logical([s.OdorA; s.OdorB; s.NoStim] .* repmat(goodTrials, 3, 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.4;
    stdDevShading = 0;
    
    ROIlist = 1:size(ROIDffAvg, 3);
%     ROIlist = [1 2 3];
    
    for iROI = ROIlist       

        nPlots = size(trialFilterVecs, 1);
        nRows = nPlots + 1;
        
        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 1000 900];
        f.Color = [1 1 1];
        
        % Plot reference image and ROI outline
        currPlane = myData.ROIdata(iROI).plane;
        xData = myData.ROIdata(iROI).xi;
        yData = myData.ROIdata(iROI).yi;
        subaxis(nRows, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(myData.refImg{currPlane}, [0 MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
        title(['Plane #', num2str(myData.ROIdata(iROI).plane)])
        
        clear ax
        yL = [];
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, trialFilterVecs(iPlot, :), iROI); % --> [volume, trial]
            annotData = annotationTypes{7}.volAnnotArr;
            currAnnotArr = annotData(trialFilterVecs(iPlot, :), :);

            subaxis(nRows, 3, ( (iPlot*3 + 1):(iPlot*3 + 3) ))
            ax(iPlot) = gca; 
            
            plot_ROI_data(ax(iPlot), currDffAvg, 'AnnotArray', currAnnotArr', ... 
                                                 'AnnotValues', annotValues', ...
                                                 'EventShading', eventShading, ...
                                                 'SingleTrials', singleTrials, ...
                                                 'SingleTrialAlpha', singleTrialAlpha, ...
                                                 'StdDevShading', stdDevShading, ...
                                                 'OutlierSD', 4);            
            yL{iPlot} = ylim(ax(iPlot));
        end
        ax(1).Title.String = plotTitle;

        % Scale y-axes to match
        yMin = min([yL{:}]);
        yMax = max([yL{:}]);
        for iPlot = 1:nPlots
            ylim(ax(iPlot), [yMin yMax]);
        end

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_', expDate];
            savefig(f, fullfile(saveDir, fileName));
            export_fig(fullfile(saveDir, fileName), '-png', f);
        end
    end
    
    
    %% PLOT dF/F WITHIN ROIs THROUGHOUT EXPERIMENT
    
    % Concatate the trial-by-trial dF/F values into one long vector
    nROIs = size(ROIDffAvg, 3);
    concatROIDff = reshape(ROIDffAvg, nVolumes * nTrials, nROIs);   % --> [volume, ROI]
    
    % Get linear arrays of event annotations
    behavVols = annotationTypes{5}.eventVolsLin;
    odorAVols = annotationTypes{2}.eventVolsLin;
    odorBVols = annotationTypes{3}.eventVolsLin;
    noStimVols = annotationTypes{4}.eventVolsLin;
    groomVols = annotationTypes{7}.eventVolsLin;
    odorVols = annotationTypes{1}.eventVolsLin;
    for iROI = 1:nROIs
        
        f = figure(iROI); clf; hold on
        f.Position = [-1850 200 1800 600];
        title(num2str(iROI))
        
        % Plot dF/F and event annotations
        plot(smooth(concatROIDff(:,iROI), 3), 'b');
        plot(odorAVols, 'r');
        plot(odorBVols, 'g');
        plot(noStimVols, 'm');
        plot(behavVols, 'k');
        plot(groomVols, 'c');
        legend({'', 'EtOH', 'CO2', 'No Stim', 'Fly Movements'});
        
        % Add trial delineators and numbers
        allVols = 1:(nTrials * nVolumes);
        yL = ylim();
        for iVol = allVols(~logical(mod(allVols, nVolumes)))
            plot([iVol, iVol], yL,  'Color', 'k')
            text(iVol + 50, 0.8 * yL(2), num2str(round(iVol / nVolumes)+1));
        end
        
        xlim([0, 3 * nVolumes]);
    end
     

end%iFold 


%% =================================================================================================
%           OTHER ANALYSES                                   
%%%=================================================================================================
for iFold = 1

%% PCA

% Pull out data for one plane
planeNum = 12;

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

clear planeNum pcaData n m d data 2D tData2D coeff score latent explained coeffReshaped 

end%iFold
