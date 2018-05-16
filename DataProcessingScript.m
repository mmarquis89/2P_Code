%==================================================================================================
%%% EXTRACT SAMPLE FRAMES FROM ANATOMY STACKS -----------------------------------------------------
%% ==================================================================================================

expDates = {...
    '2018_04_26_exp_1'
            }

targetPlanes = [ 200 ...
    ];

fileStr = '*Stack_*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    targetPlane = targetPlanes(iExp);
    
    dirPath = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    disp(['Extracting sample frames from anatomy stacks...']);
    extract_sample_frames(dirPath, fileStr, targetPlane);
    writeToLog(sprintf('%s sample frames extracted in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Extracting sample frames took ', num2str(round(toc/60, 1)) ' min']);
end
clear expDates targetPlanes fileStr dirPath
% -------------------------------------------------------------------------------------------------

%% ===================================================================================================
%%% SELECT ROI FOR OPTIC FLOW CALCULATION
%% ===================================================================================================

expDates = {...
    '2018_04_26_exp_1'
    '2018_04_26_exp_2'
            }

sids = [ ...
    0 ...
    0
    ];
     

 FRAME_RATE = 25;

 for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     sid = sids(iExp)
     
     %%% Define ROIs for optic flow combined vids
     parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
     select_video_ROIs(parentDir, sid);
     writeToLog(sprintf('%s optic flow ROIs defined', expDate));
 end
 
%% ==================================================================================================
%%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
%% ==================================================================================================

FRAME_RATE = 25;
trialDuration = 20;
expDate = '2018_04_20_exp_2';
sid = 0;

parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_Annotation.txt'];
% annotationFileName = [expDate, '_sid_', num2str(sid), '_Annotation.txt'];

tic
process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', num2str(round(toc/60, 1)) ' min']);

clear parentDir saveDir annotationFileName
% -------------------------------------------------------------------------------------------------

%% ===================================================================================================
%%% PROCESS EVENT DATA BASED ON TRIAL INTERACTIONS
%% ===================================================================================================

% Required inputs:
%       A .mat file with variables 'annotationTypes' and 'annotationTypeSummary'
%
% Outputs:
%       A .mat file with the following variables:
%           alignEventSummary
%           filterEventSummary
%           primaryEventNames
%           eventLists
%           nEventTypes
%           condNames
%           onsetFilterVecs
%           offsetFilterVecs
%===================================================================================================

volumeRate = 6.44;

saveFileName = 'EventData_Align_OdorAB_Filter_AnyMove.mat';

% Load annotation type data
[annotTypeFileName, parentDir] = uigetfile('B:\Dropbox (HMS)\2P Data\Imaging Data');
load(fullfile(parentDir, annotTypeFileName)) % Vars 'annotationTypes', 'annotationTypeSummary'

% Display available annotation event types
disp(annotationTypeSummary)

% Choose active alignment events
activeEventTypes = [2 3];
alignEventSummary = annotationTypeSummary(activeEventTypes, 2);
disp(alignEventSummary)

% Choose active filter events
activeFilterTypes = [8];
filterEventSummary = annotationTypeSummary(activeFilterTypes, 2);
disp(filterEventSummary)

% --------- ALIGNMENT EVENTS -------------
analysisWindows = []; overshoots = []; filterMatches = [];

% Odor A
analysisWindows(end+1,:) = [ 2  3 ];
overshoots(end+1)        = 1;
filterMatches{end+1} = [];

% Odor B
analysisWindows(end+1,:) = [ 2  3 ];
overshoots(end+1)        = 1;
filterMatches{end+1} = [];

% % All behavior
% analysisWindows(end+1,:) = [ 1  2 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [2];

% % Locomotion
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [];

% % Isolated Movement
% analysisWindows(end+1,:) = [ 2  2 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [3];

% % Grooming
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [];

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

% % Odor A
% withOdorA =  [ 0  1  0 ];
% noOdorA =    [-1 -1  0 ];
% anyOdorA =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  1  0  ];
% allFilts{end+1} = [withOdorA; noOdorA]; %
% allFiltNames{end+1} = {'WithOdorA', 'NoOdorA'}; %
% 
% % Odor B
% withOdorB =  [ 0  1  0 ];
% noOdorB =    [-1 -1  0 ];
% anyOdorB =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  1  0  ];
% allFilts{end+1} = [withOdorB; noOdorB]; %
% allFiltNames{end+1} = {'WithOdorB', 'NoOdorB'}; %

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

% All behavior
startMove = [-1  1  0 ];
endMove =   [ 1  0 -1 ];
contMove =  [ 1  1  0 ];
noMove =    [-1 -1 -1 ];
anyMove =   [ 0  0  0 ];
filtWindows(end+1,:)     = [ 0  0 ];
allFilts{end+1} = [noMove; startMove; contMove]; %endMove; anyMove;
allFiltNames{end+1} = {'NoMove', 'StartMove', 'ContMove'}; %'EndMove', , 'AnyMove'

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

% Save data and metadata in a .mat file
save(fullfile(parentDir, saveFileName), 'alignEventSummary', 'filterEventSummary', 'primaryEventNames', 'eventLists', ...
    'nEventTypes', 'condNames', 'onsetFilterVecs', 'offsetFilterVecs', 'analysisWindows');

clear withOdor noOdor anyOdor startMove endMove contMove noMove anyMove withLaser noLaser anyLaser currType matchedFilterTypes currFilts currFiltNames

%% ==================================================================================================
%%% ARCHIVE FILES-----------------------------------------------------------------------------------
%% ==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2018_03_16_exp_2';
parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 0;

%%% Archive raw anatomy stacks
archiveName = 'AnatomyStacks';
filterString = '*Stack*';
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% Archive raw imaging data files
archiveName = ['TrialData_sid_', num2str(sid)];
filterString = ['*sid_', num2str(sid), '_t*'];
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% Archive raw video frames
parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = ['sid_', num2str(sid), '_RawFrames'];
system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);

clear parentDir archiveName filterString
% -------------------------------------------------------------------------------------------------