function [annotationTypes, annotationTypeSummary] = process_annotation_types(analysisMetadata, skipTrials)

% Copy variables for convenience
nVolumes = analysisMetadata.nVolumes;
if ~isempty(analysisMetadata.nFrames)
    nFrames = analysisMetadata.nFrames;
else
    nFrames = nVolumes;
end
nTrials = analysisMetadata.nTrials;
goodTrials = analysisMetadata.goodTrials(1:nTrials);

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
analysisMetadata.stimSepTrials.odorTrials = logical(zeros(nTrials, 1));
for iStim = 1:numel(odorStims)
    analysisMetadata.stimSepTrials.odorTrials(analysisMetadata.stimSepTrials.(odorStims{iStim})) = 1;
end
for iTrial = 1:nTrials
    disp(analysisMetadata)
    onsetFrame = round(analysisMetadata.stimOnsetTimes(iTrial)) * analysisMetadata.FRAME_RATE;
    offsetFrame = (round(analysisMetadata.stimOnsetTimes(iTrial)) + round(analysisMetadata.stimDurs(iTrial))) * analysisMetadata.FRAME_RATE;
    if analysisMetadata.stimSepTrials.odorTrials(iTrial)
        odorAnnotArr(iTrial, onsetFrame:offsetFrame) = 4;   %--> [trial, frame]
    else
        noStimAnnotArr(iTrial, onsetFrame:offsetFrame) = 4; %--> [trial, frame]
    end
end
goodOdorTrials = logical(analysisMetadata.stimSepTrials.odorTrials .* goodTrials');
odorAnnotArr(~goodOdorTrials, :) = 0;
noStimAnnotArr(~goodTrials', :) = 0;

% All odor events
odorAnnotations = annotationType(analysisMetadata, odorAnnotArr, skipTrials, 'odor');
odorAnnotations = get_event_vols(odorAnnotations, '04', '40');
annotationTypes{end + 1} = odorAnnotations;

% Odor A events
annotArr_OdorA = odorAnnotArr;
annotArr_OdorA(~analysisMetadata.stimSepTrials.OdorA, :) = 0;
odorAnnotations_A = annotationType(analysisMetadata, annotArr_OdorA, skipTrials, 'odor_A');
odorAnnotations_A = get_event_vols(odorAnnotations_A, '04', '40');
annotationTypes{end + 1} = odorAnnotations_A;

% Odor B events
annotArr_OdorB = odorAnnotArr;
annotArr_OdorB(~analysisMetadata.stimSepTrials.OdorB, :) = 0;
odorAnnotations_B = annotationType(analysisMetadata, annotArr_OdorB, skipTrials, 'odor_B');
odorAnnotations_B = get_event_vols(odorAnnotations_B, '04', '40');
annotationTypes{end + 1} = odorAnnotations_B;

% No stim "lack of events"
annotArr_NoStim = noStimAnnotArr;
annotArr_NoStim(~analysisMetadata.stimSepTrials.NoStim, :) = 0;
odorAnnotations_NoStim = annotationType(analysisMetadata, annotArr_NoStim, skipTrials, 'NoStim');
odorAnnotations_NoStim = get_event_vols(odorAnnotations_NoStim, '04', '40');
annotationTypes{end + 1} = odorAnnotations_NoStim;


% % Carrier stream stopping events
% annotArr_CarrierStreamStop = odorAnnotArr;
% annotArr_CarrierStreamStop(~analysisMetadata.stimSepTrials.CarrierStreamStop, :) = 0;
% odorAnnotations_CarrierStreamStop = annotationType(analysisMetadata, annotArr_CarrierStreamStop, skipTrials, 'carrier_stream_stop');
% odorAnnotations_CarrierStreamStop = get_event_vols(odorAnnotations_CarrierStreamStop, '04', '40');
% annotationTypes{end + 1} = odorAnnotations_CarrierStreamStop;

% ----------------------------------------------------------------------------------------------
% Behavior
% ----------------------------------------------------------------------------------------------

% Add behavior annotations
behaviorAnnotArr = zeros(nTrials, nFrames);
if ~isempty(analysisMetadata.trialAnnotations)
    annotTrials = 1:nTrials;
    for iTrial = annotTrials(goodTrials) % Skip any trials with dropped frames
        behaviorAnnotArr(iTrial, :) = analysisMetadata.trialAnnotations{iTrial}.actionNums;           %--> [trial, frame]
    end
end

% Locomotion events
locomotionAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, skipTrials, 'locomotion');
locomotionAnnotations = get_event_vols(locomotionAnnotations, '[034]2', '2[034]');
annotationTypes{end + 1} = locomotionAnnotations;

% Isolated movement events
isoMoveAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, skipTrials, 'isoMove');
isoMoveAnnotations = get_event_vols(isoMoveAnnotations, '[023]4', '4[023]');
annotationTypes{end + 1} = isoMoveAnnotations;

% Grooming events
groomingAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, skipTrials, 'groom');
groomingAnnotations = get_event_vols(groomingAnnotations, '[024]3', '3[024]');
annotationTypes{end + 1} = groomingAnnotations;

% All behavior events
onsetRegExpStr = '[034]2|[023]4|[024]3';
offsetRegExpStr = '2[034]|4[023]|3[024]';
behaviorAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, skipTrials, 'move');
behaviorAnnotations = get_event_vols(behaviorAnnotations, onsetRegExpStr, offsetRegExpStr);
annotationTypes{end + 1} = behaviorAnnotations;

% ----------------------------------------------------------------------------------------------
% IR laser
% ----------------------------------------------------------------------------------------------

% % Add IR laser stim frames
% laserAnnotArr = zeros(nTrials, nFrames);
% for iTrial = 1:nTrials
%     if analysisMetadata.stimSepTrials.Laser(iTrial)
%         onsetFrame = round(analysisMetadata.stimOnsetTimes(iTrial)) * FRAME_RATE;
%         offsetFrame = (round(analysisMetadata.stimOnsetTimes(iTrial)) + round(analysisMetadata.stimDurs(iTrial))) * FRAME_RATE;
%         laserAnnotArr(iTrial, onsetFrame:offsetFrame) = 3;
%     end
% end
% laserAnnotArr(~goodTrials, :) = 0;
%
% laserAnnotations = annotationType(analysisMetadata, laserAnnotArr, skipTrials, 'laser');
% laserAnnotations = get_event_vols(laserAnnotations, '03', '30');
% annotationTypes{end + 1} = laserAnnotations;

% ==================================================================================================

% Make list of all annotation type names
for iType = 1:numel(annotationTypes)
    annotationTypeNames{iType} = annotationTypes{iType}.name;
end
annotationTypeSummary = table((1:numel(annotationTypeNames))', annotationTypeNames', 'VariableNames', {'Index', 'AnnotationType'});

end