function stimTypeData = sep_stim_types(infoStruct, combineStimTrials)
%================================================================================================
% 
% Separates data structure according to trial stimulus type. Can separate data into individual 
% trial types or combine all wind trials vs. control trials.
%
% INPUTS:
%       infoStruct = main experiment data structure. Must have fields "stimTypes", "wholeSession", 
%                    and "stimSepTrials". wholeSession data must be in the form: 
%                    [x, y, plane, volume, trial], and stimSepTrials must contain a "windTrials" 
%                    field.
%
%       combineStimTrials = boolean specifying whether to combine different stimulus types
%
% OUTPUTS:
%       stimTypeData = cell array of separated data in the form {[x, y, plane, volume, trial]}
%
%================================================================================================

stimTypes = infoStruct.stimTypes;

stimTypeData = [];
if ~combineStimTrials
    % Divide data by stim type
    for iStim = 1:length(stimTypes)
        stimTypeData{iStim} = infoStruct.wholeSession(:,:,:,:, infoStruct.stimSepTrials.(stimTypes{iStim})); % --> [x, y, plane, volume, trial, stimType]
    end
else
    % Divide into combined stim trials vs. control trials
    stimTypeData{1} = infoStruct.wholeSession(:,:,:,:, infoStruct.stimSepTrials.windTrials);                 % --> [x, y, plane, volume, trial, stimType]
    stimTypeData{2} = infoStruct.wholeSession(:,:,:,:, ~infoStruct.stimSepTrials.windTrials);                % --> [x, y, plane, volume, trial, stimType]
end

end