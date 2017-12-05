function stimTypeData = sep_stim_types(infoStruct, combineStimTrials)
%======================================================================================================
% 
% Separates data structure according to trial stimulus type. Can separate data into individual 
% trial types or combine all wind trials vs. control trials.
%
% INPUTS:
%       infoStruct = main experiment data structure. Must have fields "stimTypes", "wholeSession", 
%                    and "stimSepTrials". wholeSession data must be in the form: 
%                    [y, x, plane, volume, trial], and stimSepTrials must contain either a 
%                    "windTrials" field and an "odorTrials" field (only one should have a sum > 1).
%
%       combineStimTrials = boolean specifying whether to combine different stimulus types
%
% OUTPUTS:
%       stimTypeData = cell array of separated data in the form [stimType]{y, x, plane, volume, trial}.
%                      If combineStimTrials = TRUE, stimTypeData{1} will be the wind trials and 
%                      stimTypeData{2} will be the control trials.
%
%======================================================================================================

stimTypes = infoStruct.stimTypes;

stimTypeData = [];
if ~combineStimTrials
    
    % Divide data by stim type
    for iStim = 1:length(stimTypes)
        stimTypeData{iStim} = infoStruct.wholeSession(:,:,:,:, infoStruct.stimSepTrials.(stimTypes{iStim})); % --> {y, x, plane, volume, trial}
    end
    
else
    
    % Divide into combined stim trials vs. control trials
    if sum(infoStruct.stimSepTrials.windTrials) > 0
        stimTypeData{1} = infoStruct.wholeSession(:,:,:,:, infoStruct.stimSepTrials.windTrials);                 % --> {y, x, plane, volume, trial}
        stimTypeData{2} = infoStruct.wholeSession(:,:,:,:, ~infoStruct.stimSepTrials.windTrials);                % --> {y, x, plane, volume, trial}
    elseif sum(infoStruct.stimSepTrials.odorTrials) > 0
        stimTypeData{1} = infoStruct.wholeSession(:,:,:,:, infoStruct.stimSepTrials.odorTrials);                 % --> {y, x, plane, volume, trial}
        stimTypeData{2} = infoStruct.wholeSession(:,:,:,:, ~infoStruct.stimSepTrials.odorTrials);                % --> {y, x, plane, volume, trial}
    else
        disp('Error: the sum of either windTrials or odorTrials must be > 0')
    end
end

end