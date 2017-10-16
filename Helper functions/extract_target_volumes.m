function [trimmedData, preStimData, stimData] = extract_target_volumes(dataArr, infoStruct, baselineDurSec, respDurSec)
%==============================================================================================================
% 
% Extract imaging data from a specified period around the wind stimulus onset, and return that and/or just the
% pre-stimulus baseline period or the period during the stimulus.
%
% INPUTS:
%       dataArr        = 5-D imaging data array in the form [x, y, plane, volume, trial]
%
%       infoStruct     = main experiment data structure. Must have fields "stimStart" and "volumeRate"
%                        
%       baselineDurSec = duration of the baseline period in seconds
%
%       respDurSec     = duration of the stimulus/response period in seconds
%
% OUTPUTS:
%       trimmedData    = session data cropped to include just the specified time window
%
%       preStimData    = just the pre-stimulus baseline
% 
%       stimData       = just the post-baseline period
% 
%==============================================================================================================


% Extract variables from infoStruct
stimStart = infoStruct.stimStart;
volumeRate = infoStruct.volumeRate;

% Calculate stim period volumes
stimStartVol = ceil(stimStart*volumeRate);
baselineStartVol = stimStartVol-floor(baselineDurSec*volumeRate);
respEndVol = floor(stimStartVol + (respDurSec * volumeRate));

% Extract data
trimmedData = dataArr(:,:,:, baselineStartVol:respEndVol);
preStimData = dataArr(:,:,:, baselineStartVol:stimStartVol-1);
stimData = dataArr(:,:,:, stimStartVol:respEndVol);

end