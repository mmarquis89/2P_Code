function dffData = calc_dFF(respData, baselineData)
%================================================================================================
% Calculates trial-averaged dF/F from baseline and response period data. Baseline period is also
% averaged across volumes to get a steady baseline. 
%
% INPUTS:
%       respData = data from the period you want to get dF/F for, as: [x, y, plane, volume, trial]
%
%       baselineData = data in the same format, but from the baseline period only
%
% OUTPUTS:
%       dffData = trial-averaged dF/F data of form [x, y, plane, volume]
%
%================================================================================================

respMean = squeeze(mean(respData, 5));                              % --> [x, y, plane, volume]
baselineMean = squeeze(mean(mean(baselineData, 5), 4));             % --> [x, y, plane]
baselineMeanRep = repmat(baselineMean, 1, 1, 1, size(respMean, 4)); % --> [x, y, plane, volume]

dffData = (respMean - baselineMeanRep) ./ baselineMeanRep;          % --> [x, y, plane, volume]

end