function plot_ROI_data(ax, ROIDffAvg, infoStruct, varargin)
%===================================================================================================
% PLOT MEAN dF/F DATA WITHIN ROI THROUGHOUT TRIAL
% 
% Creates a plot of trial-averaged dF/F data for an entire trial in a specified set of axes. There 
% are various optional additions such as single trial plotting and standard deviation/stimulus event
% shading. Data is slightly smoothed by default but this can be disabled.
%
% REQUIRED INPUTS:
%       ax = the handle to the axes where the plot should be created
% 
%       ROIDffAvg = the array of dF/F data with format [volume, trial]
% 
%       infoStruct = the data structure for the current experiment. Specifically, must contain the 
%                    fields "nVolumes" and "volumeRate". 
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'OutlierSD'      = (default: 5) the number of standard deviations to use as an outlier 
%                           exclusion threshold.
% 
%       'SingleTrials'   = (default: 1) a boolean specifying whether to plot all individual trials 
%                           in the background of the mean dF/F line.
% 
%       'StdDevShading'  = (default: 1) boolean specifying whether to shade 1 standard deviation 
%                           above and below the mean dF/F line.
% 
%       'EventShading'   = (default: 1) boolean specifying whether to shade event times (e.g. odor 
%                           stims). Requires an 'AnnotationType' to work properly.
% 
%       'AnnotationType' = (default: []) an annotationType object with the onset and offset volumes 
%                           for an annotation that occurred in every trial.
%       
%       'SmoothWinSize'  = (default = 3) the width in volumes of the moving average smoothing window
%                           that will be applied to BOTH the mean dF/F plot and the individual trial 
%                           lines. Use a window size of "1" to skip smoothing entirely.
%
%===================================================================================================   

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'EventShading', 1);
addParameter(p, 'AnnotationType', []);
addParameter(p, 'SmoothWinSize', 3);
parse(p, varargin{:});
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
eventShading = p.Results.EventShading;
annotationType = p.Results.annotationType;
smoothWin = p.Results.SmoothWinSize;

% Setup variables
volTimes = [1:1:infoStruct.nVolumes] ./ infoStruct.volumeRate;

hold on
ax.YLabel.String =  'dF/F';
ax.XLabel.String = 'Time (s)';
xlim(ax, [0, max(volTimes)]);

trialAvgDff = mean(ROIDffAvg, 2);
stdDev = std(ROIDffAvg, 0, 2);

% Discard any trials that are >5 SDs from mean
outliers = zeros(1, size(ROIDffAvg, 2));
for iTrial = 1:size(ROIDffAvg, 2)
    if sum(abs(ROIDffAvg(:, iTrial) - trialAvgDff) > (outlierSD * stdDev))
        outliers(iTrial) = 1;
    end
end
ROIDffAvg(:, logical(outliers)) = [];
trialAvgDff = mean(ROIDffAvg, 2);
stdDev = std(ROIDffAvg, 0, 2);

% Plot individual trials in background
if singleTrials
    for iTrial = 1:size(ROIDffAvg, 2)
        currData = ROIDffAvg(:, iTrial);
        plot(ax, volTimes, smooth(currData, smoothWin), 'color', [0 0 1 0.25], 'LineWidth', 0.1)
    end
end

% Shade one SD above and below mean
if shadeSDs
    upper = smooth(trialAvgDff, 3) + stdDev;
    lower = smooth(trialAvgDff, 3) - stdDev;
    jbfill(volTimes, upper', lower', 'b', 'b', 1, 0.2);
end

% Plot mean response line
plot(ax, volTimes, smooth(trialAvgDff, smoothWin), 'LineWidth', 2, 'Color', 'b');

% Plot annotationTiming
if eventShading
    if ~isempty(annotationType)
        eventOnsetTimes = round(find(annotationType.onsetVols(1,:)) ./ infoStruct.volumeRate);
        eventOffsetTimes = round(find(annotationType.offsetVols(1,:)) ./ infoStruct.volumeRate);
        yL = ylim();
        for iEvent = 1:numel(eventOnsetTimes)
            onsetTime = eventOnsetTimes(iEvent);
            offsetTime = eventOffsetTimes(iEvent);
            fill(ax, [onsetTime, onsetTime, offsetTime, offsetTime], [-100, 100, 100, -100], 'r', 'facealpha', 0.20, 'edgealpha', 0); % Huge numbers so the bars don't get cut off if I increase the ylims later
        end
        ylim(yL);
    else
        disp('Warning: no annotation type was provided. Skipping annotation timing shading.')
    end
end
end