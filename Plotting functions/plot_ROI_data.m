function plot_ROI_data(ax, ROIDffAvg, infoStruct, varargin)
%===================================================================================================
% PLOT MEAN dF/F DATA WITHIN ROI THROUGHOUT TRIAL
% 
% Creates a plot of trial-averaged dF/F data for an entire trial in a specified set of axes. There 
% are various optional additions such as single trial plotting and standard deviation/stimulus event
% shading. Data is slightly smoothed by default but this can be disabled.
%
% INPUTS:
%
%       ax          = the handle to the axes where the plot should be created
% 
%       ROIDffAvg   = the array of dF/F data with format: [volume, trial]
% 
%       infoStruct  = <OPTIONAL> the data structure for the current experiment. Specifically, must 
%                     contain the fields "nVolumes", "volumeRate", "stimOnsetTimes", and "stimDurs". 
%                     Optional if you're providing your own volumes and omitting stim shading.
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'TrialGroups'    = <OPTIONAL> 1 x nTrials numeric array specifying two or more trial groups 
%                          (numbers starting at one) to split the data into. Each group will be assigned
%                          a color (the individual trials will no longer be plotted chronologically),
%                          and instead of plotting the average of all trials each group will 
%                          have its own trial-averaging line.
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
%       'EventShading'   = (default: []) two element vector specifying times in seconds to shade to
%                          indicate stimulus presentation
%       
%       'SmoothWinSize'  = (default: 3) the width in volumes of the moving average smoothing window
%                           that will be applied to BOTH the mean dF/F plot and the individual trial 
%                           lines. Use a window size of "1" to skip smoothing entirely.
% 
%       'VolumeRate'     = (default: []) rate that imaging volumes were acquired at. 
%
%       'VolOffset'      = (default: 0) a scalar number to add to the volumes for the purpose of 
%                          labeling the X-axis. For example, if you are imaging at 10 volumes/sec 
%                          and want your X-axis labeling to start at -2 seconds, use -20.
%
%       'SingleTrialAlpha' = (default: 1) a value from 0-1 specifying the transparency of the single
%                          trial lines.
%
%===================================================================================================   

% Parse optional arguments
p = inputParser;
% addParameter(p, 'InfoStruct', []);
addParameter(p, 'TrialGroups', []);
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'EventShading', []);
addParameter(p, 'SmoothWinSize', 3);
addParameter(p, 'VolumeRate', []);
addParameter(p, 'VolOffset', 0);
addParameter(p, 'SingleTrialAlpha', 1);
parse(p, varargin{:});
% infoStruct = p.Results.InfoStruct;
trialGroups = p.Results.TrialGroups;
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
eventShading = p.Results.EventShading;
smoothWin = p.Results.SmoothWinSize;
volumeRate = p.Results.VolumeRate;
volOffset = p.Results.VolOffset;
singleTrialAlpha = p.Results.SingleTrialAlpha;
if ~isempty(infoStruct)
    volumeRate = infoStruct.volumeRate;
end

% Setup variables
nVolumes = size(ROIDffAvg, 1);
volTimes = ([1:1:nVolumes] + volOffset) ./ volumeRate;
trialAvgDff = mean(ROIDffAvg, 2);
stdDev = std(ROIDffAvg, 0, 2);
nTrials = size(ROIDffAvg, 2);
if ~isempty(trialGroups)
   nGroups = numel(unique(trialGroups(trialGroups ~= 0))); 
else
   nGroups = 1;
   trialGroups = ones(nTrials, 1);
end

% Format axes
hold on
ax.YLabel.String =  'dF/F';
ax.XLabel.String = 'Time (s)';
xlim(ax, [min(volTimes), max(volTimes)]);

% Discard any trials that are >5 SDs from mean
outliers = zeros(1, size(ROIDffAvg, 2));
for iTrial = 1:size(ROIDffAvg, 2)
    if sum(abs(ROIDffAvg(:, iTrial) - trialAvgDff) > (outlierSD * stdDev))
        outliers(iTrial) = 1;
        
    end
end
if sum(outliers) > 0
    disp(['Omitting ' num2str(sum(outliers)), ' outlier trials'])
end
ROIDffAvg(:, logical(outliers)) = [];
trialGroups(logical(outliers)) = [];

if nGroups == 1
    cm = jet(nTrials);
else
    % Create custom colormap if using trial groups
    cm = [rgb('Blue'); ...
          rgb('Green'); ...
          rgb('Red'); ...
          rgb('Magenta'); ...
          rgb('Cyan'); ...
          rgb('Gold'); ...
          rgb('DarkRed'); ...
          rgb('Yellow'); ...
          rgb('Lime') ...
          ];
end

for iGroup = 1:nGroups
    
    % Select plotting colors
    if nGroups == 1
        shadeColor = [0 0 1];
        meanLineColor = [0 0 0];
    else
        shadeColor = cm(iGroup, :);
        meanLineColor = shadeColor;
    end
    
    % Separate data from current trial group and calculate average dF/F
    groupAvgDff = mean(ROIDffAvg(:, trialGroups == iGroup), 2);
    groupDff = ROIDffAvg(:, trialGroups == iGroup);
    groupStdDev = std(groupDff, 0, 2);
    
    % Plot individual trials in background
    if singleTrials
        for iTrial = 1:size(groupDff, 2)
            currData = groupDff(:, iTrial);
            if nGroups == 1
                currColor = cm(iTrial, :);
            else
                currColor = cm(iGroup, :);
            end
            plt = plot(ax, volTimes, smooth(currData, smoothWin), 'color', currColor, 'LineWidth', 0.1);
            plt.Color(4) = singleTrialAlpha; 
        end
    end
    
    % Shade one SD above and below mean
    if shadeSDs
        upper = smooth(groupAvgDff, 3) + groupStdDev;
        lower = smooth(groupAvgDff, 3) - groupStdDev;
        jbfill(volTimes, upper', lower', shadeColor, shadeColor, 1, 0.2);
    end
    
    % Plot mean response line
    plot(ax, volTimes, smooth(groupAvgDff, smoothWin), 'LineWidth', 2, 'Color', meanLineColor * 1);
    
end

% Plot alignment line if applicable
yL = ylim();
if volTimes(1) < 0
   plot(ax, [0 0], [-100 100], 'Color', 'k', 'LineWidth', 2); % Huge Y-range so it doesn't get cut off if I increase the ylims later
end
ylim(yL);

% Plot event timing
if ~isempty(eventShading)
    onsetTime = eventShading(1);
    offsetTime = eventShading(2);
    fill(ax, [onsetTime, onsetTime, offsetTime, offsetTime], [-100, 100, 100, -100], 'r', 'facealpha', 0.1, 'edgealpha', 0); % Huge numbers so the bars don't get cut off if I increase the ylims later
    ylim(yL);
end
end