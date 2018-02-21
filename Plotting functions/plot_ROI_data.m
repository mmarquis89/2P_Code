function plot_ROI_data(ax, ROIDffAvg, infoStruct, annotationType, varargin)
    



% Parse optional arguments
p = inputParser;
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'SmoothWinSize', 3);
parse(p, varargin{:});
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
smoothWin = p.Results.SmoothWinSize;

% Setup variables
volTimes = [1:1:infoStruct.nVolumes] ./ infoStruct.volumeRate;
eventOnsetTimes = round(find(annotationType.onsetVols(1,:)) ./ infoStruct.volumeRate);

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
yL = ylim();
for iOnset = 1:numel(eventOnsetTimes)
    currTime = eventOnsetTimes(iOnset);
    fill(ax, [currTime, currTime, currTime + 1, currTime + 1], [-100, 100, 100, -100], 'r', 'facealpha', 0.20, 'edgealpha', 0); % Huge numbers so the bars don't get cut off if I increase the ylims later
end
ylim(yL);

end