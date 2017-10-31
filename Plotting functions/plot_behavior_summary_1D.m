function [plotHandle, plotAxes, plotFig] = plot_behavior_summary_1D(infoStruct, annotDataArr, plotAxes, titleStr)
%=======================================================================================================
% PLOT A 1D SUMMARY FIGURE OF BEHAVIOR ANNOTATION DATA THROUGHOUT EXPERIMENT
% 
% Plots the total number of trials when the fly was moving for each frame of the behavior video. Can 
% create a new figure for the plot or plot in a specific axes object. Plotted data is gently smoothed
% (each point averaged with the two adjacent points).
%
% INPUTS:
%       infoStruct    = the main imaging data structure containing metadata for the experiment. 
%                         Specifically, must contain the fields "nFrames", and "trialDuration".                         
%       annotationArr = 1 x nFrames array of behavioral annotation data to be plotted
%       plotAxes      = <OPTIONAL> the handle to the axes you want the figure to be plotted in. Pass 
%                          [] to create a new figure and axes for the plot.
%       titleStr      = <OPTIONAL> string to be used as the title of the plot.
%       
% OUTPUTS:
%       plotHandle  = the handle to the plot that was created by the function
%       plotAxes    = the axes that the figure was plotted in
%       plotfig     = handle to the new figure, if one was created (otherwise returns [])
%========================================================================================================

% Create or select figure and axes depending on whether an axes handle was provided
if isempty(plotAxes)
    plotFig = figure(1); clf; 
    plotAxes = axes();
    plotFig.Position = [200 100 1120 840];
    plotFig.Color = [1 1 1];
else
    axes(plotAxes)
    plotFig = [];
end

% Plot data
plotHandle = plot(smooth(annotDataArr,3));

% Format axes
plotHandle.LineWidth = 1;
plotAxes.FontSize = 12;
plotAxes.XLabel.String = 'Time (sec)';
plotAxes.XLabel.FontSize = 16;
plotAxes.YLabel.String = 'Total movement frames';
plotAxes.YLabel.FontSize = 14;
plotAxes.Title.String = titleStr;
plotAxes.XTick = [0:(1/sum(infoStruct.trialDuration)):1]*infoStruct.nFrames;
plotAxes.XTickLabel = [0:(1/sum(infoStruct.trialDuration)):1] * sum(infoStruct.trialDuration);
plotAxes.XLim = [0 infoStruct.nFrames];

end