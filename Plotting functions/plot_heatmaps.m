function [f, plotAxes] = plot_heatmaps(dataArr, infoStruct, cLimRange, plotTitles, figPos, cMapName)
%===================================================================================================
% PLOT GCaMP dF/F HEATMAPS IN SEPARATE FIGURE WINDOWS FOR EACH PLANE
% 
% Creates a figure for each plane in the dataStruct and divides in into 2-6 subplots, with the 
% plane's mean reference image in the first subplot and a dF/F heatmap of each trial type or 
% condition distributed throughout the rest of the subplots. Returns handles to all the figures and 
% axes that were created.
% 
% INPUTS:
%       dataArr    =  a 4-D numeric array with dims [x, y, plane, plot] containing the plotting 
%                     data. The "plot" dimension will control how many subplots are created for each
%                     plane and cannot be >5.
%
%       infoStruct =  the main imaging data structure containing metadata for the experiment. 
%                     Specifically, must contain the fields "nPlanes", "refImg", and "MAX_INTENSITY".
%
%       cLimRange  =  a two-element vector [min max] to set the min and max colormap values.
%
%       plotTitles =  a cell array of strings to be used as titles for each plot. The length of the 
%                     cell array must equal size(dataArr, 4).
%
%       figPos     =  <OPTIONAL> position values [x y width height] for the figure window. Pass [] 
%                     to use the default figPos of [50 45 1800 950]
%
%       cMapName   =  <OPTIONAL> name of a colormap function to use in the heatmaps. Pass [] to use 
%                     the default colormap 'bluewhitered'.
%
% OUTPUTS:
%       f          =  Cell array with size(dataArr, 3) of the figure handles that were created.
%
%       plotAxes   =  m x(n+1) cell array containing handles to each axes that was created, where m 
%                     is the number of planes and n is the number of plot types/trial types.
%===================================================================================================

% Use default figure position and colormap if none were provided
if isempty(figPos)
    figPos = [50 45 1800 950];
end
if isempty(cMapName)
   cMapName = 'bluewhitered'; 
end

% Determine necessary subplot arrangement
nPlots = size(dataArr, 4);
if nPlots == 1;
    subplotDims = [1 2]; % Side-by-side
elseif nPlots == 2
    subplotDims = [1 3]; % Side-by-side
elseif nPlots == 3
    subplotDims = [2 2]; % Square
elseif nPlots == 4 || nPlots == 5
    subplotDims = [2 3]; % Two rows of three
else
    error('ERROR: invalid number of plots');
end

f = []; plotAxes = [];
for iPlane = infoStruct.nPlanes:-1:1 % Figure windows arranged dorsal --> ventral
    
    % Create fig
    f{iPlane} = figure(iPlane);clf
    f{iPlane}.Position = figPos; % [x y width height]
    
    % Plot reference image for the current plane
    plotAxes{iPlane, 1} = subplot(subplotDims(1), subplotDims(2), 1);
    imshow(infoStruct.refImg{iPlane}, [0 infoStruct.MAX_INTENSITY]);
    
    % Plot dF/F heatmaps for each of the other trial types
    for iPlot = 1:nPlots
        plotAxes{iPlane, iPlot+1} = subplot(subplotDims(1), subplotDims(2), iPlot+1);
        imagesc(dataArr(:,:,iPlane,iPlot));
        caxis(cLimRange);
        colormap(plotAxes{iPlane, iPlot+1}, cMapName);
        axis equal; axis off;
        title(plotTitles{iPlot});
        colorbar
    end
end

