function [f, plotAxes] = plot_heatmaps(dataArr, infoStruct, cLimRange, plotTitles, figPos, cMapName, ... 
                                       sigma, makeVid, saveDir, fileName)
%===================================================================================================
% PLOT GCaMP dF/F HEATMAPS IN SEPARATE FIGURE WINDOWS FOR EACH PLANE
% 
% Creates a figure for each plane in the dataStruct and divides it into 2-6 subplots, with the 
% plane's mean reference image in the first subplot and a dF/F heatmap of each trial type or 
% condition distributed throughout the rest of the subplots. Returns handles to all the figures and 
% axes that were created. Can optionally also save each figure as a frame in a video file.
%
% There are also buttons that allow the user to draw an ROI on any figure (after first clicking on
% the figure to select it), which will then be plotted on all the images...this is intended to 
% facilitate comparisons of a specific anatomical region across different plots. The ROI locations 
% can also be saved for use in future analyses.
% 
% INPUTS:
%       dataArr    =  a 4-D numeric array with dims [y, x, plane, plotNum] containing the plotting 
%                     data. The "plot" dimension will control how many subplots are created for each
%                     plane and cannot be >5.
%
%       infoStruct =  the main imaging data structure containing metadata for the experiment. 
%                     Specifically, must contain the fields "nPlanes", "refImg", "expDate", 
%                     and "MAX_INTENSITY" (and "sid" if no save dir is provided).
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
%       sigma      =  <OPTIONAL> the degree of gaussian smoothing to be applied to the heatmaps, 
%                     specified as a standard deviation. Pass [] to skip filtering.
%
%       makeVid    =  <OPTIONAL> boolean value indicating whether to also create and save a video 
%                     file with each frame corresponding to the contents of one plane figure.
%
%       saveDir    =  <OPTIONAL> if makeVid == 1, the directory to save the video file in. Can pass
%                     [] even if makeVid == 1 to prompt user to select a save directory.
%
%       fileName   =  <OPTIONAL> if makeVid == 1, the name for the video file to be created. Can be 
%                     [] if makeVid == 0.
%
% OUTPUTS:
%       f          =  Cell array with size(dataArr, 3) of the figure handles that were created.
%
%       plotAxes   =  m x (n+1) cell array containing handles to each axes that was created, where m
%                     is the number of planes and n is the number of plot types/trial types.
%===================================================================================================


% Prompt user for saveDir if makeVid == true and none was provided
if makeVid
    
    if isempty(saveDir)
        startDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', infoStruct.expDate, '\sid_', num2str(infoStruct.sid), '\Analysis'];
        if ~isdir(startDir)
           mkdir(startDir) 
        end
        saveDir = uigetdir(startDir, 'Select a save directory');
        if saveDir == 0
            % Display error message if user canceled without choosing a directory
            disp('ERROR: you must select a save directory or provide one as an argument');
            return
        end
    else
        % Create save dir if it doesn't already exist
        if ~isdir(saveDir)
            mkdir(saveDir)
        end
    end
    
    % Warn user and offer to cancel save if this video will overwrite an existing file
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.avi']), 'file') ~= 0
        dlgAns = questdlg('Creating this video will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    if overwrite
        % Create video writer
        myVid = VideoWriter(fullfile(saveDir, fileName));
        myVid.FrameRate = 2;
        open(myVid)
    end
    
end%if

% Use default figure position and colormap if none were provided
if isempty(figPos)
    figPos = [50 45 1800 950];
end
if isempty(cMapName)
    cMapName = 'bluewhitered';
end

% Figure
if ~isempty(findobj('Tag', 'mainHeatmapFig'));
    delete(findobj('Tag', 'mainHeatmapFig'));
end
f = figure('position', figPos, 'Name', 'Plane heatmaps', 'Tag', 'mainHeatmapFig', 'NumberTitle', ... 
            'off', 'WindowButtonDownFcn', {@figure_WindowButtonDownFcn}); clf;        

% Initialize global vars
ROIplots = [];
selected = 0;
ROIdata = [];
indexROI = 1;

% Create tab group
tabGroup = uitabgroup(f, 'Units', 'Normalized', 'Position', [0 0 1 0.99], 'Tag', 'tabGroup');

% Create DrawROI button
drawROIButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Draw ROI', 'Units', 'Normalized', ...
    'Position', [0.4, 0.965, 0.05, 0.035], 'FontSize', 12, 'Callback', {@drawROIButton_Callback}, ...
    'Tag', 'drawROIButton');

% Create clearROI button
clearROIButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Clear ROIs', 'Units', 'Normalized', ...
    'Position', [0.455, 0.965, 0.05, 0.035], 'FontSize', 12, 'Callback', {@clearROIButton_Callback}, ...
    'Tag', 'clearROIButton');

% Create saveROI button
ROISaveButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Save ROIs', 'Units', 'Normalized', ...
    'Position', [0.51, 0.965, 0.05, 0.035], 'FontSize', 12, 'Callback', {@ROISaveButton_Callback}, ...
    'Tag', 'ROISaveButton');
% --------------------------------------------------------------------------------------------------
% Note that I'm putting these buttons in this weird location because it will make them easier to 
% crop out when I save these heatmap figures as movie frames.
% --------------------------------------------------------------------------------------------------


plotAxes = []; frameStruct = []; tabs = [];
for iPlane = 1:infoStruct.nPlanes%:-1:1 % Figure windows arranged dorsal --> ventral
    
    % Create tab for current plane
    tabs{iPlane} = uitab(tabGroup);
    tabs{iPlane}.Title = ['Plane #', num2str(iPlane)];
    tabs{iPlane}.Tag = ['Plane #', num2str(iPlane)];
    
    % Determine necessary subplot arrangement
    nSubplots = size(dataArr, 4) + 1;
    if nSubplots == 3
        subplotPos = [2 2]; % numSubplots() would return [1 3] for this but that doesn't look good
    else
        subplotPos = numSubplots(nSubplots);
    end
    
    % Calculate position of each axes
    axesPos = [];
    padSize = 0.05;
    rightMargin = 0.1;
    for iAxis = 1:subplotPos(1) % rows
        for jAxis = 1:subplotPos(2) % cols
            xLen = (1 - rightMargin - (subplotPos(1) + 1) * padSize) / subplotPos(1);
            yLen = (1 - (subplotPos(2) + 1) * padSize) / subplotPos(2);
            xPos = padSize + (iAxis - 1) * (xLen + padSize);
            yPos = padSize + (jAxis - 1) * (yLen + padSize);
            axesPos(end+1,:) = [xPos, yPos, xLen, yLen];
        end
    end
    
    % Sort axes positions to start plotting in upper left corner
    axesPosSort = sortrows(axesPos, [1, -2]);
    
    % Plot reference image for the current plane
    plotAxes{iPlane, 1} = axes(tabs{iPlane}, 'Units', 'Normalized', 'Position', axesPosSort(1,:));
    hold on;
    refImg = imshow(infoStruct.refImg{iPlane}, [0 infoStruct.MAX_INTENSITY]);
    refImg.Parent.Tag = 'refImg';
    refImg.ButtonDownFcn = {@image_ButtonDownFcn};
    refImg.Parent.Title.String = regexprep(infoStruct.expDate, '_', '\\_');   
    
    % Plot dF/F heatmaps for each of the other trial types
    dffPlots = [];
    for iPlot = 1:(nSubplots - 1)
        plotAxes{iPlane, iPlot+1} = axes(tabs{iPlane}, 'Units', 'Normalized', 'Position', axesPosSort(iPlot+1,:), 'Tag', ['ax', num2str(iPlot)]);
        if ~isempty(sigma)
            dffPlots{iPlot} = imagesc(imgaussfilt(dataArr(:,:, iPlane, iPlot), sigma)); 
        else
            dffPlots{iPlot} = imagesc(dataArr(:,:, iPlane, iPlot));
        end
        dffPlots{iPlot}.ButtonDownFcn = {@image_ButtonDownFcn};
        caxis(cLimRange);
        colormap(plotAxes{iPlane, iPlot + 1}, cMapName);
        colorbar(plotAxes{iPlane, iPlot + 1});
        axis equal; axis off;
        title(plotTitles{iPlot});
        hold on;
    end
    
    % Grab video frame if makeVid == true
    if makeVid
        tabGroup.SelectedTab = tabs{iPlane}; % So the tab that was just created is shown
        drawnow();        
        frameStruct{iPlane} = getframe(f);
        frameStruct{iPlane}.cdata = frameStruct{iPlane}.cdata(40:end,:,:); % Crop out tabs at top of window
    end
    
end%for

tabGroup.SelectedTab = tabs{1};

% Create video and close videowriter if necessary
if makeVid
    for iFrame = 1:length(frameStruct)
        writeVideo(myVid, frameStruct{iFrame});
    end
    close(myVid);
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%   CALLBACK FUNCTIONS                                                                             %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

function image_ButtonDownFcn(src, ~)
    
    % Append [SELECTED] to the title of an image when it is clicked
    src.Parent.Title.String = [src.Parent.Title.String, ' [SELECTED]'];
    selected = 1; % For tracking purposes when drawing ROI
end

%---------------------------------------------------------------------------------------------------

function figure_WindowButtonDownFcn(~, ~)

    % Resets all the axes titles whenever the window is clicked
    for iPlane = 1:size(plotAxes, 1)
        currObj = findobj(tabs{iPlane}, 'Tag', 'refImg');
        currObj.Title.String = regexprep(infoStruct.expDate, '_', '\\_');   
        for iAx = 2:size(plotAxes, 2)
            plotAxes{iPlane, iAx}.Title.String = plotTitles{iAx - 1};
        end
    end
    selected = 0; % For tracking purposes when drawing ROI
end

%---------------------------------------------------------------------------------------------------
function drawROIButton_Callback(~, ~)
    
if selected
    
        % Clear any previous ROIs on this tab
        currTab = tabGroup.SelectedTab;
        oldROIs = findobj(currTab, 'Tag', 'ROIplot');
        if ~isempty(oldROIs)
            indexROI = indexROI - 1;  % Remove this from the count of total ROIs drawn
        end
        for iROI = 1:length(oldROIs)
            delete(oldROIs(iROI));
        end
        
        % Prompt user to create a polygon ROI
        [ROIdata.mask{indexROI}, ROIdata.xi{indexROI}, ROIdata.yi{indexROI}] = roipoly;
        
        % Plot the ROI on top of each image in the figure
        ROIplots = [];
        currTabAxes  = findobj(currTab, 'Type', 'Axes');
        for iROI = 1:length(currTabAxes)
            ROIplots{iROI} = plot(currTabAxes(iROI), ROIdata.xi{indexROI}, ROIdata.yi{indexROI}, 'linewidth', 2, 'color', rgb('LimeGreen'), 'tag', 'ROIplot');
        end
        
        % Record other identifying info about the ROI
        ROIdata.plane{indexROI} = tabGroup.SelectedTab.Tag;
        ROIdata.expDate{indexROI} = infoStruct.expDate;
    else
        disp('Must click to select a plot before drawing an ROI')
    end
    
    indexROI = indexROI + 1; % Track total # of ROIs that have been drawn
    
end

%---------------------------------------------------------------------------------------------------
function clearROIButton_Callback(~, ~)
    
    % Clear all existing ROIs and plots
    allROIs = findobj(tabGroup, 'tag', 'ROIplot');
    for iROI = 1:length(allROIs)
        delete(allROIs(iROI));
    end
    ROIdata = [];
    indexROI = 1; % Reset global count of total # of ROIs drawn
    drawnow()
end
%---------------------------------------------------------------------------------------------------
function ROISaveButton_Callback(~,~)

    % Prompt user for save directory
    saveDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select a save directory');
    if saveDir == 0
        % Throw error if user canceled without choosing a directory
        disp('ERROR: you must select a save directory or provide one as an argument');
        return
    end

    % Prompt user for file name
    fileName = inputdlg('Please choose a file name', 'Save ROI data', 1, {'ROI_Data'});
    fileName = fileName{:};

    % Warn user and offer to cancel save if this video will overwrite an existing file
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.mat']), 'file') ~= 0
        dlgAns = questdlg('Creating this data will overwrite an existing file in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end

    % Save ROI data
    if overwrite
        save(fullfile(saveDir, [fileName, '.mat']), 'ROIdata');
        disp('ROI data saved!')
    end

end%function
%---------------------------------------------------------------------------------------------------

end%function

