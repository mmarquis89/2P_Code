function ROIselectionGui()




% myData.wholeSession = [x, y, plane, volume, trial]

close all

%++++++++++++++++++++++++ INITIALIZATION TASKS +++++++++++++++++++++++++++++++++++++++++++++++++++++

% Prompt user for ref images data file
[dataFile, pathName, ~] = uigetfile('*.mat', 'Select reference image data file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');

if dataFile == 0
    
    % Prompt user for imaging data file to get ref images
    myData = load_imaging_data();
    refImages = myData.refImg;
    nPlanes = myData.nPlanes;
else
    % Load ref images data file
    disp(['Loading ' dataFile, '...'])
    imgData = load([pathName, dataFile]); % 1 x nPlanes cell array
    refImages = imgData.refImages;
    nPlanes = length(refImages);
    disp([dataFile, ' loaded'])
end

% Create hardcoded parameters
maxInts = [];
for iPlane = 1:nPlanes
    maxInts(iPlane) = max(refImages{iPlane}(:));
end
maxIntDefault = num2str(max(maxInts));
myData.MAX_INTENSITY = str2double(inputdlg('Enter max intensity value for reference images', '', 1, {maxIntDefault}));
MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots

myData.ROIs = [];
indexROI = 1;
roiPlaneAxes = [];
selected = 0;


%----------CONSTRUCTING COMPONENTS----------

% Figure
f = figure('position', [50 45, 1800, 950], 'Tag', 'figure', 'WindowButtonDownFcn', ...
    {@figure_WindowButtonDownFcn}, 'Name', 'Imaging Analysis GUI', 'NumberTitle', 'off');

%%%% TOP LEVEL TAB GROUP %%%%
baseTabGroup = uitabgroup(f, 'Units', 'Pixels', 'Position', [5 5 1800 940],...
    'Tag', 'baseTabGroup');

% ROI selection tab
roiTab = uitab(baseTabGroup, 'Title', 'ROI selection', 'Tag', 'roiTab');

% ROI subtab group
roiSubtabGroup = uitabgroup(roiTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99],...
    'CreateFcn', {@roiSubtabGroup_CreateFcn}, 'Tag', 'roiSubtabGroup');

% Change "Units" to normalized so components resize automatically with window
f.Units = 'normalized';
baseTabGroup.Units = 'normalized';
roiSubtabGroup.Units = 'normalized';

%---------------------------------------------------------------------------------------------------
    function roiSubtabGroup_CreateFcn(src, ~)
        
        % Create tabs for each plane and for combined ref images
        roiSelectTabs = [];
        roiRefImages = [];
        
        % Create subTabs
        for iTab = 1:nPlanes+1
            roiSelectTabs{iTab} = uitab(src);
            if iTab == 1
                
                %-----Plot combined image of all planes-----
                roiSelectTabs{1}.Title = 'All planes';
                roiSelectTabs{1}.Tag = 'All planes';
                
                % Calculate necessary axis dimensions
                subplotDim1 = ceil(sqrt(nPlanes));
                subplotDim2 = floor(sqrt(nPlanes));
                padSize = 0.01;
                rightMargin = 0.12;
                axWidth = (1 - rightMargin - ((subplotDim1 + 1) * padSize)) / subplotDim1;
                axHeight = (1 - ((subplotDim2 + 1) * padSize)) / subplotDim2;
                
                % Place axes in correct locations
                roiRefImgAxes = [];
                for yDim = 1:subplotDim2
                    for xDim = 1:subplotDim1
                        xPos = (xDim * padSize) + ((xDim-1) * axWidth);
                        yPos = (yDim * padSize) + ((yDim-1) * axHeight);
                        roiRefImgAxes{yDim, xDim} = axes(roiSelectTabs{1}, 'Units', 'Normalized', 'Position', ...
                            [xPos, yPos, axWidth, axHeight]);
                    end
                end
                
                % Order handles by z-planes, displayed from L-R, top-bottom
                roiRefImgAxes = reshape(flipud(roiRefImgAxes)', 1, nPlanes);
                
                % Plot reference frames for each plane
                myImages = [];
                for iPlane = 1:nPlanes
                    currAxes = roiRefImgAxes{iPlane};
                    axes(currAxes);
                    axis image; hold on
                    myImages{iPlane} = imshow(refImages{iPlane}, [0 MAX_INTENSITY], ...
                        'InitialMagnification', 'fit', 'Parent', currAxes);
                    myImages{iPlane}.ButtonDownFcn = {@image_ButtonDownFcn};
                    currAxes.Title.String = ['Plane #' num2str(iPlane)];
                    currAxes.Tag = ['Plane #' num2str(iPlane)];
                    axis image
                end
            else
                %-----Plot larger individual plane reference images-----
                roiSelectTabs{iTab}.Title = ['Plane #', num2str(iTab-1)];
                roiPlaneAxes{iTab-1} = axes(roiSelectTabs{iTab}, 'Position', [0.05 0.1 0.8 0.8]);
                axis image; hold on
                roiRefImages{iTab-1} = imshow(refImages{iTab-1}, [0 MAX_INTENSITY], ...
                    'InitialMagnification', 'fit', 'Parent', roiPlaneAxes{iTab-1});
                roiRefImages{iTab-1}.ButtonDownFcn = {@image_ButtonDownFcn};
                roiPlaneAxes{iTab-1}.Title.String = ['Plane #' num2str(iTab-1)];
                roiPlaneAxes{iTab-1}.Tag = ['Plane #' num2str(iTab-1)];
            end
        end
        
        % Create DrawROI button
        drawROIButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Draw ROI', 'Units', 'Normalized', ...
            'Position', [0.89, 0.6, 0.1, 0.05], 'FontSize', 12, 'Callback', {@drawROIButton_Callback}, ...
            'Tag', 'drawROIButton');
        
        % Create clearROI button
        clearROIButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Clear ROIs', 'Units', 'Normalized', ...
            'Position', [0.89, 0.5, 0.1, 0.05], 'FontSize', 12, 'Callback', {@clearROIButton_Callback}, ...
            'Tag', 'clearROIButton');
        
        % Create saveROI button
        ROISaveButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Save current plots', ...
            'Units', 'Normalized', 'Position', [0.89, 0.4, 0.1, 0.05], 'FontSize', 12, ...
            'Callback', {@ROISaveButton_Callback}, 'Tag', 'ROISaveButton');
        
    end

%---------------------------------------------------------------------------------------------------
    function figure_WindowButtonDownFcn(~, ~)
        
        % Resets all the axes titles in the ROI selection tab whenever the window is clicked
        currObj = findobj(roiSubtabGroup, 'Tag', 'All planes');
        allAxes = findobj(currObj, 'Type', 'axes');
        for iAx = 1:length(allAxes)
            allAxes(iAx).Title.String = allAxes(iAx).Tag;
            roiPlaneAxes{iAx}.Title.String = roiPlaneAxes{iAx}.Tag;
        end
        selected = 0;
    end
%---------------------------------------------------------------------------------------------------
    function image_ButtonDownFcn(src, ~)
        % Append [SELECTED] to the title of a clicked image on the "All planes" tab
        src.Parent.Title.String = [src.Parent.Title.String(end-2:end), ' [SELECTED]'];
        selected = 1;
    end
%---------------------------------------------------------------------------------------------------
    function drawROIButton_Callback(~, ~)
        
        % ---------- Draw ROI and save relevant information about it ----------        
        if selected
            cm = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Cyan'); rgb('Purple'); rgb('Brown'); ...
                rgb('Indigo'); rgb('DarkRed') ; rgb('Magenta') ; rgb('Gold')];
            currcolor = cm(mod(indexROI,size(cm,1))+1,:); % indexROI starts at 1 when gui initializes
            
            % Prompt user to create a polygon ROI
            [myData.ROIs.masks(:,:,indexROI), xi, yi] = roipoly; % --> [x, y, ROInum]
            currAxes = gca;
            
            % Save other useful information about the ROI
            numLoc = strfind(currAxes.Tag, '#');
            myData.ROIs.planes{indexROI} = currAxes.Tag(numLoc+1:end);
            myData.ROIs.plots{indexROI} = plot(xi, yi, 'linewidth', 2, 'color', currcolor);
            myData.ROIs.plotData{indexROI} = [xi, yi];
            myData.ROIs.nums{indexROI} = text(mean(xi),mean(yi),num2str(indexROI),'Color',currcolor, ...
                'FontSize',12);
            myData.ROIs.color{indexROI} = currcolor;
            
            indexROI = indexROI + 1; % Track total # of ROIs that have been drawn
        else
            disp('Must click to select a plot before drawing an ROI')
        end
    end
%---------------------------------------------------------------------------------------------------
    function clearROIButton_Callback(~, ~)
        % Clear all existing ROIs and plots
        for iROI = 1:length(myData.ROIs.plots)
            delete(myData.ROIs.plots{iROI})
            delete(myData.ROIs.nums{iROI})
        end
        myData.ROIs = [];
        indexROI = 1; % Reset global count of total # of ROIs drawn
        drawnow()
        disp('ROIs cleared')
    end
%---------------------------------------------------------------------------------------------------
    function ROISaveButton_Callback(~,~)
        
        ROIdata = myData.ROIs;
        
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



end