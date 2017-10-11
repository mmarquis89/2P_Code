function ImagingAnalysisGui()
%
%
%
%
%
%
%
% myData.wholeSession = [x, y, plane, volume, trial]

close all

%++++++++++++++++++++++++ INITIALIZATION TASKS +++++++++++++++++++++++++++++++++++++++++++++++++++++

% Prompt user for data file and load data
myData = load_imaging_data();
if ~isempty(myData) % Abort initialization if no data was loaded
    
    expDate = myData.expDate;
    nPlanes = myData.nPlanes;
    nVolumes = myData.nVolumes;
    refImg = myData.refImg;
    nFrames = myData.nFrames;
    stimTypes = myData.stimTypes;
    
    % Create global/hardcoded variables
    myData.maxIntensity = 1800; maxIntensity = myData.maxIntensity; % To control brightness of ref images
    myData.volumeRate = 6.5; volumeRate = myData.volumeRate;
    myData.trialDuration = 16; trialDuration = myData.trialDuration;
    myData.stimDuration = [4 trialDuration-7]; stimDuration = myData.stimDuration; % [stim start time, stim length] in seconds
    myData.stimStart = myData.stimDuration(1); stimStart = myData.stimStart;
    myData.stimEnd = sum(myData.stimDuration); stimEnd = myData.stimEnd;
    myData.frameRate = 25; frameRate = 25;
    
    myData.ROIs = [];
    myData.dffData = [];
    plotAxes2D = [];
    plotAxes1D_wind = [];
    plotAxes1D_noWind = [];
    index = 1;
    roiPlaneAxes = [];
    dffPlotCounter = 1;
    
    
    
    %----------CONSTRUCTING COMPONENTS----------
    
    % Figure
    f = figure('position', [50 45, 1800, 950], 'Tag', 'figure', 'WindowButtonDownFcn', ...
        {@figure_WindowButtonDownFcn}, 'Name', 'Imaging Analysis GUI');
    
    %%%% TOP LEVEL TAB GROUP %%%%
    baseTabGroup = uitabgroup(f, 'Units', 'Pixels', 'Position', [5 5 1800 940],...
        'CreateFcn', {@baseTabGroup_CreateFcn}, 'Tag', 'baseTabGroup');
    
    % ROI selection tab
    roiTab = uitab(baseTabGroup, 'Title', 'ROI selection', 'Tag', 'roiTab');
    
        % ROI subtab group
        roiSubtabGroup = uitabgroup(roiTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99],...
            'CreateFcn', {@roiSubtabGroup_CreateFcn}, 'Tag', 'roiSubtabGroup');
    
    % Behavior data summary tab
    behavSummaryTab = uitab(baseTabGroup, 'Title', 'Behavior summary', 'Tag', 'behavSummaryTab');
    
        % Behavior summary subtab group
        behavSummarySubtabGroup = uitabgroup(behavSummaryTab, 'Units', 'Normalized', ...
            'Position', [0 0 1 0.99], 'CreateFcn', {@behavSummarySubtabGroup_CreateFcn}, ...
            'Tag', 'behavSummarySubtabGroup');
    
    % Analyze stim responses tab
    stimRespTab = uitab(baseTabGroup, 'Title', 'Analyze stim responses', 'Tag', 'stimRespTab');
    
        % Wind response subtab group
        stimRespSubtabGroup = uitabgroup(stimRespTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99], ...
            'CreateFcn', {@stimRespSubtabGroup_CreateFcn}, 'Tag', 'stimRespSubtabGroup');
    
    % Analyze behavior responses Tab
    behavRespTab = uitab(baseTabGroup, 'Title', 'Analyze behavior responses', 'Tag', 'behavResptab');
    
        % Behavior response subtab group
        behavRespSubtabGroup = uitabgroup(behavRespTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99], ...
            'CreateFcn', {@behavRespSubtabGroup_CreateFcn}, 'Tag', 'behavRespSubtabGroup');
    
    
    
    
    
    
    
    %
    %     % Create "Plot dF/F" button
    %     plotDffButton = uicontrol('Style', 'pushbutton', 'String', 'Plot dF/F',...
    %         'Position', [f.Position(3)-200, f.Position(4)-300, 180, 50], 'FontSize', 12, 'Callback', ...
    %         {@plotDffButton_Callback}, 'Tag', 'plotDffButton', 'Units', 'Normalized');
    %
    %     % Create "Save dF/F" button
    %     saveDffButton = uicontrol('Style', 'pushbutton', 'String', 'Save dF/F data',...
    %         'Position', [f.Position(3)-200, f.Position(4)-400, 180, 50], 'FontSize', 12, 'Callback', ...
    %         {@saveDffButton_Callback}, 'Tag', 'saveDffButton', 'Units', 'Normalized');
    
    % Should maybe have a "load ROIs" button as well?
    
    % Change "Units" to normalized so components resize automatically with window
    f.Units = 'normalized';
    baseTabGroup.Units = 'normalized';
    roiSubtabGroup.Units = 'normalized';
    behavSummarySubtabGroup.Units = 'normalized';
    stimRespSubtabGroup.Units = 'normalized';
    behavRespSubtabGroup.Units = 'normalized';
    
end%if


%+++++++++++++ CALLBACK FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
                    myImages{iPlane} = imshow(refImg{iPlane}, [0 maxIntensity], ...
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
                roiRefImages{iTab-1} = imshow(refImg{iTab-1}, [0 maxIntensity], ...
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
        
    end
%---------------------------------------------------------------------------------------------------
    function behavSummarySubtabGroup_CreateFcn(src, ~)
        
        if ~isempty(myData.trialAnnotations)
            
            %----- Create and populate 2D summary tab -----
            summary2D = uitab(src, 'Title', '2D summary', 'Tag', 'summary2D');
            plotAxes2D = axes(summary2D, 'Tag', 'plotAxes2D', 'OuterPosition', [0.01, 0.01, 0.85, 1]);
            
            % Create array of annotation data (row = trial, col = frame)
            annotationArr = [];
            annotTrials = 1:myData.nTrials;
            for iTrial = annotTrials(myData.goodTrials) % Skip any trials with dropped frames
                annotationArr(iTrial,:) = myData.trialAnnotations{iTrial}.actionNums; %--> [trial, frame]
            end
            
            [plot2D,~,~] = plot_behavior_summary_2D(myData, annotationArr, plotAxes2D, [], []);
            
            % ----- Create and populate 1D summary tab -----
            summary1D = uitab(src, 'Title', 'Trial-averaged summary', 'Tag', 'summary1D');
            plotAxes1D_wind = axes(summary1D, 'Tag', 'plotAxes1D_wind', ...
                'OuterPosition',[0, 0.525, 0.9, 0.45]);
            plotAxes1D_noWind = axes(summary1D, 'Tag', 'plotAxes1D_noWind', ...
                'OuterPosition', [0, 0.015, 0.9, 0.45]);
            
            % Get trial-averaged movment data for wind and control trials
            annotArrLog = annotationArr ~= 0;
            annotArrSum_wind = sum(annotArrLog(myData.stimSepTrials.windTrials,:), 1);
            annotArrSum_noWind = sum(annotArrLog(~myData.stimSepTrials.windTrials,:), 1);
            
            % Make plots
            plot_behavior_summary_1D(myData, annotArrSum_wind, ...
                plotAxes1D_wind, 'Wind Trials');
            plot_behavior_summary_1D(myData, annotArrSum_wind, ...
                plotAxes1D_noWind, 'Control Trials');
            
            % Add shading during stimulus presentation
            yL = ylim(plotAxes1D_wind);
            rectPos = [stimStart*frameRate, yL(1), (stimEnd-stimStart)*frameRate, diff(yL)]; % [x y width height]
            rectangle(plotAxes1D_wind, 'Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');
            ylim(yL);
            
            % Create customization options?
            
            % Create save button
            behavSummarySaveButton = uicontrol(behavSummaryTab, 'Style', 'pushbutton', ...
                'String', 'Save plot(s)', 'Units', 'Normalized', 'Position', [0.89, 0.5, 0.1, 0.05], ...
                'FontSize', 12, 'Callback', {@behavSummarySaveButton_Callback}, ...
                'Tag', 'behavSummarySaveButton');
        end
    end
%---------------------------------------------------------------------------------------------------
    function stimRespSubtabGroup_CreateFcn(src, ~)
        
        %----- CREATE PLOTTING PARAMETERS TAB -----
        stimRespParams = uitab(src, 'Title', 'Plotting Parameters', 'Tag', 'stimRespParams');
        
        
        %%%% dF/F Heatmaps %%%%
        
        % Create stim type selection button group
        stimTypeButtonGroup = uibuttongroup(stimRespParams, 'Position', [0.01 0.8 0.18 0.15], ...
            'Tag', 'stimTypeButtonGroup');
        trialTypeRadio = uicontrol(stimTypeButtonGroup, 'Style', 'radiobutton', 'String', ...
            'Plot each trial type separately', 'Units', 'Normalized', 'Tag', 'trialTypeRadio', ...
            'Position', [0.05 0.6 0.95 0.2], 'FontSize', 12);
        stimControlRadio = uicontrol(stimTypeButtonGroup, 'Style', 'radiobutton', 'String', ...
            'Plot stim trials vs. control trials', 'Units', 'Normalized', 'Tag', 'stimControlRadio', ...
            'Position', [0.05, 0.2, 0.95, 0.2], 'FontSize', 12);
        
        % Create baseline and response duration text boxes and labels
        baselineDurBox = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'edit', ...
            'Position', [0.33 0.89 0.025 0.05], 'FontSize', 12, 'Tag', 'baselineDurBox');
        baselineDurBoxLabel = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'Text', ...
            'Position', [0.19 0.88 0.135 0.05], 'String', 'Pre-stim baseline duration (sec)', ...
            'HorizontalAlignment', 'Right', 'FontSize', 12, 'Tag', 'baselineDurBoxLabel');
        respDurBox = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'edit', ...
            'Position', [0.33 0.82 0.025 0.05], 'FontSize', 12, 'Tag', 'respDurBox');
        respDurBoxLabel = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'Text', ...
            'Position', [0.19 0.81 0.135 0.05], 'String', 'Response duration (sec)', ...
            'HorizontalAlignment', 'Right', 'FontSize', 12, 'Tag', 'respDurBoxLabel');
        
        
        % Create label with overall trial duration parameters
        
        
        % Create "Plot heatmaps" button
        plotStimHeatmapsButton = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'pushbutton', ...
            'Position', [0.38 0.855 0.135 0.05], 'FontSize', 12, 'String', 'Plot dF/F Heatmaps', ...
            'Callback', {@plotStimHeatmapsButton_Callback}, 'Tag', 'plotStimHeatmapsButton');
        
        % Create "Plot ROI stim responses" button
        plotROIstimResponsesButton = uicontrol(stimRespParams, 'Units', 'Normalized', 'Style', 'pushbutton', ...
            'Position', [0.54 0.855 0.135 0.05], 'FontSize', 12, 'String', 'Plot stim responses in ROIs', ...
            'Callback', {@plotROIstimResponsesButton_Callback}, 'Tag', 'plotROIstimResponsesButton');
        
        
    end
%---------------------------------------------------------------------------------------------------
    function figure_WindowButtonDownFcn(~, ~)
        
        % Resets all the axes titles in the ROI selection tab whenever the window is clicked
        currObj = findobj('Tag', 'All planes');
        allAxes = currObj.Children;
        for iAx = 1:length(allAxes)
            allAxes(iAx).Title.String = allAxes(iAx).Tag;
            roiPlaneAxes{iAx}.Title.String = roiPlaneAxes{iAx}.Tag;
        end
    end
%---------------------------------------------------------------------------------------------------
    function image_ButtonDownFcn(src, ~)
        % Append [SELECTED] to the title of a clicked image on the "All planes" tab
        src.Parent.Title.String = [src.Parent.Title.String(end-2:end), ' [SELECTED]'];
    end
%---------------------------------------------------------------------------------------------------
    function plotStimHeatmapsButton_Callback(~,~)
        
        % Find parameter objects
        stimTypeButtonGroup = findobj('Tag', 'stimTypeButtonGroup');
        baselineDurBox = findobj('Tag', 'baselineDurBox');
        respDurBox = findobj('Tag', 'respDurBox');
        
        % Delete previous dF/F heatmap tab if it exists
        if ~isempty(findobj('Tag', 'stimHeatmapTab'))
            delete(findobj('Tag', 'stimHeatmapTab'));
            drawnow();
        end
        
        % Create new tab for heatmap plots
        stimHeatmapsTab = uitab(stimRespSubtabGroup, 'Title', 'dF/F Heatmaps', 'Tag', 'stimHeatmapTab');
        
        % Create dFF heatmap subtab group
        stimHeatmapSubtabGroup = uitabgroup(stimHeatmapsTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99],...
            'Tag', 'stimHeatmapSubtabGroup');
        
        % Create save button
        stimHeatmapSaveButton = uicontrol(stimHeatmapsTab, 'Style', 'pushbutton', 'String', 'Save current plots', ... 
            'Units', 'Normalized', 'Position', [0.89, 0.6, 0.1, 0.05], 'FontSize', 12, ...
            'Callback', {@heatmapSaveButton_Callback}, 'Tag', 'stimHeatmapSaveButton'); 
        
        %----- Calculate mean dF/F for baseline and response periods -----
        
        % Determine which volumes are in the baseline and response periods
        baselineDurSec = str2num(baselineDurBox.String);
        baselineStart = floor(stimStart*volumeRate)-floor(baselineDurSec*volumeRate);
        baselineEnd = floor(stimStart*volumeRate);
        respDur = str2num(respDurBox.String);
        respStart = ceil(stimStart*volumeRate);
        respEnd = floor((stimStart+respDur)*volumeRate);
        
        combineStimTrials = ~strcmp(stimTypeButtonGroup.SelectedObject.Tag, 'trialTypeRadio');
        if ~combineStimTrials
            % Calculate dF/F for each trial type separately
            for iStim = 1:length(stimTypes)
                % Get trial averaged data
                wholeTrialAvg(:,:,:,:,iStim) = mean(myData.wholeSession(:,:,:,:,myData.stimSepTrials.(stimTypes{iStim})), 5); % --> [x, y, plane, volume, StimType]
                baselineAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:, baselineStart:baselineEnd, iStim), 4);                   % --> [x, y, plane, StimType]
                stimAvg(:,:,:,iStim) = mean(wholeTrialAvg(:,:,:, respStart:respEnd, iStim), 4);                               % --> [x, y, plane, StimType]
            end%for
        else
            % Get trial averaged data for all combined stim trials vs. control trials
            wholeTrialAvg(:,:,:,:,1) = squeeze(mean(myData.wholeSession(:,:,:,:,myData.stimSepTrials.windTrials), 5));   % --> [x, y, plane, volume, StimType]
            baselineAvg(:,:,:,1) = squeeze(mean(wholeTrialAvg(:,:,:, baselineStart:baselineEnd, 1), 4));                 % --> [x, y, plane, StimType]
            stimAvg(:,:,:,1) = squeeze(mean(wholeTrialAvg(:,:,:, respStart:respEnd,1), 4));                              % --> [x, y, plane, StimType]
            wholeTrialAvg(:,:,:,:,2) = squeeze(mean(myData.wholeSession(:,:,:,:,~myData.stimSepTrials.windTrials), 5));  % --> [x, y, plane, volume, StimType]
            baselineAvg(:,:,:,2) = squeeze(mean(wholeTrialAvg(:,:,:, baselineStart:baselineEnd, 2), 4));                 % --> [x, y, plane, StimType]
            stimAvg(:,:,:,2) = squeeze(mean(wholeTrialAvg(:,:,:, respStart:respEnd,2), 4));                              % --> [x, y, plane, StimType]           
        end%if
        
        % Calculate dF/F
        dffAvg = (stimAvg - baselineAvg) ./ baselineAvg; % --> [x, y, plane, StimType]

        % Calculate absolute max dF/F value across all planes and stim types
        range = calc_range(dffAvg, []);
        
        %----- Create subTabs and plot heatmap data -----
        if combineStimTrials
            nTabs = nPlanes + 1;
        else
            nTabs = nPlanes;
        end
        for iTab = 1:nTabs
            stimHeatmapSubtabs{iTab} = uitab(stimHeatmapSubtabGroup);
            if iTab == 1 && combineStimTrials
                
                %----- Plot all planes in a single tab for combined stim trials only -----
                stimHeatmapSubtabs{1}.Title = 'All planes';
                stimHeatmapSubtabs{1}.Tag = 'All planes';
                
                % Calculate necessary axis dimensions
                subplotDim1 = ceil(sqrt(nPlanes));
                subplotDim2 = floor(sqrt(nPlanes));
                padSize = 0.01;
                rightMargin = 0.12;
                axWidth = (1 - rightMargin - ((subplotDim1 + 1) * padSize)) / subplotDim1;
                axHeight = (1 - ((subplotDim2 + 1) * padSize)) / subplotDim2;

                % Place axes in correct locations
                allPlanesAxes = [];
                for yDim = 1:subplotDim2
                    for xDim = 1:subplotDim1
                        xPos = (xDim * padSize) + ((xDim-1) * axWidth);
                        yPos = (yDim * padSize) + ((yDim-1) * axHeight);
                        allPlanesAxes{yDim, xDim} = axes(stimHeatmapSubtabs{1}, 'Units', 'Normalized', ...
                            'Position', [xPos, yPos, axWidth, axHeight]);
                    end
                end

                % Order handles by z-planes, displayed from L-R, top-bottom
                allPlanesAxes = reshape(flipud(allPlanesAxes)', 1, nPlanes);
                
                % Plot heatmaps for each plane (stim trials only)
                allPlanesHeatmapPlots = [];
                for iPlane = 1:nPlanes
                    currAxes = allPlanesAxes{iPlane};
                    axes(currAxes);
                    axis image; hold on
                    allPlanesHeatmapPlots{iPlane} = imagesc(dffAvg(:,:,iPlane, 1));
                    caxis(range);
                    colormap(currAxes, 'bluewhitered');
                    axis image; axis off;
                    currAxes.Title.String = ['Plane #' num2str(iPlane)];
                    currAxes.Tag = ['Plane #' num2str(iPlane)];
                end
            else
                %----- Plot individual plane heatmaps -----
                if combineStimTrials
                    jTab = iTab - 1;
                else
                    jTab = iTab;
                end
                stimHeatmapSubtabs{iTab}.Title = ['Plane #', num2str(jTab)];
                stimHeatmapSubtabs{iTab}.Tag = ['plane', num2str(jTab)];
                
                % Determine necessary subplot arrangement
                nSubplots = size(dffAvg, 4) + 1;
                sizeRatio = stimHeatmapSubtabs{iTab}.Position(3) / stimHeatmapSubtabs{iTab}.Position(4);
                if nSubplots == 3
                    subplotPos = [2 2]; % numSubplots() would return [1 3] for this
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
                
                % Plot reference image for the current plane in the first axis
                stimHeatmapAxes{1} = axes(stimHeatmapSubtabs{iTab}, 'Units', 'Normalized', ...
                    'Position', axesPosSort(1,:), 'Tag', 'heatmapRefImg');
                imshow(myData.refImg{jTab}, [0 maxIntensity]);
                title(['Baseline duration = ', baselineDurBox.String, ',  Response duration = ', respDurBox.String])
                stimHeatmapAxes{1}.Tag = 'heatmapRefImg';
                
                % Plot dF/F heatmaps for each of the trial types in the other axes
                for iStim = 1:nSubplots-1
                    stimHeatmapAxes{iStim+1} = axes(stimHeatmapSubtabs{iTab}, 'Units', 'Normalized',...
                        'Position', axesPosSort(iStim+1,:));
                    imagesc(dffAvg(:,:,jTab,iStim));
                    caxis(range);
                    colormap(stimHeatmapAxes{iStim+1}, 'bluewhitered');
                    axis image; axis off;
                    if combineStimTrials
                        plotTitles = {'Wind trials', 'Control trials'};
                        stimHeatmapAxes{iStim+1} = plotTitles{iStim};
                    else
                        stimHeatmapAxes{iStim+1}.Title.String = stimTypes{iStim};
                    end
                    stimHeatmapAxes{iStim+1}.Tag = 'stimHeatmap';
                end                    
            end%if
        end%for
    end%function
%---------------------------------------------------------------------------------------------------
    function plotROIstimResponsesButton_Callback(~,~)
        
        if isempty(myData.ROIs)
            errordlg('Error: no ROIs defined!', 'Error')
        else
            
            % Find parameter objects
            stimTypeButtonGroup = findobj('Tag', 'stimTypeButtonGroup');
            baselineDurBox = findobj('Tag', 'baselineDurBox');
            respDurBox = findobj('Tag', 'respDurBox');
            
            % Delete previous ROI stim plotting tab if it exists
            if ~isempty(findobj('Tag', 'stimROITab'))
                delete(findobj('Tag', 'stimROITab'));
                drawnow();
            end
            
            % Create new tab for heatmap plots
            stimROITab = uitab(stimRespSubtabGroup, 'Title', 'ROI Stim Responses', 'Tag', 'stimROITab');
            
            % Create ROI plotting subtab group
            stimROISubtabGroup = uitabgroup(stimROITab, 'Units', 'Normalized', 'Position', [0 0 1 0.99],...
                'Tag', 'stimROISubtabGroup');
            
            % Create save button
            stimROISaveButton = uicontrol(stimROITab, 'Style', 'pushbutton', 'String', 'Save current plots', ...
                'Units', 'Normalized', 'Position', [0.89, 0.6, 0.1, 0.05], 'FontSize', 12, ...
                'Callback', {@ROISaveButton_Callback}, 'Tag', 'stimROISaveButton');
            
            disp('Creating dF/F plots for selected ROIs...')
            stimROItabs = []; dffAxes = []; smImgAxes = []; dffRefImages = []; dffROIplots = []; dffROItext = []; dffAxes = []; dffCaptures = [];
            for iROI = 1:length(myData.ROIs.plots)
                
                % Create new tab for current ROI
                stimROItabs{iROI} = uitab(stimROISubtabGroup, 'Title', ['ROI #', num2str(iROI)]);
                                
                %----- Calculate mean dF/F for baseline and response periods -----
                
                % Separate out data from the current ROI's plane
                currData = squeeze(myData.wholeSession(:,:,str2double(myData.ROIs.planes{iROI}),:,:)); % --> [x, y, volume, trial]
                
                % Determine which volumes are in the baseline and response periods
                baselineDurSec = str2num(baselineDurBox.String);
                baselineStart = floor(stimStart*volumeRate)-floor(baselineDurSec*volumeRate);
                baselineEnd = floor(stimStart*volumeRate);
                respDur = str2num(respDurBox.String);
                respStart = ceil(stimStart*volumeRate);
                respEnd = floor((stimStart+respDur)*volumeRate);
                
                combineStimTrials = ~strcmp(stimTypeButtonGroup.SelectedObject.Tag, 'trialTypeRadio');
                wholeTrialAvg = []; baselineAvg = [];
                if ~combineStimTrials
                    % Calculate dF/F for each trial type separately
                    for iStim = 1:length(stimTypes)
                        % Get trial averaged data
                        wholeTrialAvg(:,:,:,iStim) = mean(currData(:,:,:,myData.stimSepTrials.(stimTypes{iStim})), 4);  % --> [x, y, volume, StimType]
                        baselineAvg(:,:,iStim) = mean(wholeTrialAvg(:,:, baselineStart:baselineEnd, iStim), 3);         % --> [x, y, StimType]
                    end%for
                else
                    % Get trial averaged data for all combined stim trials vs. control trials
                    wholeTrialAvg(:,:,:,1) = squeeze(mean(currData(:,:,:,myData.stimSepTrials.windTrials), 4));         % --> [x, y, volume, StimType]
                    baselineAvg(:,:,1) = squeeze(mean(wholeTrialAvg(:,:, baselineStart:baselineEnd, 1), 3));            % --> [x, y, StimType]
                    wholeTrialAvg(:,:,:,2) = squeeze(mean(currData(:,:,:,~myData.stimSepTrials.windTrials), 4));        % --> [x, y, volume, StimType]
                    baselineAvg(:,:,2) = squeeze(mean(wholeTrialAvg(:,:, baselineStart:baselineEnd, 2), 3));            % --> [x, y, StimType]              
                end%if
                
                trialAvg = wholeTrialAvg(:,:,baselineStart:respEnd,:);                                                  % --> [x, y, volume, stimType]
                
                % Zero baseline and trial averaged data outside of ROI
                currVolNum = size(trialAvg, 3);
                repMask = repmat(myData.ROIs.masks(:,:,iROI), [1 1 currVolNum]);                       % Expand ROI mask to cover all volumes --> [x, y, volume]
                ROIdata = []; baselineAvgROI = [];
                for iStim = 1:size(wholeTrialAvg, 4)
                    baselineAvgROI(:,:,iStim) = baselineAvg(:,:,iStim) .* myData.ROIs.masks(:,:,iROI); % --> [x, y, stimType]
                    baselineRepROI = permute(repmat(baselineAvgROI, [1 1 1, currVolNum]), [1 2 4 3]);  % --> [x, y, volume, stimType] (permute makes sure the dimensions stay consistent)
                    ROIdata(:,:,:,iStim) = trialAvg(:,:,:,iStim) .* repMask;                           % --> [x, y, volume, stimType]
                end
                
                % Calculate mean dF/F within ROIs
                dffRaw = (ROIdata - baselineRepROI) ./ baselineRepROI;                                 % --> [x, y, volume, stimType]
                dffMean = squeeze(mean(mean(dffRaw, 1, 'omitnan'), 2, 'omitnan'));                     % --> [volume, stimType]
                
                % Offset all mean dF/F values so the average value from the baseline period is zero
                preStimMean = mean(dffMean(baselineStart:baselineEnd, :), 1);                          % --> [stimType]
                dffOffset = dffMean - repmat(preStimMean, [currVolNum, 1]);                            % --> [volume, stimType]
 
                % Calculate min and max dF/F values across all stim types for setting yLims.
                dffTrim = dffOffset(2:end-1, :); % Exclude the first and last volumes
                yL = calc_ylims(dffTrim, []);                
                
                
                % ----- Make refImg and dF/F plots -----
                
                % Plot miniature reference image with current ROI
                smImgAxes{iROI} = axes(stimROItabs{iROI}, 'Position', [0 0.59 0.4 0.4]);
                axis image; hold on
                dffRefImages{iROI} = imshow(refImg{str2double(myData.ROIs.planes{iROI})}, [0 maxIntensity], 'InitialMagnification', 'fit', 'Parent', smImgAxes{iROI});
                dffROIplots{iROI} = plot(myData.ROIs.plotData{iROI}(:,1), myData.ROIs.plotData{iROI}(:,2),'linewidth', 2, 'color', myData.ROIs.color{iROI});
                currText = myData.ROIs.nums{iROI};
                dffROItext{iROI} = text(currText.Position(1), currText.Position(2), num2str(iROI), 'Color', myData.ROIs.color{iROI});
                
                % Plot mean ROI dF/F for each stim type
                nPlots = size(dffOffset, 2);
                dffTime = (1:nVolumes) ./ volumeRate;
                dffAxes = [];
                plotHeight = (0.9 - (0.05 * nPlots)) / nPlots;
                plotWidth = 1 - 0.45 - 0.15; % To stay out of the way of the refImg and the save button
                for iPlot = 1:nPlots
                    
                    % Calculate plot location
                    firstPlotPos = [0.45, 0.06, plotWidth, plotHeight];
                    if iPlot == 0
                        plotPos = firstPlotPos;
                    else
                        plotPos = firstPlotPos + [0, (iPlot - 1) * (0.05 + plotHeight), 0, 0];
                    end
                    
                    % Plot data for current stim type
                    if combineStimTrials
                        stimTypeNames = {'WindTrials', 'ControlTrials'};
                    else
                        stimTypeNames = stimTypes;
                    end
                    dffAxes{iROI}.(stimTypeNames{iPlot}) = axes(stimROItabs{iROI}, 'Position', plotPos);
                    hold on
                    plot(dffTime(baselineStart:respEnd), dffOffset(:, iPlot));
                    title(stimTypeNames{iPlot},'FontSize', 10, 'FontWeight', 'normal')
                    ylabel('dF/F');
                    ylim(yL);
                    if iPlot == 1
                        xlabel('Time (sec)', 'FontSize', 12);
                    else
                        dffAxes{iROI}.XTickLabel = '';
                        set(gca, 'XTickLabel', '');
                    end
                    
                    % Add shading during stimulus presentation
                    if ~combineStimTrials || strcmp(stimTypeNames{iPlot}, 'WindTrials')
                        xL = xlim();
                        rectPos = [stimDuration(1), yL(1), stimDuration(2), diff(yL)];
                        rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');
                        ylim(yL);
                        xlim(xL);
                    end                    
                end%for
               
                myData.dffData(iROI).dffRaw = dffRaw;
                myData.dffData(iROI).dffOffset = dffOffset;               
            end%for
            disp('Plotting complete')
        end%if
    end%function
%---------------------------------------------------------------------------------------------------
    function drawROIButton_Callback(~, ~)
        
        % ---------- Draw ROI and save relevant information about it ----------
        cm = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Cyan'); rgb('Purple'); rgb('Brown'); ...
            rgb('Indigo'); rgb('DarkRed') ; rgb('Magenta') ; rgb('Gold')];
        currcolor = cm(mod(index,size(cm,1))+1,:); % index starts at 1 when gui initializes
        
        % Prompt user to create a polygon ROI
        [myData.ROIs.masks(:,:,index), xi, yi] = roipoly; % --> [x, y, ROInum]
        currAxes = gca;
        
        % Save other useful information about the ROI
        numLoc = strfind(currAxes.Tag, '#');
        myData.ROIs.planes{index} = currAxes.Tag(numLoc+1:end);
        myData.ROIs.plots{index} = plot(xi, yi, 'linewidth', 2, 'color', currcolor);
        myData.ROIs.plotData{index} = [xi, yi];
        myData.ROIs.nums{index} = text(mean(xi),mean(yi),num2str(index),'Color',currcolor, ...
            'FontSize',12);
        myData.ROIs.color{index} = currcolor;
        
        index = index + 1; % Track total # of ROIs that have been drawn
    end
%---------------------------------------------------------------------------------------------------
    function clearROIButton_Callback(~, ~)
        % Clear all existing ROIs and plots
        for iROI = 1:length(myData.ROIs.plots)
            delete(myData.ROIs.plots{iROI})
            delete(myData.ROIs.nums{iROI})
        end
        myData.ROIs = [];
        index = 1; % Reset global count of total # of ROIs drawn
        drawnow()
        disp('ROIs cleared')
    end
%---------------------------------------------------------------------------------------------------
    function behavSummarySaveButton_Callback(~,~)
        
        % Prompt user to specify a save directory
        saveDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select a save directory');
        if saveDir == 0
            disp('Saving cancelled')
        end
        
        baseFileName = 'Behavior_Summary';
        
        % Let user modify base file name if desired
        fileName = inputdlg('Enter base file name:', 'Save behavior summary plots as', ...
            1, {baseFileName});
        if isempty(fileName)
            disp('Save cancelled')
        else
            fileName = fileName{:};
            
            % Warn user and offer to cancel save if this will overwrite existing files
            overwrite = 1;
            if exist(fullfile(saveDir, [fileName, '_2D.png']), 'file') ~= 0 || ...
                    exist(fullfile(saveDir, [fileName, '_1D_wind.png']), 'file') ~= 0 || ...
                    exist(fullfile(saveDir, [fileName, '_1D_noWind.png']), 'file') ~= 0
                dlgAns = questdlg('Saving this plot will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
                if strcmp(dlgAns, 'No')
                    overwrite = 0;
                    disp('Saving cancelled')
                end
            end
            
            % Save plots
            if overwrite
                
                % Have to copy plots over into a new figure to save just a single axes
                tempFig = figure('Visible', 'Off', 'Color', [1 1 1]);
                tempFig.Units = 'normalized';
                tempFig.OuterPosition = [0 0 0.85 0.85];
                copyobj(plotAxes2D, tempFig);
                export_fig(fullfile(saveDir, [fileName, '_2D.png']), '-png');
                clf; copyobj(plotAxes1D_wind, tempFig);
                export_fig(fullfile(saveDir, [fileName, '_1D_wind.png']), '-png');
                clf; copyobj(plotAxes1D_noWind, tempFig);
                export_fig(fullfile(saveDir, [fileName, '_1D_noWind.png']), '-png');
                close(tempFig);
                disp('Plots saved')
            end
        end%if
    end%function
%---------------------------------------------------------------------------------------------------
    function heatmapSaveButton_Callback(~,~)
        
        % Find currently active stim heatmap tab
        base = findobj('Tag', 'baseTabGroup');
        L1 = base.SelectedTab;
        G1 = findobj(L1, 'Type', 'uitabgroup', '-depth', 1); % e.g. stimRespSubtabGroup
        L2 = G1.SelectedTab;
        G2 = findobj(L2, 'Type', 'uitabgroup', '-depth', 1); % e.g. stimHeatmapSubtabGroup
        currTab = G2.SelectedTab;
        planeNum = currTab.Tag(end);
        
        % Identify plots to be saved
        currTabHeatmaps = findobj(currTab.Children, 'Tag', 'stimHeatmap');
        
        
        % Get baseline and response period durations from refImg title
        planeRefImg = findobj(currTab.Children, 'Tag', 'heatmapRefImg');
        titleStr = planeRefImg.Title.String;
        baselineDur = titleStr(strfind(titleStr, ',')-1);
        titleStrShort = titleStr(strfind(titleStr, ','):end);
        respDur = titleStrShort(strfind(titleStrShort, '=') + 2);
        
        % Prompt user to specify a save directory
        saveDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select a save directory');
        if saveDir == 0
            disp('Saving cancelled')
        else
            
            % Generate base file name
            baseFileName = ['dFF_Heatmap_Plane_', planeNum, '_Base_', baselineDur, '_Resp_', respDur];
            
            % Let user modify base file name if desired
            fileName = inputdlg('Enter base file name:', 'Save dF/F heatmap plots as', ...
                [1 50], {baseFileName});
            if isempty(fileName)
                disp('Save cancelled')
            else
                fileName = fileName{:};
                
                % Warn user and offer to cancel save if this will overwrite existing files
                overwrite = 1;
                for iAxes = 1:length(currTabHeatmaps)
                    if exist(fullfile(saveDir, [fileName, '_', currTabHeatmaps(iAxes).Title.String, '.png']), 'file') ~= 0
                        dlgAns = questdlg('Saving this plot will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
                        if strcmp(dlgAns, 'No')
                            overwrite = 0;
                            disp('Saving cancelled')
                        end
                        break
                    end
                end
                
                % Save plots as .png files
                if overwrite
                    
                    % Have to copy plots over into a new figure to save just a single axes
                    tempFig = figure('Visible', 'Off', 'Color', [1 1 1]);
                    tempFig.Units = 'normalized';
                    tempFig.OuterPosition = [0 0 0.85 0.85];
                    for iAxes = 1:length(currTabHeatmaps)
                        copyobj(currTabHeatmaps(iAxes), tempFig);
                        export_fig(fullfile(saveDir, [fileName, '_', currTabHeatmaps(iAxes).Title.String]), '-png');
                        clf;
                    end
                    close(tempFig);
                    disp('Plots saved')
                end
            end%if
        end%if
    end%function

end%function


function myData = load_imaging_data()
% Prompts user for input data file(s) and loads data

% Prompt user for imaging data file
[dataFile, pathName, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
if dataFile == 0
    disp('Initialization cancelled')
    myData = []; % Skip loading if user clicked "Cancel"
else
    disp(['Loading ' dataFile, '...'])
    imgData = load([pathName, dataFile]); % Fields are: 'regProduct','trialType','origFileNames','tE_sec', 'expDate'
    disp([dataFile, ' loaded'])
    
    % Prompt user for behavioral annotation data file
    [annotDataFile, pathName, ~] = uigetfile('*.mat', 'Select a behavioral annotation data file if desired', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
    if annotDataFile == 0
        disp('No behavioral annotation data selected')
        annotData.behaviorLabels = [];
        annotData.goodTrials = [];
        annotData.trialAnnotations = [];
    else
        disp(['Loading ' annotDataFile, '...'])
        annotData = load([pathName, annotDataFile]);
        disp([annotDataFile, ' loaded'])
    end
    
    myData = setstructfields(imgData, annotData); % Combine imaging data and annotation data into one structure
    
    
    % Load .mat file containing trial metadata
    
    
    % Process raw data structure
    if ~isfield(myData, 'wholeSession')
        myData.wholeSession = myData.regProduct; % --> [x, y, plane, volume, trial]
    end
    myData.nTrials = size(myData.wholeSession, 5);
    singleTrial = squeeze(myData.wholeSession(:,:,:,:,1));  % --> [x, y, plane, volume]
    myData.nPlanes = size(singleTrial, 3);
    myData.nVolumes = size(singleTrial, 4);
    myData.stimTypes = sort(unique(myData.trialType));
    if ~isempty(annotData.trialAnnotations)
        myData.nFrames = max(cellfun(@height, myData.trialAnnotations));
    else
        myData.nFrames = [];
    end
    
    % For compatibility with early experiments
    if ~isfield(myData, 'expDate')
        cellDate = inputdlg('Please enter experiment date in YYYY_MM_DD format');
        myData.expDate = cellDate{1};
    end
    
    % Extract session number
    origFileName = myData.origFileNames{1};
    sidLoc = strfind(origFileName, 'sid_');
    myData.sid = origFileName(sidLoc+4);
    
    % Create mean reference image for each plane
    myData.refImg = [];
    for iPlane = 1:myData.nPlanes
        myData.refImg{iPlane} = squeeze(mean(mean(myData.wholeSession(:,:,iPlane,:,:),4),5)); % --> [x, y]
    end
    
    % Separate out wind trials
    myData.stimSepTrials = [];
    for iStim = 1:length(myData.stimTypes)
        myData.stimSepTrials.(myData.stimTypes{iStim}) = logical(cellfun(@(x) ...
            strcmp(x, myData.stimTypes{iStim}), myData.trialType));
    end
    myData.stimSepTrials.windTrials = logical(myData.stimSepTrials.LeftWind);
    
end%if
end%function




