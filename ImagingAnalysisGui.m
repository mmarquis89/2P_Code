function ImagingAnalysisGui()
%
%
%
%
%
%
%
% myData.sessionData = [x, y, plane, volume, trial]

close all

%----------INITIALIZATION TASKS----------

% Prompt user for data file and load data
myData = loadData();
if ~isempty(myData) % Abort initialization if no data was loaded
    nPlanes = myData.nPlanes;
    nVolumes = myData.nVolumes;
    refImg = myData.refImg;
    
    % Create global variables
    myData.maxIntensity = 700; maxIntensity = myData.maxIntensity;
    myData.volumeRate = 6.5; volumeRate = myData.volumeRate;
    myData.stimDuration = [4 2]; stimDuration = myData.stimDuration; % [start time, length] in seconds
    myData.ROIs= [];
    myData.dffData = [];
    index = 1;
    planeAxes = [];
    dffPlotCounter = 1;
    
    %----------CONSTRUCTING COMPONENTS----------
    
    % Figure
    f = figure('position', [50 45, 1800, 950], 'Tag', 'figure',...
        'WindowButtonDownFcn', {@figure_WindowButtonDownFcn}, 'Name', 'Imaging Analysis GUI');
    
    % Image tabs
    tabGroup = uitabgroup(f, 'Units', 'Pixels', 'Position', [10 10 1570 920],...
        'CreateFcn', {@tabGroup_CreateFcn}, 'Tag', 'tabGroup');
    
    % DrawROI button
    drawROIButton = uicontrol('Style', 'pushbutton', 'String', 'Draw ROI',...
        'Position', [f.Position(3)-200, f.Position(4)-100, 180, 50],...
        'FontSize', 12, 'Callback', {@drawROIButton_Callback}, 'Tag', 'drawROIButton');
    
    % Create clearROI button
    clearROIButton = uicontrol('Style', 'pushbutton', 'String', 'Clear ROIs',...
        'Position', [f.Position(3)-200, f.Position(4)-200, 180, 50],...
        'FontSize', 12, 'Callback', {@clearROIButton_Callback}, 'Tag', 'clearROIButton');
    
    % Create "Plot dF/F" button
    plotDffButton = uicontrol('Style', 'pushbutton', 'String', 'Plot dF/F',...
        'Position', [f.Position(3)-200, f.Position(4)-300, 180, 50],...
        'FontSize', 12, 'Callback', {@plotDffButton_Callback}, 'Tag', 'plotDffButton');
    
    % Create "Save dF/F" button
    saveDffButton = uicontrol('Style', 'pushbutton', 'String', 'Save dF/F data',...
        'Position', [f.Position(3)-200, f.Position(4)-400, 180, 50],...
        'FontSize', 12, 'Callback', {@saveDffButton_Callback}, 'Tag', 'saveDffButton');
    
% Should have a "load ROIs" button as well?
    
    % Change "Units" to normalized so components resize automatically.
    f.Units = 'normalized';
    tabGroup.Units = 'normalized';
    drawROIButton.Units = 'normalized';
    clearROIButton.Units = 'normalized';
    plotDffButton.Units = 'normalized';
    for iTab = 1:length(myTabs)
        myTabs{iTab}.Units = 'normalized';
    end
    for iAxes = 1:length(myAxes)
        myAxes{iAxes}.Units = 'normalized';
    end
    
end%if

%------------- CALLBACK FUNCTIONS ------------------------------------------------------------------


%---------------------------------------------------------------------------------------------------
    function tabGroup_CreateFcn(src, ~)
        
        %Create image tabs for each plane and for combined images
        myTabs = [];
        singleImages = [];
        for iTab = 1:nPlanes+1
            myTabs{iTab} = uitab(src);
            if iTab == 1
                
                %-----Plot combined image of all planes-----
                myTabs{1}.Title = 'All planes';
                
                % Calculate necessary axis dimensions
                subplotDim1 = ceil(sqrt(nPlanes));
                subplotDim2 = floor(sqrt(nPlanes));
                axWidth = floor((src.Position(3)- 40 - 10*(subplotDim1-1)) / subplotDim1);
                axHeight = floor((src.Position(4) - 40 - 10*(subplotDim2-1)) / subplotDim2);
                
                % Place axes in correct locations
                myAxes = [];
                for yDim = 1:subplotDim2
                    for xDim = 1:subplotDim1
                        myAxes{yDim, xDim} = axes(myTabs{1}, 'Units', 'Pixels', 'Position', [xDim*10+(xDim-1)*axWidth, yDim*10+(yDim-1)*axHeight, axWidth, axHeight]);
                    end
                end
                myAxes = reshape(flipud(myAxes)', 1, nPlanes); % Now handles are ordered by z-planes, displayed from L-R, top-bottom
                
                % Plot reference frames
                myImages = [];
                for iPlane = 1:nPlanes
                    currAxes = myAxes{iPlane};
                    axes(currAxes);
                    axis image; hold on
                    myImages{iPlane} = imshow(refImg{iPlane}, [0 maxIntensity], 'InitialMagnification', 'fit', 'Parent', currAxes);
                    myImages{iPlane}.ButtonDownFcn = {@image_ButtonDownFcn};
                    currAxes.Title.String = ['Plane  #' num2str(iPlane)];
                    currAxes.Tag = ['Plane  #' num2str(iPlane)];
                end
            else
                %-----Plot individual plane reference images-----
                myTabs{iTab}.Title = ['Plane #', num2str(iTab-1)];
                planeAxes{iTab-1} = axes(myTabs{iTab}, 'Position', [0.1 0.1 0.8 0.8]);
                axis image; hold on
                singleImages{iTab-1} = imshow(refImg{iTab-1}, [0 maxIntensity], 'InitialMagnification', 'fit', 'Parent', planeAxes{iTab-1});
                singleImages{iTab-1}.ButtonDownFcn = {@image_ButtonDownFcn};
                planeAxes{iTab-1}.Title.String = ['Plane  #' num2str(iTab-1)];
                planeAxes{iTab-1}.Tag = ['Plane  #' num2str(iTab-1)];
            end
        end        
    end

%---------------------------------------------------------------------------------------------------
    function figure_WindowButtonDownFcn(src, ~)
        % Resets all the axes titles on the "All planes" tab whenever the window is clicked
        allAxes = src.Children(end).Children(1).Children;
        for iAx = 1:length(allAxes)
            allAxes(iAx).Title.String = allAxes(iAx).Tag;
            planeAxes{iAx}.Title.String = planeAxes{iAx}.Tag;
        end
    end

%---------------------------------------------------------------------------------------------------
    function image_ButtonDownFcn(src, ~)
        % Append [SELECTED] to the title of a clicked image on the "All planes" tab
        src.Parent.Title.String = [src.Parent.Title.String(end-2:end), ' [SELECTED]'];
    end

%---------------------------------------------------------------------------------------------------
    function drawROIButton_Callback(~, ~)
        
        % Draw ROI and save relevant information about it
        cm = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Cyan'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ; rgb('Magenta') ; rgb('Gold')];
        currcolor = cm(mod(index,size(cm,1))+1,:);
        
        [myData.ROIs.masks{index}, xi, yi] = roipoly;
        currAxes = gca;
        
        numLoc = strfind(currAxes.Tag, '#');
        myData.ROIs.planes{index} = currAxes.Tag(numLoc+1:end);
        myData.ROIs.plots{index} = plot(xi,yi, 'linewidth', 2, 'color', currcolor);
        myData.ROIs.plotData{index} = [xi, yi]; 
        myData.ROIs.nums{index} = text(mean(xi),mean(yi),num2str(index),'Color',currcolor,'FontSize',12);
        myData.ROIs.color{index} = currcolor;
        
        index = index + 1; % Track total # of ROIs drawn
    end

%---------------------------------------------------------------------------------------------------
    function clearROIButton_Callback(~, ~)
        % Clear all existing ROIs and plots
        for iROI = 1:length(myData.ROIs.plots)
           delete(myData.ROIs.plots{iROI}) 
           delete(myData.ROIs.nums{iROI})
        end
        myData.ROIs = []; 
       index = 1; % Reset ROI count for drawing new ROIs
       drawnow()
       disp('ROIs cleared')
    end

%---------------------------------------------------------------------------------------------------
    function saveDffButton_Callback(~, ~)
        % Saves the current raw dF/F data for each ROI in a .mat file
        disp('Saving dF/F data for selected ROIs...')
        
        % Get save directory from user
        saveDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\');
        if saveDir == 0
            disp('Save cancelled')
        else            
            % Extract session number
            origFileName = myData.origFileNames{1};
            sidLoc = strfind(origFileName, 'sid_');
            sessionName = origFileName(sidLoc:sidLoc+4);
            baseFileName = ['dff_data_', sessionName];
            
            % Let user modify file name if desired
            fileName = inputdlg('Enter file name (without .mat extension):', 'Save dF/F data as', 1, {baseFileName});
            if isempty(fileName)
                disp('Save cancelled')
            else
                % Save file
                save(fullfile(saveDir,fileName{:}), '-struct', 'myData', 'dffData', 'ROIs');
                disp('Saving complete')
            end
        end
    end

%---------------------------------------------------------------------------------------------------
    function plotDffButton_Callback(~, ~)
    % Create a new tab group containing a subtab for each ROI, then plot the reference image and
    % ROI along with the mean dF/F on the page for each stimulus and each ROI.
        
        if isempty(myData.ROIs)
            disp('Error: no ROIs defined!')
        else
            disp('Creating dF/F plots for selected ROIs...')
            myData.dffData = [];
            % Create one more tab for plotting data
            plottingTab = uitab(tabGroup);
            plottingTab.Tag = 'plottingTab';
            if dffPlotCounter > 1
                plottingTab.Title = ['dF/F Plotting #', num2str(dffPlotCounter)];
            else
                plottingTab.Title = 'dF/F Plotting';
            end
            
            % Create a subgroup of tabs for this page
            subTabGroup = uitabgroup(plottingTab, 'Units', 'Pixels', 'Position', [5 1 1555 880], 'Tag', 'subTabGroup');
            ROItabs = []; dffAxes = []; smImgAxes = []; dffRefImages = []; dffROIplots = []; dffROItext = []; dffAxes = []; dffCaptures = [];
            for iROI = 1:length(myData.ROIs.plots)
                
                ROItabs{iROI} = uitab(subTabGroup, 'Title', ['ROI #', num2str(iROI)]);

                % Plot miniature reference image with current ROI
                smImgAxes{iROI} = axes(ROItabs{iROI}, 'Position', [0 0.59 0.4 0.4]);
                axis image; hold on
                dffRefImages{iROI} = imshow(refImg{str2double(myData.ROIs.planes{iROI})}, [0 maxIntensity], 'InitialMagnification', 'fit', 'Parent', smImgAxes{iROI});
                dffROIplots{iROI} = plot(myData.ROIs.plotData{iROI}(:,1), myData.ROIs.plotData{iROI}(:,2),'linewidth', 2, 'color', myData.ROIs.color{iROI});
                currText = myData.ROIs.nums{iROI};
                dffROItext{iROI} = text(currText.Position(1), currText.Position(2), num2str(iROI), 'Color', myData.ROIs.color{iROI});
                
                
                %-----Plot mean dF/F across all trials for the current ROI-----
                stimTypes = sort(unique(myData.trialType)); myData.dffData(iROI).stimTypes = stimTypes;
                currData = squeeze(myData.wholeSession(:,:,str2double(myData.ROIs.planes{iROI}),:,:)); % --> [x, y, volume, trial]

                % For each stimulus type...
                stimSepTrials = []; trialAvg = []; baselineAvg = []; ROIdata = []; baselineF = []; dff = []; dffRaw = [];
                for iStim = 1:length(stimTypes)
                    
                    % Separate trials by stimulus type
                    stimSepTrials.(stimTypes{iStim}) = logical(cellfun(@(x) strcmp(x, stimTypes{iStim}), myData.trialType));
                    
                    % Get trial averaged and baseline data
                    trialAvg{iStim} = mean(currData(:,:,:,stimSepTrials.(stimTypes{iStim})),4); % --> [x, y, volume]
                    baselineAvg{iStim} = mean(trialAvg{iStim},3); % --> [x, y] %(:,:,1:(2*volumeRate))
                    
                    % Zero baseline and trial averaged data outside of ROI
                    baselineAvg{iStim}(~myData.ROIs.masks{iROI}) = 0;
                    volumeMask = repmat(myData.ROIs.masks{iROI}, [1 1 nVolumes]); % Expand ROI mask to cover all volumes
                    ROIdata{iStim} = trialAvg{iStim};
                    ROIdata{iStim}(~volumeMask) = 0;
                    
                    % Calculate dF/F within ROIs
                    baselineF{iStim} = sum(baselineAvg{iStim}(:));
                    ROIdataF{iStim} = squeeze(sum(sum(ROIdata{iStim},1),2));
                    dffRaw{iStim} = (ROIdataF{iStim} - baselineF{iStim}) ./ baselineF{iStim};
                    
                    % Offset all mean dF/F values so that the average value from the
                    % first two seconds of the trial is at zero
                    preStimMean = mean(dffRaw{iStim}(1 : (2*volumeRate) - 1));
                    dff{iStim} = dffRaw{iStim} - preStimMean;
                    
                end%for                                
                
                % Calculate min and max dF/F values across all stim types for setting yLims. 
                dffVals = [dff{1}(2:end), dff{2}(2:end), dff{3}(2:end)]; % Dropping the first volume of each set.
                yMax = max(dffVals(:));
                yMin = min(dffVals(:));
                padFactor = 0.1 * diff([yMin, yMax]); % To give a little extra padding
                yL = [yMin - padFactor, yMax + padFactor];
                
                % Plot mean ROI dF/F for each stim type
                nPlots = length(stimTypes);
                dffTime = (1:nVolumes) ./ volumeRate;
                dffAxes = [];
                plotHeight = (0.8 - (0.05 * nPlots)) / nPlots;
                plotWidth = 0.54;
                for iPlot = 1:nPlots
                    
                    % Calculate plot location
                    firstPlotPos = [0.45, 0.06, plotWidth, plotHeight];
                    if iPlot == 0
                        plotPos = firstPlotPos;
                    else
                        plotPos = firstPlotPos + [0, (iPlot - 1) * (0.05 + plotHeight), 0, 0];
                    end
                
                    % Plot data for current stim type
                    dffAxes{iROI}.(stimTypes{iPlot}) = axes(ROItabs{iROI}, 'Position', plotPos);
                    hold on
                    plot(dffTime(2:end), dff{iPlot}(2:end)); % Dropping first volume of each trial
                    title(stimTypes{iPlot},'FontSize', 10, 'FontWeight', 'normal')
                    ylabel('dF/F');
                    ylim(yL);
                    if iPlot == 1
                        xlabel('Time (sec)', 'FontSize', 12);
                    else
                        dffAxes{iROI}.XTickLabel = '';
                        set(gca, 'XTickLabel', '');
                    end%if
                    
                    % Add shading during stimulus presentation
                    rectPos = [stimDuration(1), yL(1), stimDuration(2), diff(yL)];
                    rectangle('Position', rectPos, 'FaceColor', [rgb('red'), 0.05], 'EdgeColor', 'none');                    
                    ylim(yL);
                end%for
                 
                % Capture the contents of the current subtab as an image
                tabGroup.SelectedTab = plottingTab;
                subTabGroup.SelectedTab = ROItabs{iROI};
                rectPos = [10 10 1560 880];
                dffCaptures = frame2im(getframe(f, rectPos));               
                
                % Save important variables to "myData" structure
                myData.dffData(iROI).stimSepTrials = stimSepTrials;
                myData.dffData(iROI).baselineF = baselineF;
                myData.dffData(iROI).ROIdataF = ROIdataF;
                myData.dffData(iROI).dffRaw = dffRaw;
                myData.dffData(iROI).dff = dff;
                myData.dffData(iROI).dffCaptures = dffCaptures;
                
            end%for
            
            subtabGroup.Units = 'normalized';
            dffPlotCounter = dffPlotCounter + 1; % For naming df/f subtabs
            disp('Plotting complete');
        end%if
    end%function
end%function


function myData = loadData()
% Prompts user for input data file and loads data

    % Prompt user for data file
    [dataFile, pathName, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
    if dataFile == 0        
        disp('Initialization cancelled')
        myData = []; % Skip loading if user clicked "Cancel"
    else
        disp(['Loading ' dataFile, '...'])
        myData = load([pathName, dataFile]);
        disp([dataFile, ' loaded'])

        % Process raw data structure
        if ~isfield(myData, 'wholeSession')
            myData.wholeSession = myData.regProduct;
        end
        myData.sessionData = squeeze(myData.wholeSession(:,:,:,:,1)); %squeeze(mean(wholeSession(:,:,:,:,1), 3)); [x y plane vol trial]
        myData.nPlanes = size(myData.sessionData, 3);
        myData.nVolumes = size(myData.sessionData, 4);

        % Create reference image for each plane
        myData.refImg = [];
        for iPlane = 1:myData.nPlanes
            myData.refImg{iPlane} = squeeze(mean(mean(myData.wholeSession(:,:,iPlane,:,:),4),5));
        end
    end%if
end%function




