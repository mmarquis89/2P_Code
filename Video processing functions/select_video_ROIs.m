function select_video_ROIs(parentDir, sid)







% Load and plot 100th frame acquired during the first trial (number choice is arbitrary)
myFolders = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*tid*']));
myFolders = myFolders([myFolders.isdir]);
firstFolder = myFolders(1);
myFrames = dir(fullfile(parentDir, firstFolder.name, '*.tif'));
plotImg = uint8(imread(fullfile(parentDir, firstFolder.name, myFrames(1).name)));

h = figure(1);clf; hold on
im = imshow(plotImg, []);

% Prompt user to create an ROI
h.Name = 'Define an ROI for fly movement';
disp('First, define an ROI for fly movement')
roiData = []; xi = []; yi = [];
[roiData(:,:,1), xi{1}, yi{1}] = roipoly; % --> [y, x, ROInum]
h.Name = 'Define an ROI for washer movement';
disp('Next, define an ROI for washer movement')
[roiData(:,:,2), xi{2}, yi{2}] = roipoly; % --> [y, x, ROInum]

% Save ROI data
saveDir = uigetdir(parentDir, 'Select a save directory');
if saveDir == 0
    % Throw error if user canceled without choosing a directory
    disp('ERROR: you must select a save directory or provide one as an argument');
else
    % Prompt user for file name
    fileName = inputdlg('Please choose a file name', 'Save ROI data', 1, {'Behavior_Vid_ROI_Data'});
    fileName = fileName{:};
    
    % Save data
    save(fullfile(saveDir, [fileName, '.mat']), 'roiData');
end
disp('ROI data saved')

end