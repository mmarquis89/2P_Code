function create_anatomy_stack(dirPath, fileStr, outputFilePrefix)
%===================================================================================================
% MAKE AVERAGED ANATOMY STACK FROM INDIVIDUAL 2P VOLUMES
% Will average together all stacks in the directory that meet the requirement specified by 'fileStr'
% and save as .tif files the mean stack as well as a max Z-projection. 
% 
% INPUTS:
%    dirPath = folder with stack files, e.g. 'B:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
%
%    fileStr = filter string for dir() selection to get anatomy stack files, e.g. '*stack_*'. Be 
%              careful to ensure that only the desired files will meet this specification.
%
%    outputFilePrefix = name to prepend to the names of output files. Can be [].
%===================================================================================================

stacks = dir(fullfile(dirPath, fileStr));
nStacks = length(stacks);

% Add underscore to prefix if one is provided
if ~isempty(outputFilePrefix)
   outputFilePrefix = [outputFilePrefix, '_']; 
end

% Throw error if the files about to be created already exist
assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']), 'file')==0 ... 
    && exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'file')==0 ...
    , 'Error: output .tif files already exist in this directory');

% Get sum of all stacks
summedStacks = [];
for iVol = 1:nStacks
    disp(['Reading ' stacks(iVol).name, '...']);
    if iVol == 1
        firstStack = uint16(read_tif(fullfile(dirPath,stacks(iVol).name)));                  % --> [y, x, plane, volume, channel]
        summedStacks = uint16(firstStack);                                                   % --> [y, x, plane, volume, channel]
    else
        summedStacks = summedStacks + uint16(read_tif(fullfile(dirPath,stacks(iVol).name))); % --> [y, x, plane, volume, channel]
    end
end
% disp(['Max summed value = ', num2str(max(summedStacks(:)))]);

% Convert to double
summedStacks = squeeze(double(summedStacks));                                                % --> [y, x, plane, channel]                                            

% Check whether data has multiple channels
nChannels = size(summedStacks, 4);

if nChannels == 1
    
    % Calculate mean by dividing by total number of stacks and scale from 0-1
    avgStack = summedStacks ./ nStacks;             % --> [y, x, plane, channel]
    avgStackScaled = avgStack ./ max(avgStack(:));  % --> [y, x, plane, channel]
    
    % Get max Z-projection of averaged stack and scale intensity range to 0-1
    maxZ = squeeze(max(avgStack, [], 3));           % --> [y, x, channel]
    maxZScaled = maxZ ./ max(maxZ(:));              % --> [y, x, channel]
    
    % Make sure files with these names don't already exist in the chosen directory
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    
    % Write averaged stack and Z-projection to .tif files
    imwrite(maxZScaled, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']));
    imwrite(avgStackScaled(:,:,1), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']));
    for iPlane = 2:size(avgStackScaled, 3)
        imwrite(avgStackScaled(:,:,iPlane), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'writemode','append');
    end
    
else % nChannels == 2
    
    % Calculate mean by dividing by total number of stacks
    avgStack = summedStacks ./ nStacks;                  % --> [y, x, plane, channel]
    
    % Separate data by channel and scale intensity range to 0-1
    avgStack_1 = avgStack(:,:,:,1);                      % --> [y, x, plane]
    avgStack_2 = avgStack(:,:,:,2);                      % --> [y, x, plane]
    avgStackScaled_1 = avgStack_1 ./ max(avgStack_1(:)); % --> [y, x, plane]
    avgStackScaled_2 = avgStack_2 ./ max(avgStack_2(:)); % --> [y, x, plane]
    
    % Get max Z-projection of averaged stack and scale intensity range to 0-1
    maxZ_1 = squeeze(max(avgStack_1, [], 3));            % --> [y, x]
    maxZ_2 = squeeze(max(avgStack_2, [], 3));            % --> [y, x]
    maxZScaled_1 = maxZ_1 ./ max(maxZ_1(:));             % --> [y, x]
    maxZScaled_2 = maxZ_2 ./ max(maxZ_2(:));             % --> [y, x]
    
    % Make sure files with these names don't already exist in the chosen directory
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_1.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_2.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack_1.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack_2.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
    
    % Write Z-projection to .tif files
    imwrite(maxZScaled_1, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_1.tif']));
    imwrite(maxZScaled_2, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_2.tif']));
    
    % Write averaged stack to .tif files
    imwrite(avgStackScaled_1(:,:,1), fullfile(dirPath, [outputFilePrefix, 'MeanStack_1.tif']));
    for iPlane = 2:size(avgStackScaled_1, 3)
        imwrite(avgStackScaled_1(:,:,iPlane), fullfile(dirPath, [outputFilePrefix, 'MeanStack_1.tif']), 'writemode','append');
    end
    imwrite(avgStackScaled_2(:,:,1), fullfile(dirPath, [outputFilePrefix, 'MeanStack_2.tif']));
    for iPlane = 2:size(avgStackScaled_2, 3)
        imwrite(avgStackScaled_2(:,:,iPlane), fullfile(dirPath, [outputFilePrefix, 'MeanStack_2.tif']), 'writemode','append');
    end
    
end

