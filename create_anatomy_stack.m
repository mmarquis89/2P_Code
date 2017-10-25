function create_anatomy_stack(dirPath, fileStr, outputFilePrefix)
%===================================================================================================
% MAKE AVERAGED ANATOMY STACK FROM INDIVIDUAL 2P VOLUMES
% Will average together all stacks in the directory that meet the requirement specified by 'fileStr'
% and save as .tif files the mean stack as well as a max Z-projection. 
% 
% INPUTS:
%    dirPath = folder with stack files, e.g. 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
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
        firstStack = uint16(read_tif(fullfile(dirPath,stacks(iVol).name)));                  % --> [x, y, plane, volume, channel]
        summedStacks = uint16(firstStack);                                                   % --> [x, y, plane, volume, channel]
    else
        summedStacks = summedStacks + uint16(read_tif(fullfile(dirPath,stacks(iVol).name))); % --> [x, y, plane, volume, channel]
    end
end
% disp(['Max summed value = ', num2str(max(summedStacks(:)))]);

% Convert to double
summedStacks = squeeze(double(summedStacks));                                                % --> [x, y, plane, channel]                                            

% Calculate mean by dividing by total number of stacks and scale from 0-1
avgStack = summedStacks ./ nStacks;                                                          % --> [x, y, plane, channel]
avgStackScaled = avgStack ./ max(avgStack(:));                                               % --> [x, y, plane, channel]



% Get max Z-projection of averaged stack and scale intensity range to 0-1
maxZ = max(avgStack, [], 3);                                                                 % --> [x, y, channel]
maxZScaled = maxZ ./ max(maxZ(:));                                                           % --> [x, y, channel]

% Make sure files with these names don't already exist in the chosen directory
assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');
assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'file')==0, 'Error: a file with this name already exists in this directory');

% Check whether data has multiple channels
nChannels = size(summedStacks, 4);

if nChannels == 1
    
    % Write averaged stack and Z-projection to .tif files
    imwrite(maxZScaled, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']));
    imwrite(avgStackScaled(:,:,1), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']));
    for iPlane = 2:size(avgStackScaled, 3)
        imwrite(avgStackScaled(:,:,iPlane), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'writemode','append');
    end
    
else % nChannels == 2
    
    % Separate data by channel if necessary
    avgStackScaled_1 = avgStackScaled(:,:,:,1);
    maxZScaled_1 = maxZScaled(:,:,1);
    if nChannels == 2
        avgStackScaled_2 = avgStackScaled(:,:,:,2);
        maxZScaled_2 = maxZScaled(:,:,2);
    end
    
    % Write Z-projection to .tif files
    imwrite(maxZScaled, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_1.tif']));
    imwrite(maxZScaled, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ_2.tif']));
    
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




